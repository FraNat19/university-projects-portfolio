#!/usr/bin/env python3
"""
🎯 TEXT CHUNKS - FIXED VERSION V2

CHANGES FROM ORIGINAL:
✅ Usa enriched.laws_structured PRIMA di semantic_metadata
✅ Gestisce correttamente TUTTE le citazioni (no dedup aggressiva)
✅ Upsert in Qdrant con stessi ID = AGGIORNAMENTO, non duplicati
✅ Opzione --delete-first per pulizia completa prima di reindicizzare
✅ Statistiche migliorate su citazioni legali
"""

import json
import hashlib
import re
from pathlib import Path
from typing import List, Dict
from tqdm import tqdm
import sys
import os
import torch
import argparse

sys.path.insert(0, os.path.dirname(__file__))
import importlib.util
spec = importlib.util.spec_from_file_location("config", "01_config.py")
config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(config)

from qdrant_client import QdrantClient
from qdrant_client.models import PointStruct
from sentence_transformers import SentenceTransformer

# ✅ FIX: Import corretti
from law_extraction_enhanced import EnhancedLawExtractor
from law_validation import LawValidator

# OPTIMAL VALUES
CHUNK_SIZE = 3000
OVERLAP = 200
MAX_CHUNKS_PER_DOC = 100

from doc_id_utils import clean_doc_id, extract_doc_id_from_metadata


def extract_references(text: str) -> Dict[str, List[int]]:
    """Extract figure and table references from text"""
    refs = {'figures': [], 'tables': []}
    
    # Figures
    fig_patterns = [
        r'\b[Ff]ig(?:ura|ure)?[\.\s]+(\d+)',
        r'\bFigure[\.\s]+(\d+)',
    ]
    for pattern in fig_patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        refs['figures'].extend([int(m) for m in matches])
    
    # Tables
    tab_patterns = [
        r'\b[Tt]ab(?:ella|elle|le)?[\.\s]+(\d+)',
        r'\bTable[\.\s]+(\d+)',
    ]
    for pattern in tab_patterns:
        matches = re.findall(pattern, text, re.IGNORECASE)
        refs['tables'].extend([int(m) for m in matches])
    
    refs['figures'] = sorted(set(refs['figures']))
    refs['tables'] = sorted(set(refs['tables']))
    
    return refs


def chunk_text_smart(text: str, chunk_size: int, overlap: int) -> List[str]:
    """Smart text chunking with overlap"""
    if not text or len(text) < 100:
        return []
    
    text = "\n".join([line.strip() for line in text.split("\n") if line.strip()])
    
    if len(text) < chunk_size:
        return [text]
    
    chunks = []
    start = 0
    
    while start < len(text):
        end = start + chunk_size
        
        if end >= len(text):
            chunk = text[start:].strip()
            if len(chunk) > 100:
                chunks.append(chunk)
            break
        
        split_point = text.rfind("\n\n", start, end)
        if split_point == -1 or split_point < start + chunk_size // 2:
            split_point = text.rfind(". ", start, end)
        
        if split_point == -1 or split_point < start + chunk_size // 2:
            split_point = text.rfind(" ", start, end)
        
        if split_point == -1 or split_point <= start:
            split_point = end
        
        chunk = text[start:split_point].strip()
        if len(chunk) > 100:
            chunks.append(chunk)
        
        start = max(split_point - overlap, start + 1)
    
    return chunks


def sample_chunks_intelligent(chunks: List[str], target: int) -> List[str]:
    """Intelligent sampling"""
    if len(chunks) <= target:
        return chunks
    
    n_intro = max(int(target * 0.25), 3)
    n_conclusion = max(int(target * 0.15), 2)
    n_middle = target - n_intro - n_conclusion
    
    sampled = []
    sampled.extend(chunks[:n_intro])
    
    middle_start = n_intro
    middle_end = len(chunks) - n_conclusion
    middle_chunks = chunks[middle_start:middle_end]
    
    if n_middle > 0 and middle_chunks:
        step = max(1, len(middle_chunks) // n_middle)
        sampled.extend(middle_chunks[::step][:n_middle])
    
    sampled.extend(chunks[-n_conclusion:])
    
    return sampled


def extract_hierarchical_chunks(
    json_path: Path, 
    source: str, 
    extractor: EnhancedLawExtractor,
    validator: LawValidator
) -> List[Dict]:
    """Extract chunks with VALIDATED citations - FIXED VERSION"""
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        if "document_content" not in data:
            return []
        
        web_meta = data.get("web_metadata", {})
        sem_meta = data.get("semantic_metadata", {})
        enriched = data.get("enriched", {})  # ✅ NUOVO: Leggi da enriched
        content = data["document_content"]
        
        # Extract and CLEAN doc_id
        doc_id = extract_doc_id_from_metadata(web_meta, json_path)
        
        text = content.get("markdown_content") or content.get("plain_text", "")
        
        if not text or len(text) < 100:
            return []
        
        chunk_metas = []
        
        # Basic metadata
        title = web_meta.get("title", "")
        abstract = web_meta.get("abstract", "")
        sintesi = sem_meta.get("sintesi", "")
        category = sem_meta.get("categoria_principale", "")
        keywords = sem_meta.get("parole_chiave", [])[:15]
        
        # =====================================
        # ✅ FIX CRITICO: Priorità a enriched.laws_structured
        # =====================================
        doc_level_laws_cleaned = []
        
        # 1. Prima prova enriched.laws_structured (dal nuovo reenrich_laws_v3)
        if enriched.get("laws_structured"):
            doc_level_laws_cleaned = enriched["laws_structured"]
        
        # 2. Fallback a semantic_metadata.laws_structured
        elif sem_meta.get("laws_structured"):
            doc_level_laws_cleaned = sem_meta["laws_structured"]
        
        # 3. Ultimo fallback: estrai al volo
        if not doc_level_laws_cleaned:
            laws_raw, conf = extractor.extract_with_confidence(text, max_laws=25)
            doc_level_laws_cleaned = validator.clean_laws(laws_raw)
        
        # ============================================
        # CHUNK 0: SUMMARY
        # ============================================
        summary_parts = []
        
        if title:
            summary_parts.append(f"# {title}")
        
        if abstract:
            summary_parts.append(f"\n## Abstract\n{abstract}")
        
        if sintesi:
            summary_parts.append(f"\n## Sintesi\n{sintesi}")
        
        if keywords:
            summary_parts.append(f"\n## Keywords\n{', '.join(keywords)}")
        
        if category:
            summary_parts.append(f"\n## Categoria\n{category}")
        
        summary_text = "\n".join(summary_parts)
        
        if summary_text and len(summary_text) > 50:
            # Extract references
            refs = extract_references(summary_text)
            
            # ✅ Cited laws (stringhe per backward compat)
            cited_laws_strings = [
                law.get('stringa_originale', '') for law in doc_level_laws_cleaned
                if law.get('stringa_originale')
            ]
            
            # ✅ Full citations (nuove, più complete)
            full_citations = [
                law.get('full_citation', '') for law in doc_level_laws_cleaned
                if law.get('full_citation')
            ]
            
            chunk_metas.append({
                "id": hashlib.md5(f"{doc_id}_summary".encode()).hexdigest(),
                "text": summary_text,
                "payload": {
                    "doc_id": doc_id,
                    "neo4j_id": doc_id,
                    "source": source,
                    "chunk_type": "summary",
                    "title": title,
                    "abstract": abstract,
                    "data_pubblicazione": web_meta.get("data_pubblicazione", ""),
                    "pdf_url": web_meta.get("pdf_url", ""),
                    "chunk_index": 0,
                    "total_chunks": 1,
                    "text": summary_text,
                    "sintesi": sintesi,
                    "categoria_principale": category,
                    
                    # ✅ LAWS: Structured + Strings + Full Citations
                    "cited_laws": cited_laws_strings,
                    "full_citations": full_citations,  # NUOVO
                    "laws_structured": doc_level_laws_cleaned,
                    
                    # ✅ Neo4j-ready IDs
                    "neo4j_law_ids": [law['law_id'] for law in doc_level_laws_cleaned if 'law_id' in law],
                    "neo4j_law_urns": [law['urn'] for law in doc_level_laws_cleaned if 'urn' in law],
                    "neo4j_law_prefixes": [law['qdrant_prefix'] for law in doc_level_laws_cleaned if 'qdrant_prefix' in law],
                    "neo4j_full_citations": full_citations,  # NUOVO
                    
                    # ✅ Metadata
                    "has_law_citations": len(doc_level_laws_cleaned) > 0,
                    "num_law_citations": len(doc_level_laws_cleaned),
                    "num_unique_laws": len(set(law['law_id'] for law in doc_level_laws_cleaned if 'law_id' in law)),
                    
                    # Keep old for compatibility
                    "articoli_legge": cited_laws_strings[:10],
                    
                    "categorie_professionali": sem_meta.get("categorie_professionali", []),
                    "parole_chiave": keywords,
                    "temi_principali": sem_meta.get("temi_principali", []),
                    "confidence_score": float(sem_meta.get("confidence", {}).get("overall_score", 0.5)),
                    "num_pages": content.get("num_pages", 0),
                    "num_tables": content.get("num_tables", 0),
                    "num_figures": content.get("num_figures", 0),
                    "has_images": len(content.get("figures", [])) > 0,
                    "has_tables": len(content.get("tables", [])) > 0,
                    
                    # REFERENCES
                    "referenced_figures": refs['figures'],
                    "referenced_tables": refs['tables'],
                    "has_references": bool(refs['figures'] or refs['tables']),
                }
            })
        
        # ============================================
        # CHUNKS 1+: CONTENT
        # ============================================
        content_chunks = chunk_text_smart(text, CHUNK_SIZE, OVERLAP)
        
        if len(content_chunks) > MAX_CHUNKS_PER_DOC - 1:
            content_chunks = sample_chunks_intelligent(content_chunks, MAX_CHUNKS_PER_DOC - 1)
        
        for idx, chunk_text in enumerate(content_chunks, start=1):
            if not chunk_text or len(chunk_text) < 100:
                continue
            
            # Extract references
            refs = extract_references(chunk_text)
            
            chunk_metas.append({
                "id": hashlib.md5(f"{doc_id}_{idx}".encode()).hexdigest(),
                "text": chunk_text,
                "payload": {
                    "doc_id": doc_id,
                    "neo4j_id": doc_id,
                    "source": source,
                    "chunk_type": "content",
                    "title": title,
                    "abstract": abstract,
                    "data_pubblicazione": web_meta.get("data_pubblicazione", ""),
                    "pdf_url": web_meta.get("pdf_url", ""),
                    "chunk_index": idx,
                    "total_chunks": len(content_chunks) + 1,
                    "text": chunk_text,
                    "sintesi": sintesi,
                    "categoria_principale": category,
                    
                    # ✅ LAWS (same structure as summary)
                    "cited_laws": cited_laws_strings,
                    "full_citations": full_citations,
                    "laws_structured": doc_level_laws_cleaned,
                    
                    # ✅ Neo4j-ready
                    "neo4j_law_ids": [law['law_id'] for law in doc_level_laws_cleaned if 'law_id' in law],
                    "neo4j_law_urns": [law['urn'] for law in doc_level_laws_cleaned if 'urn' in law],
                    "neo4j_law_prefixes": [law['qdrant_prefix'] for law in doc_level_laws_cleaned if 'qdrant_prefix' in law],
                    "neo4j_full_citations": full_citations,
                    
                    "has_law_citations": len(doc_level_laws_cleaned) > 0,
                    "num_law_citations": len(doc_level_laws_cleaned),
                    "num_unique_laws": len(set(law['law_id'] for law in doc_level_laws_cleaned if 'law_id' in law)),
                    
                    "articoli_legge": cited_laws_strings[:10],
                    
                    "categorie_professionali": sem_meta.get("categorie_professionali", []),
                    "parole_chiave": keywords,
                    "temi_principali": sem_meta.get("temi_principali", []),
                    "confidence_score": float(sem_meta.get("confidence", {}).get("overall_score", 0.5)),
                    "num_pages": content.get("num_pages", 0),
                    "num_tables": content.get("num_tables", 0),
                    "num_figures": content.get("num_figures", 0),
                    "has_images": len(content.get("figures", [])) > 0,
                    "has_tables": len(content.get("tables", [])) > 0,
                    
                    # REFERENCES
                    "referenced_figures": refs['figures'],
                    "referenced_tables": refs['tables'],
                    "has_references": bool(refs['figures'] or refs['tables']),
                }
            })
        
        # Update summary total_chunks
        if chunk_metas:
            chunk_metas[0]["payload"]["total_chunks"] = len(chunk_metas)
        
        return chunk_metas
    
    except Exception as e:
        print(f"\n⚠️  Error processing {json_path.name}: {str(e)[:100]}")
        return []


def index_text_chunks(client: QdrantClient, model: SentenceTransformer, delete_first: bool = False):
    """Index text chunks with VALIDATED citations - FIXED VERSION"""
    collection_name = config.COLLECTIONS["technical_text"]
    
    print("\n" + "="*80)
    print(f"🚀 INDEXING TEXT (FIXED LAWS): {collection_name}")
    print(f"⚙️  CHUNK_SIZE={CHUNK_SIZE}, OVERLAP={OVERLAP}, MAX_CHUNKS={MAX_CHUNKS_PER_DOC}")
    if delete_first:
        print("🗑️  DELETE FIRST MODE: Will clear collection before indexing")
    print("="*80 + "\n")
    
    # ✅ Opzione per cancellare prima
    if delete_first:
        try:
            info = client.get_collection(collection_name)
            old_count = info.points_count
            print(f"⚠️  Deleting {old_count:,} existing points from {collection_name}...")
            
            # Ricrea la collection
            from qdrant_client.models import VectorParams, Distance
            client.delete_collection(collection_name)
            client.create_collection(
                collection_name=collection_name,
                vectors_config=VectorParams(
                    size=config.EMBEDDING_DIM,
                    distance=Distance.COSINE
                )
            )
            print(f"✅ Collection recreated\n")
        except Exception as e:
            print(f"⚠️  Could not delete collection: {e}\n")
    
    print("📦 STEP 1/3: Extracting chunks with fixed citations...\n")
    
    # ✅ CREATE VALIDATOR + EXTRACTOR (once!)
    validator = LawValidator()
    extractor = EnhancedLawExtractor()
    
    all_chunk_metas = []
    stats = {
        "truncated_docs": 0, 
        "total_docs": 0,
        "validation_errors": 0,
        "validation_warnings": 0,
        "total_citations": 0,
        "total_unique_laws": 0,
        "docs_with_citations": 0,
    }
    
    print("📘 Processing INAIL...")
    inail_files = sorted(list(config.INAIL_ENRICHED.glob("*.json")))
    
    for json_file in tqdm(inail_files, desc="INAIL"):
        chunks = extract_hierarchical_chunks(json_file, source="INAIL", extractor=extractor, validator=validator)
        all_chunk_metas.extend(chunks)
        stats["total_docs"] += 1
        if len(chunks) >= MAX_CHUNKS_PER_DOC:
            stats["truncated_docs"] += 1
        # Track citation stats
        if chunks and chunks[0]["payload"]["num_law_citations"] > 0:
            stats["docs_with_citations"] += 1
            stats["total_citations"] += chunks[0]["payload"]["num_law_citations"]
            stats["total_unique_laws"] += chunks[0]["payload"]["num_unique_laws"]
    
    inail_count = len(all_chunk_metas)
    print(f"✅ INAIL: {inail_count:,} chunks\n")
    
    print("📗 Processing OSHA...")
    osha_files = sorted(list(config.OSHA_ENRICHED.glob("*.json")))
    
    for json_file in tqdm(osha_files, desc="OSHA"):
        chunks = extract_hierarchical_chunks(json_file, source="OSHA", extractor=extractor, validator=validator)
        all_chunk_metas.extend(chunks)
        stats["total_docs"] += 1
        if len(chunks) >= MAX_CHUNKS_PER_DOC:
            stats["truncated_docs"] += 1
        if chunks and chunks[0]["payload"]["num_law_citations"] > 0:
            stats["docs_with_citations"] += 1
            stats["total_citations"] += chunks[0]["payload"]["num_law_citations"]
            stats["total_unique_laws"] += chunks[0]["payload"]["num_unique_laws"]
    
    osha_count = len(all_chunk_metas) - inail_count
    print(f"✅ OSHA: {osha_count:,} chunks\n")
    
    # Validation stats
    stats["validation_errors"] = len(validator.errors)
    stats["validation_warnings"] = len(validator.warnings)
    
    # Statistics
    avg_chunks = len(all_chunk_metas) / stats["total_docs"] if stats["total_docs"] > 0 else 0
    summary_count = sum(1 for cm in all_chunk_metas if cm["payload"]["chunk_type"] == "summary")
    with_refs = sum(1 for cm in all_chunk_metas if cm["payload"]["has_references"])
    with_citations = sum(1 for cm in all_chunk_metas if cm["payload"]["num_law_citations"] > 0)
    
    print(f"📊 Statistics:")
    print(f"   Total docs: {stats['total_docs']:,}")
    print(f"   Avg chunks/doc: {avg_chunks:.1f}")
    print(f"   Truncated docs (≥{MAX_CHUNKS_PER_DOC} chunks): {stats['truncated_docs']}")
    print(f"   Summary chunks: {summary_count:,}")
    print(f"   Content chunks: {len(all_chunk_metas) - summary_count:,}")
    print(f"   Chunks with fig/tab references: {with_refs:,}")
    print()
    print(f"📜 LAW CITATION Statistics:")
    print(f"   Docs with law citations: {stats['docs_with_citations']:,} ({stats['docs_with_citations']/stats['total_docs']*100:.1f}%)")
    print(f"   Total citations: {stats['total_citations']:,}")
    print(f"   Total unique laws referenced: {stats['total_unique_laws']:,}")
    print(f"   Avg citations/doc (with citations): {stats['total_citations']/stats['docs_with_citations']:.1f}" if stats['docs_with_citations'] > 0 else "")
    print()
    print(f"   Validation errors: {stats['validation_errors']}")
    print(f"   Validation warnings: {stats['validation_warnings']}\n")
    
    if validator.errors:
        print("⚠️  First 5 validation errors:")
        for err in validator.errors[:5]:
            print(f"   - {err}")
        print()
    
    if not all_chunk_metas:
        print("\n❌ NO CHUNKS!\n")
        return 0
    
    print(f"🚀 STEP 2/3: Encoding {len(all_chunk_metas):,} chunks (GPU)...\n")
    
    texts = [cm["text"] for cm in all_chunk_metas]
    embeddings = model.encode(
        texts,
        batch_size=128,
        show_progress_bar=True,
        normalize_embeddings=True
    )
    
    print(f"✅ Encoded {len(embeddings):,} vectors\n")
    
    print(f"⬆️  STEP 3/3: Uploading to Qdrant (UPSERT = update if exists)...\n")
    
    for i in tqdm(range(0, len(all_chunk_metas), config.BATCH_SIZE), desc="Upload"):
        batch = all_chunk_metas[i:i+config.BATCH_SIZE]
        batch_emb = embeddings[i:i+config.BATCH_SIZE]
        
        points = [
            PointStruct(
                id=meta["id"],
                vector=emb.tolist(),
                payload=meta["payload"]
            )
            for meta, emb in zip(batch, batch_emb)
        ]
        
        # ✅ UPSERT: aggiorna se esiste, inserisce se nuovo
        client.upsert(collection_name=collection_name, points=points)
    
    print(f"\n✅ Uploaded {len(all_chunk_metas):,} chunks!\n")
    
    return len(all_chunk_metas)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Index text chunks with fixed law citations')
    parser.add_argument('--delete-first', action='store_true', 
                        help='Delete all existing points before indexing (clean reindex)')
    args = parser.parse_args()
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    client = QdrantClient(host=config.QDRANT_HOST, port=config.QDRANT_PORT)
    
    print(f"\n🔧 Loading model: {config.MODEL_NAME}")
    print(f"🖥️  Device: {device.upper()}")
    
    model = SentenceTransformer(config.MODEL_NAME, device=device)
    print(f"✅ Model loaded (dim={config.EMBEDDING_DIM})\n")
    
    count = index_text_chunks(client, model, delete_first=args.delete_first)
    
    print(f"\n{'='*80}")
    print(f"✅ TEXT INDEXING COMPLETE! ({count:,} vectors)")
    print(f"{'='*80}\n")
