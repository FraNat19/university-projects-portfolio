#!/usr/bin/env python3



"""

📸 IMAGES INDEXING - FILTERED VERSION



IMPROVEMENTS:

1. Filter out logos, covers, and low-quality images

2. Detect meaningful technical content

3. Better confidence scoring

4. Handle truncated descriptions

"""



import json

import hashlib

import re

from pathlib import Path

from typing import List, Dict, Tuple, Optional

from tqdm import tqdm

import sys

import os

import torch



sys.path.insert(0, os.path.dirname(__file__))

import importlib.util

spec = importlib.util.spec_from_file_location("config", "01_config.py")

config = importlib.util.module_from_spec(spec)

spec.loader.exec_module(config)



from qdrant_client import QdrantClient

from qdrant_client.models import PointStruct

from sentence_transformers import SentenceTransformer

from doc_id_utils import clean_doc_id

# Paths

INAIL_IMAGES = Path("/mnt/beegfs/proj/dss.dmaia/INAILProj/images")

OSHA_IMAGES = Path("/mnt/beegfs/proj/dss.dmaia/OSHAProj/images")

INAIL_JSON = Path("/mnt/beegfs/proj/dss.dmaia/INAILProj/enriched_json")

OSHA_JSON = Path("/mnt/beegfs/proj/dss.dmaia/OSHAProj/enriched_json")





# 🚨 FILTER KEYWORDS - immagini da ESCLUDERE

NOISE_KEYWORDS = [

    # Loghi e branding

    'logo', 'inail', 'osha', 'brand', 'trademark',

    # Copertine e frontespizi

    'cover page', 'title page', 'front cover', 'copertina',

    # Elementi generici

    'document title', 'header', 'footer', 'page number',

    # Contenuto vuoto/generico

    'blank', 'empty', 'white background', 'plain',

    # Testo generico (non tecnico)

    'simple text', 'paragraph', 'just text',

]



# ✅ QUALITY KEYWORDS - immagini da INCLUDERE

QUALITY_KEYWORDS = [

    # Elementi tecnici

    'diagram', 'chart', 'graph', 'schematic', 'blueprint',

    'technical drawing', 'floor plan', 'layout',

    # Sicurezza

    'hazard', 'warning', 'safety', 'equipment', 'ppe',

    'emergency', 'fire', 'chemical', 'radiation',

    # Visualizzazioni dati

    'table', 'matrix', 'flowchart', 'process',

    'statistics', 'data', 'measurement',

    # Foto tecniche

    'machinery', 'tools', 'equipment', 'workplace',

    'construction', 'industrial',

]





def clean_caption(caption: str) -> str:

    """Clean malformed captions"""

    if not caption:

        return ""

    

    caption = re.sub(r'<!--.*?-->', '', caption, flags=re.DOTALL)

    caption = caption.strip()

    

    if len(caption) < 20:

        return ""

    

    if re.match(r'^[-•*\d+\.]\s', caption):

        if len(caption.split()) < 5:

            return ""

    

    if caption.isupper() and len(caption) < 50:

        return ""

    

    return caption





def assess_caption_quality(caption: str) -> str:

    """Assess caption quality"""

    cleaned = clean_caption(caption)

    

    if not cleaned:

        return "none"

    

    if len(cleaned) > 100:

        return "high"

    elif len(cleaned) > 40:

        return "medium"

    else:

        return "low"





def is_technical_content(vlm_description: str, caption: str) -> bool:

    """

    🎯 FILTER: Determine if image contains useful technical content

    

    Returns True if image should be indexed, False if it's noise

    """

    text = (vlm_description + " " + caption).lower()

    

    # Check for noise keywords (logos, covers, etc.)

    noise_score = sum(1 for kw in NOISE_KEYWORDS if kw in text)

    

    # Check for quality keywords (technical content)

    quality_score = sum(1 for kw in QUALITY_KEYWORDS if kw in text)

    

    # 🚨 REJECT if:

    # - High noise and low quality

    if noise_score >= 2 and quality_score == 0:

        return False

    

    # - Only mentions "logo" or "title"

    if any(word in text for word in ['only logo', 'just a logo', 'title page only']):

        return False

    

    # - Very generic description (< 100 chars and no quality keywords)

    if len(vlm_description) < 100 and quality_score == 0:

        return False

    

    # ✅ ACCEPT if:

    # - Has quality keywords

    if quality_score >= 2:

        return True

    

    # - Long description (likely detailed technical content)

    if len(vlm_description) > 200:

        return True

    

    # - Has good caption

    if len(caption) > 50 and assess_caption_quality(caption) in ['medium', 'high']:

        return True

    

    # Default: reject if unsure

    return False





def calculate_image_confidence(

    vlm_description: str,

    caption: str,

    caption_quality: str,

    page_confidence: str

) -> Tuple[float, str]:

    """

    Calculate overall confidence score for image usefulness

    Returns: (score 0-1, label)

    """

    score = 0.0

    

    # VLM description quality (40%)

    if len(vlm_description) > 500:

        score += 0.4

    elif len(vlm_description) > 200:

        score += 0.25

    elif len(vlm_description) > 100:

        score += 0.15

    

    # Caption quality (30%)

    caption_scores = {'high': 0.3, 'medium': 0.2, 'low': 0.1, 'none': 0}

    score += caption_scores.get(caption_quality, 0)

    

    # Technical keywords presence (20%)

    text = (vlm_description + " " + caption).lower()

    quality_keywords_found = sum(1 for kw in QUALITY_KEYWORDS if kw in text)

    score += min(quality_keywords_found * 0.05, 0.2)

    

    # Page estimation confidence (10%)

    page_scores = {'medium_high': 0.1, 'medium': 0.07, 'low': 0.04, 'very_low': 0.02, 'none': 0}

    score += page_scores.get(page_confidence, 0)

    

    # Determine label

    if score >= 0.7:

        label = "high"

    elif score >= 0.5:

        label = "medium"

    elif score >= 0.3:

        label = "low"

    else:

        label = "very_low"

    

    return (round(score, 2), label)





def estimate_page_from_sequence(

    figure_num: int, 

    total_figures: int, 

    num_pages: int

) -> Tuple[int, str]:

    """Estimate page from figure sequence"""

    if num_pages == 0 or total_figures == 0:

        return (0, "none")

    

    ratio = figure_num / total_figures

    estimated_page = int(ratio * num_pages)

    estimated_page = max(1, min(estimated_page, num_pages))

    

    if total_figures < 5:

        confidence = "very_low"

    elif total_figures < 20:

        confidence = "low"

    elif total_figures < 50:

        confidence = "medium"

    else:

        confidence = "medium_high"

    

    return (estimated_page, confidence)





def parse_vlm_txt(txt_path: Path) -> Dict:

    """Parse VLM description file - handle truncated descriptions"""

    try:

        with open(txt_path, 'r', encoding='utf-8') as f:

            content = f.read()

        

        metadata = {}

        lines = content.split('\n')

        

        for line in lines:

            if 'Model:' in line:

                metadata['vlm_model'] = line.split('Model:')[1].strip()

            elif 'Resolution:' in line:

                metadata['vlm_resolution'] = line.split('Resolution:')[1].strip()

        

        if 'VLM ANALYSIS:' in content:

            description = content.split('VLM ANALYSIS:')[-1].strip()

        else:

            description = content.strip()

        

        # Remove separator lines and clean

        description = re.sub(r'-{20,}', '', description)

        description = description.strip('"\'').strip()

        

        # Check if truncated

        is_truncated = description.endswith('..') or len(description) > 1500

        

        return {

            'description': description,

            'model': metadata.get('vlm_model', 'SmolVLM-256M-Instruct'),

            'resolution': metadata.get('vlm_resolution', '512px'),

            'is_truncated': is_truncated

        }

    

    except Exception as e:

        return {

            'description': '',

            'model': 'SmolVLM-256M-Instruct',

            'resolution': '512px',

            'is_truncated': False

        }





def extract_figure_metadata(json_path: Path, figure_num: int) -> Dict:

    """Extract figure caption + document stats"""

    try:

        with open(json_path, 'r', encoding='utf-8') as f:

            data = json.load(f)

        

        doc_content = data.get('document_content', {})

        figures = doc_content.get('figures', [])

        

        num_pages = doc_content.get('num_pages', 0)

        total_figures = len(figures)

        

        for fig in figures:

            if fig.get('figure_id', '').endswith(f"_{figure_num}"):

                return {

                    'caption': fig.get('caption', ''),

                    'position': fig.get('position', 0),

                    'num_pages': num_pages,

                    'total_figures': total_figures

                }

        

        return {

            'caption': '',

            'position': 0,

            'num_pages': num_pages,

            'total_figures': total_figures

        }

    

    except:

        return {

            'caption': '',

            'position': 0,

            'num_pages': 0,

            'total_figures': 0

        }





def process_image_directory(

    doc_dir: Path, 

    source: str, 

    json_dir: Path

) -> Tuple[List[Dict], Dict]:

    """

    Process all images in a document directory

    Returns: (chunks, stats)

    """

    chunks = []

    stats = {

        'total': 0,

        'filtered_out': 0,

        'indexed': 0,

        'reasons': {

            'noise_content': 0,

            'low_quality': 0,

            'no_description': 0,

        }

    }

    

    doc_id = clean_doc_id(doc_dir.name, max_len=80)

    

    enriched_json = None

    for json_file in json_dir.glob(f"{doc_id}*.json"):

        enriched_json = json_file

        break

    

    png_files = sorted(doc_dir.glob("*.png"))

    

    for png_file in png_files:

        stats['total'] += 1

        txt_file = png_file.with_suffix('.txt')

        

        if not txt_file.exists():

            stats['filtered_out'] += 1

            stats['reasons']['no_description'] += 1

            continue

        

        vlm_data = parse_vlm_txt(txt_file)

        

        if not vlm_data['description']:

            stats['filtered_out'] += 1

            stats['reasons']['no_description'] += 1

            continue

        

        match = re.search(r'figure_(\d+)', png_file.stem)

        figure_num = int(match.group(1)) if match else 0

        

        fig_meta = {'caption': '', 'position': 0, 'num_pages': 0, 'total_figures': 0}

        if enriched_json:

            fig_meta = extract_figure_metadata(enriched_json, figure_num)

        

        caption_cleaned = clean_caption(fig_meta['caption'])

        caption_quality = assess_caption_quality(fig_meta['caption'])

        

        # 🎯 APPLY FILTER

        if not is_technical_content(vlm_data['description'], caption_cleaned):

            stats['filtered_out'] += 1

            stats['reasons']['noise_content'] += 1

            continue

        

        estimated_page, page_confidence = estimate_page_from_sequence(

            figure_num,

            fig_meta['total_figures'],

            fig_meta['num_pages']

        )

        

        # Calculate confidence score

        confidence_score, confidence_label = calculate_image_confidence(

            vlm_data['description'],

            caption_cleaned,

            caption_quality,

            page_confidence

        )

        

        # Filter by minimum confidence

        if confidence_score < 0.25:  # Soglia minima

            stats['filtered_out'] += 1

            stats['reasons']['low_quality'] += 1

            continue

        

        caption_text = f"Caption: {caption_cleaned}" if caption_cleaned else "Caption: Non disponibile"

        

        # Use more of the VLM description (up to 1500 chars for embedding)

        vlm_desc_for_embedding = vlm_data['description'][:1500]

        

        embedding_text = f"""

Figura {figure_num} dal documento {doc_id}.

Fonte: {source}



Descrizione visiva (VLM):

{vlm_desc_for_embedding}



{caption_text}



Posizione stimata: pagina ~{estimated_page} di {fig_meta['num_pages']} (figura {figure_num} di {fig_meta['total_figures']})

"""

        

        image_id = f"{doc_id}_figure_{figure_num}"

        

        payload = {

            # IDs

            'image_id': image_id,

            'parent_doc_id': doc_id,

            'neo4j_doc_id': doc_id,

            'figure_num': figure_num,

            

            # Files

            'filename': png_file.name,

            'image_path': str(png_file),

            'vlm_txt_path': str(txt_file),

            

            # VLM data

            'vlm_description': vlm_data['description'][:3000],  # Aumentato

            'vlm_model': vlm_data['model'],

            'vlm_resolution': vlm_data['resolution'],

            'vlm_truncated': vlm_data['is_truncated'],

            

            # Caption

            'caption': caption_cleaned[:500],

            'caption_quality': caption_quality,

            

            # Position

            'position': fig_meta['position'],

            'estimated_page': estimated_page,

            'page_estimation_confidence': page_confidence,

            'page_estimation_method': 'figure_sequence_proportion',

            'total_pages': fig_meta['num_pages'],

            'total_figures': fig_meta['total_figures'],

            

            # Quality scores

            'confidence_score': confidence_score,

            'confidence_label': confidence_label,

            

            # Source

            'source': source,

            

            # Text for retrieval

            'text': embedding_text

        }

        

        chunks.append({

            'id': hashlib.md5(image_id.encode()).hexdigest(),

            'text': embedding_text,

            'payload': payload

        })

        

        stats['indexed'] += 1

    

    return chunks, stats





def index_images(client: QdrantClient, model: SentenceTransformer):

    """Index all images with VLM descriptions - FILTERED VERSION"""

    collection_name = config.COLLECTIONS["technical_images"]

    

    print("\n" + "="*80)

    print(f"📸 INDEXING IMAGES (FILTERED): {collection_name}")

    print("="*80 + "\n")

    

    print("📦 STEP 1/3: Extracting and filtering images...\n")

    

    all_chunks = []

    total_stats = {

        'total': 0,

        'filtered_out': 0,

        'indexed': 0,

        'reasons': {

            'noise_content': 0,

            'low_quality': 0,

            'no_description': 0,

        }

    }

    

    # INAIL

    print("📘 Processing INAIL images...")

    

    for doc_dir in tqdm(sorted([d for d in INAIL_IMAGES.iterdir() if d.is_dir()]), desc="INAIL"):

        chunks, stats = process_image_directory(doc_dir, "INAIL", INAIL_JSON)

        all_chunks.extend(chunks)

        

        for key in total_stats:

            if key == 'reasons':

                for reason_key in total_stats['reasons']:

                    total_stats['reasons'][reason_key] += stats['reasons'][reason_key]

            else:

                total_stats[key] += stats[key]

    

    inail_count = len(all_chunks)

    print(f"✅ INAIL: {inail_count:,} images indexed\n")

    

    # OSHA

    print("📗 Processing OSHA images...")

    

    for doc_dir in tqdm(sorted([d for d in OSHA_IMAGES.iterdir() if d.is_dir()]), desc="OSHA"):

        chunks, stats = process_image_directory(doc_dir, "OSHA", OSHA_JSON)

        all_chunks.extend(chunks)

        

        for key in total_stats:

            if key == 'reasons':

                for reason_key in total_stats['reasons']:

                    total_stats['reasons'][reason_key] += stats['reasons'][reason_key]

            else:

                total_stats[key] += stats[key]

    

    osha_count = len(all_chunks) - inail_count

    print(f"✅ OSHA: {osha_count:,} images indexed\n")

    

    # Print filtering stats

    print("\n" + "="*80)

    print("📊 FILTERING STATISTICS")

    print("="*80)

    print(f"Total images found: {total_stats['total']:,}")

    print(f"Images indexed: {total_stats['indexed']:,} ({100*total_stats['indexed']/max(total_stats['total'],1):.1f}%)")

    print(f"Images filtered out: {total_stats['filtered_out']:,} ({100*total_stats['filtered_out']/max(total_stats['total'],1):.1f})")

    print(f"\nFiltering reasons:")

    print(f"  - Noise content (logos, covers): {total_stats['reasons']['noise_content']:,}")

    print(f"  - Low quality/confidence: {total_stats['reasons']['low_quality']:,}")

    print(f"  - No VLM description: {total_stats['reasons']['no_description']:,}")

    print("="*80 + "\n")

    

    if not all_chunks:

        print("❌ NO IMAGES TO INDEX!\n")

        return 0

    

    # Encoding

    print(f"🚀 STEP 2/3: Encoding {len(all_chunks):,} images (GPU)...\n")

    

    texts = [c['text'] for c in all_chunks]

    embeddings = model.encode(

        texts,

        batch_size=128,

        show_progress_bar=True,

        normalize_embeddings=True

    )

    

    print(f"✅ Encoded {len(embeddings):,} vectors\n")

    

    # Upload

    print(f"⬆️ STEP 3/3: Uploading to Qdrant...\n")

    

    for i in tqdm(range(0, len(all_chunks), 100), desc="Upload"):

        batch = all_chunks[i:i+100]

        batch_emb = embeddings[i:i+100]

        

        points = [

            PointStruct(

                id=chunk['id'],

                vector=emb.tolist(),

                payload=chunk['payload']

            )

            for chunk, emb in zip(batch, batch_emb)

        ]

        

        client.upsert(collection_name=collection_name, points=points)

    

    print(f"\n✅ Uploaded {len(all_chunks):,} images!\n")

    

    return len(all_chunks)





if __name__ == "__main__":

    device = "cuda" if torch.cuda.is_available() else "cpu"

    

    client = QdrantClient(host=config.QDRANT_HOST, port=config.QDRANT_PORT)

    

    print(f"\n🔧 Loading model: {config.MODEL_NAME}")

    print(f"🖥️ Device: {device.upper()}")

    

    model = SentenceTransformer(config.MODEL_NAME, device=device)

    print(f"✅ Model loaded\n")

    

    count = index_images(client, model)

    

    print(f"\n{'='*80}")

    print(f"✅ IMAGES INDEXING COMPLETE! ({count:,} vectors)")

    print(f"{'='*80}\n")
