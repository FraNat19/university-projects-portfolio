#!/usr/bin/env python3
"""
✅ LAYER 3 v2: SEMANTIC BRIDGE - ALIGNED WITH ENRICHMENT

MIGLIORAMENTI rispetto a v1:
1. ✅ Usa ESATTAMENTE le stesse categorie di enrichment_inail
2. ✅ Mapping professioni → settori SPECIFICO (non generico)
3. ✅ Preserva la granularità delle professioni
4. ✅ Legge da enriched.laws_structured (priorità) + fallback sem_meta
5. ✅ Opzione --clean per cancellare SOLO relazioni layer 3
"""

import json
import sys
import os
from pathlib import Path
from neo4j import GraphDatabase
from collections import defaultdict
from tqdm import tqdm
import argparse

sys.path.insert(0, os.path.dirname(__file__))

# Try to import doc_id_utils, fallback if not available
try:
    from doc_id_utils import clean_doc_id, extract_doc_id_from_metadata, validate_doc_id
    HAS_DOC_ID_UTILS = True
except ImportError:
    HAS_DOC_ID_UTILS = False
    print("⚠️  doc_id_utils not found, using fallback")

# ======================== CONFIG ========================

NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "thesis2024")

INAIL_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/INAIL_enriched")
OSHA_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/OSHA_enriched")

# Fallback paths
if not INAIL_ENRICHED.exists():
    INAIL_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/INAILProj/enriched_json")
if not OSHA_ENRICHED.exists():
    OSHA_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/OSHAProj/enriched_json")

BATCH_SIZE = 100

driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))


# ======================== DOMAIN KNOWLEDGE ========================
# ESATTAMENTE come in enrichment_inail.py

DOCUMENT_CATEGORIES = [
    "Rischio chimico e cancerogeno",
    "Rischio biologico", 
    "Rischio fisico (rumore, vibrazioni, radiazioni)",
    "Prevenzione incendi ed esplosioni",
    "Sicurezza macchine e attrezzature",
    "Edilizia e cantieri",
    "Agricoltura e silvicoltura",
    "Trasporti e logistica",
    "Settore sanitario",
    "DPI e dispositivi di protezione",
    "Formazione e addestramento",
    "Sorveglianza sanitaria",
    "Normativa e adempimenti",
    "Valutazione e gestione del rischio",
    "Infortuni e malattie professionali",
    "Altro"
]

KNOWN_PROFESSIONS = [
    "Operai edili", "Carpentieri", "Muratori", "Ponteggiatori", 
    "Elettricisti", "Idraulici", "Saldatori",
    "Meccanici industriali", "Operatori macchine utensili", "Manutentori",
    "Agricoltori", "Operatori forestali", "Potatori",
    "Infermieri", "Medici", "OSS", "Tecnici di laboratorio", 
    "Radiologi", "Operatori sanitari",
    "Autisti professionali", "Magazzinieri", "Mulettisti",
    "Tecnici della prevenzione", "RSPP", "ASPP", "HSE Manager",
    "Coordinatori per la sicurezza", "RLS", "Preposti alla sicurezza",
    "Chimici", "Biologi", "Tecnici di laboratorio chimico",
    "Lavoratori esposti ad agenti chimici pericolosi",
    "Lavoratori esposti ad agenti cancerogeni o mutageni",
    "Lavoratori esposti ad agenti biologici",
    "Lavoratori esposti a polveri e fibre minerali",
    "Lavoratori esposti a rumore e vibrazioni",
    "Lavoratori in ambienti confinati",
    "Lavoratori in quota",
    "Addetti antincendio"
]

# ✅ MAPPING SPECIFICO professione → categoria (allineato con enrichment)
PROFESSION_TO_CATEGORY = {
    # Edilizia
    "Operai edili": "Edilizia e cantieri",
    "Carpentieri": "Edilizia e cantieri",
    "Muratori": "Edilizia e cantieri",
    "Ponteggiatori": "Edilizia e cantieri",
    "Lavoratori in quota": "Edilizia e cantieri",
    "Coordinatori per la sicurezza": "Edilizia e cantieri",
    
    # Impianti/Manutenzione
    "Elettricisti": "Sicurezza macchine e attrezzature",
    "Idraulici": "Sicurezza macchine e attrezzature",
    "Saldatori": "Sicurezza macchine e attrezzature",
    "Meccanici industriali": "Sicurezza macchine e attrezzature",
    "Operatori macchine utensili": "Sicurezza macchine e attrezzature",
    "Manutentori": "Sicurezza macchine e attrezzature",
    
    # Agricoltura
    "Agricoltori": "Agricoltura e silvicoltura",
    "Operatori forestali": "Agricoltura e silvicoltura",
    "Potatori": "Agricoltura e silvicoltura",
    
    # Sanitario
    "Infermieri": "Settore sanitario",
    "Medici": "Settore sanitario",
    "OSS": "Settore sanitario",
    "Tecnici di laboratorio": "Settore sanitario",
    "Radiologi": "Settore sanitario",
    "Operatori sanitari": "Settore sanitario",
    
    # Trasporti
    "Autisti professionali": "Trasporti e logistica",
    "Magazzinieri": "Trasporti e logistica",
    "Mulettisti": "Trasporti e logistica",
    
    # Sicurezza/Prevenzione
    "Tecnici della prevenzione": "Valutazione e gestione del rischio",
    "RSPP": "Valutazione e gestione del rischio",
    "ASPP": "Valutazione e gestione del rischio",
    "HSE Manager": "Valutazione e gestione del rischio",
    "RLS": "Valutazione e gestione del rischio",
    "Preposti alla sicurezza": "Valutazione e gestione del rischio",
    
    # Chimico
    "Chimici": "Rischio chimico e cancerogeno",
    "Biologi": "Rischio biologico",
    "Tecnici di laboratorio chimico": "Rischio chimico e cancerogeno",
    "Lavoratori esposti ad agenti chimici pericolosi": "Rischio chimico e cancerogeno",
    "Lavoratori esposti ad agenti cancerogeni o mutageni": "Rischio chimico e cancerogeno",
    "Lavoratori esposti ad agenti biologici": "Rischio biologico",
    "Lavoratori esposti a polveri e fibre minerali": "Rischio chimico e cancerogeno",
    "Lavoratori esposti a rumore e vibrazioni": "Rischio fisico (rumore, vibrazioni, radiazioni)",
    
    # Ambienti speciali
    "Lavoratori in ambienti confinati": "Valutazione e gestione del rischio",
    "Addetti antincendio": "Prevenzione incendi ed esplosioni",
}

# Severity mapping per categoria
CATEGORY_SEVERITY = {
    "Rischio chimico e cancerogeno": "HIGH",
    "Rischio biologico": "HIGH",
    "Rischio fisico (rumore, vibrazioni, radiazioni)": "MEDIUM",
    "Prevenzione incendi ed esplosioni": "HIGH",
    "Sicurezza macchine e attrezzature": "MEDIUM",
    "Edilizia e cantieri": "HIGH",
    "Agricoltura e silvicoltura": "MEDIUM",
    "Trasporti e logistica": "MEDIUM",
    "Settore sanitario": "MEDIUM",
    "DPI e dispositivi di protezione": "LOW",
    "Formazione e addestramento": "LOW",
    "Sorveglianza sanitaria": "LOW",
    "Normativa e adempimenti": "LOW",
    "Valutazione e gestione del rischio": "MEDIUM",
    "Infortuni e malattie professionali": "HIGH",
    "Altro": "LOW"
}


# ======================== HELPER FUNCTIONS ========================

def fallback_clean_doc_id(raw_id: str) -> str:
    """Fallback se doc_id_utils non disponibile"""
    import re
    # Rimuovi timestamp YYYYMMDD_HHMMSS
    cleaned = re.sub(r'_\d{8}_\d{6}', '', raw_id)
    cleaned = cleaned.strip('_')
    if len(cleaned) > 80:
        cleaned = cleaned[:80].rsplit('_', 1)[0]
    return cleaned


def fallback_extract_doc_id(web_meta: dict, json_path: Path) -> str:
    """Fallback se doc_id_utils non disponibile"""
    # Priorità: pdf_filename > title > filename
    pdf_url = web_meta.get('pdf_url', '')
    if pdf_url:
        pdf_name = pdf_url.split('/')[-1].replace('.pdf', '')
        return fallback_clean_doc_id(pdf_name)
    
    title = web_meta.get('title', '')
    if title:
        import re
        clean_title = re.sub(r'[^\w\s-]', '', title)
        clean_title = re.sub(r'\s+', '_', clean_title)
        return fallback_clean_doc_id(clean_title[:80])
    
    return fallback_clean_doc_id(json_path.stem)


def get_doc_id(web_meta: dict, json_path: Path) -> str:
    """Get doc_id using doc_id_utils or fallback"""
    if HAS_DOC_ID_UTILS:
        return extract_doc_id_from_metadata(web_meta, json_path)
    return fallback_extract_doc_id(web_meta, json_path)


def is_valid_doc_id(doc_id: str) -> bool:
    """Validate doc_id"""
    if HAS_DOC_ID_UTILS:
        return validate_doc_id(doc_id)
    return doc_id and len(doc_id) >= 5 and len(doc_id) <= 200


# ======================== SCHEMA ========================

def ensure_schema():
    """Setup schema constraints"""
    
    print("\n🔧 Ensuring schema constraints...")
    
    with driver.session() as session:
        constraints = [
            ("TechnicalDocument_doc_id", 
             "CREATE CONSTRAINT technical_document_doc_id IF NOT EXISTS FOR (d:TechnicalDocument) REQUIRE d.doc_id IS UNIQUE"),
            ("RiskCategory_name", 
             "CREATE CONSTRAINT RiskCategory_name IF NOT EXISTS FOR (r:RiskCategory) REQUIRE r.name IS UNIQUE"),
            ("ProfessionalTag_categoria", 
             "CREATE CONSTRAINT ProfessionalTag_categoria IF NOT EXISTS FOR (p:ProfessionalTag) REQUIRE p.categoria IS UNIQUE"),
            ("Topic_name", 
             "CREATE CONSTRAINT Topic_name IF NOT EXISTS FOR (t:Topic) REQUIRE t.name IS UNIQUE"),
        ]
        
        for name, query in constraints:
            try:
                session.run(query)
                print(f"   ✅ {name}")
            except Exception as e:
                if "already exists" in str(e).lower():
                    print(f"   ⭐ {name} (exists)")
    
    print()


# ======================== CLEANUP ========================

def clean_layer3_only():
    """
    ✅ Cancella SOLO nodi e relazioni create da Layer 3
    NON tocca: CanonicalLaw, Article, TechnicalDocument, CITES, HAS_ARTICLE
    """
    
    print("\n🧹 Cleaning Layer 3 data ONLY...")
    
    with driver.session() as session:
        # 1. Cancella relazioni Layer 3
        relations_to_delete = [
            ("ADDRESSES", "TechnicalDocument → RiskCategory"),
            ("TARGETS", "TechnicalDocument → ProfessionalTag"),
            ("DISCUSSES", "TechnicalDocument → Topic"),
            ("GOVERNED_BY", "RiskCategory → CanonicalLaw"),
            ("REGULATED_BY", "ProfessionalTag → CanonicalLaw"),
        ]
        
        for rel_type, desc in relations_to_delete:
            result = session.run(f"""
                MATCH ()-[r:{rel_type}]->()
                DELETE r
                RETURN count(r) as deleted
            """)
            cnt = result.single()['deleted']
            print(f"   🗑️  {rel_type}: {cnt:,} deleted ({desc})")
        
        # 2. Cancella nodi Layer 3 (orfani)
        nodes_to_delete = [
            ("RiskCategory", "Risk categories"),
            ("ProfessionalTag", "Professional tags"),
            ("Topic", "Topics"),
        ]
        
        for node_type, desc in nodes_to_delete:
            result = session.run(f"""
                MATCH (n:{node_type})
                DELETE n
                RETURN count(n) as deleted
            """)
            cnt = result.single()['deleted']
            print(f"   🗑️  {node_type}: {cnt:,} deleted ({desc})")
    
    print("\n✅ Layer 3 cleanup complete!\n")


# ======================== ENTITY EXTRACTOR ========================

class SemanticEntityExtractorV2:
    """
    ✅ Estrattore migliorato che usa le stesse categorie di enrichment_inail
    """
    
    def __init__(self):
        self.risk_categories = defaultdict(int)
        self.professional_tags = defaultdict(int)
        self.topics = defaultdict(int)
        self.law_to_entities = defaultdict(lambda: {'risks': set(), 'professions': set()})
        
        # Stats
        self.processed_count = 0
        self.error_count = 0
        self.invalid_docid_count = 0
    
    def extract_from_json(self, json_path: Path, source: str) -> dict:
        """
        ✅ Estrae entità semantiche preservando la granularità di enrichment
        """
        
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except Exception as e:
            self.error_count += 1
            return None
        
        web_meta = data.get('web_metadata', {})
        sem_meta = data.get('semantic_metadata', {})
        enriched = data.get('enriched', {})
        
        # Get doc_id
        doc_id = get_doc_id(web_meta, json_path)
        
        if not is_valid_doc_id(doc_id):
            self.invalid_docid_count += 1
            return None
        
        entities = {
            'doc_id': doc_id,
            'source': source,
            'risk_categories': [],
            'professional_tags': [],
            'topics': [],
            'cited_laws': []
        }
        
        # ============================================
        # 1. RISK CATEGORIES (categoria_principale + temi)
        # ============================================
        
        main_cat = sem_meta.get('categoria_principale', '')
        
        # ✅ Valida che sia una categoria nota
        if main_cat and main_cat in DOCUMENT_CATEGORIES:
            severity = CATEGORY_SEVERITY.get(main_cat, 'MEDIUM')
            entities['risk_categories'].append({
                'name': main_cat,
                'type': 'main_category',
                'severity': severity,
                'source': source
            })
            self.risk_categories[main_cat] += 1
        elif main_cat and len(main_cat) > 3:
            # Categoria custom (non standard)
            entities['risk_categories'].append({
                'name': main_cat.strip(),
                'type': 'custom_category',
                'severity': 'MEDIUM',
                'source': source
            })
            self.risk_categories[main_cat.strip()] += 1
        
        # Temi principali come categorie secondarie
        for tema in sem_meta.get('temi_principali', [])[:5]:
            if tema and len(tema) > 3:
                # Inferisci severity dal nome del tema
                severity = self._infer_severity(tema)
                entities['risk_categories'].append({
                    'name': tema.strip(),
                    'type': 'theme',
                    'severity': severity,
                    'source': source
                })
                self.risk_categories[tema.strip()] += 1
        
        # ============================================
        # 2. PROFESSIONAL TAGS (con mapping specifico)
        # ============================================
        
        for prof in sem_meta.get('categorie_professionali', [])[:10]:
            if prof and len(prof) > 3:
                # ✅ Usa mapping specifico invece di inferenza generica
                related_category = PROFESSION_TO_CATEGORY.get(
                    prof.strip(), 
                    main_cat or "Altro"
                )
                
                entities['professional_tags'].append({
                    'categoria': prof.strip(),
                    'related_category': related_category,
                    'source': source
                })
                self.professional_tags[prof.strip()] += 1
        
        # ============================================
        # 3. TOPICS (parole chiave)
        # ============================================
        
        for kw in sem_meta.get('parole_chiave', [])[:15]:
            if kw and len(kw) > 2:
                entities['topics'].append({
                    'name': kw.lower().strip(),
                    'source': source
                })
                self.topics[kw.lower().strip()] += 1
        
        # ============================================
        # 4. CITED LAWS (priorità a enriched.laws_structured)
        # ============================================
        
        # Priorità: enriched.laws_structured > sem_meta.laws_structured
        laws_structured = enriched.get('laws_structured', [])
        if not laws_structured:
            laws_structured = sem_meta.get('laws_structured', [])
        
        for law in laws_structured[:30]:
            law_id = law.get('law_id')
            if not law_id:
                continue
            
            entities['cited_laws'].append(law_id)
            
            # Track entity-law associations
            for risk in entities['risk_categories']:
                self.law_to_entities[law_id]['risks'].add(risk['name'])
            for prof in entities['professional_tags']:
                self.law_to_entities[law_id]['professions'].add(prof['categoria'])
        
        entities['cited_laws'] = list(set(entities['cited_laws']))
        
        self.processed_count += 1
        return entities
    
    def _infer_severity(self, text: str) -> str:
        """Inferisci severity dal testo"""
        text_lower = text.lower()
        
        high_keywords = [
            'morte', 'grave', 'caduta', 'esplosion', 'incendi', 
            'elettric', 'crollo', 'infortuni', 'tossic', 'cancerogen',
            'amianto', 'chimico', 'biologico'
        ]
        medium_keywords = [
            'lesion', 'danno', 'esposizion', 'stress', 'disturb',
            'rumore', 'vibrazion', 'ergonomic'
        ]
        
        if any(kw in text_lower for kw in high_keywords):
            return 'HIGH'
        elif any(kw in text_lower for kw in medium_keywords):
            return 'MEDIUM'
        return 'LOW'
    
    def print_stats(self):
        """Stampa statistiche estrazione"""
        print(f"\n📊 Extraction Stats:")
        print(f"   Processed: {self.processed_count:,}")
        print(f"   Errors: {self.error_count:,}")
        print(f"   Invalid doc_id: {self.invalid_docid_count:,}")


# ======================== NEO4J OPS ========================

def create_semantic_nodes_batch(session, entities_batch, entity_type):
    """Crea nodi semantici in batch"""
    
    if not entities_batch:
        return 0
    
    if entity_type == 'risk':
        query = """
        UNWIND $entities as entity
        MERGE (r:RiskCategory {name: entity.name})
        ON CREATE SET
            r.type = entity.type,
            r.severity_level = entity.severity,
            r.source = entity.source,
            r.frequency = 1,
            r.layer = 3,
            r.created_at = datetime()
        ON MATCH SET
            r.frequency = r.frequency + 1
        RETURN count(r) as cnt
        """
    
    elif entity_type == 'profession':
        query = """
        UNWIND $entities as entity
        MERGE (p:ProfessionalTag {categoria: entity.categoria})
        ON CREATE SET
            p.related_category = entity.related_category,
            p.source = entity.source,
            p.frequency = 1,
            p.layer = 3,
            p.created_at = datetime()
        ON MATCH SET
            p.frequency = p.frequency + 1
        RETURN count(p) as cnt
        """
    
    elif entity_type == 'topic':
        query = """
        UNWIND $entities as entity
        MERGE (t:Topic {name: entity.name})
        ON CREATE SET
            t.source = entity.source,
            t.frequency = 1,
            t.layer = 3,
            t.created_at = datetime()
        ON MATCH SET
            t.frequency = t.frequency + 1
        RETURN count(t) as cnt
        """
    
    result = session.run(query, entities=entities_batch)
    return result.single()['cnt']


def link_document_to_entities(session, doc_id, risks, professions, topics):
    """
    ✅ Link documento a entità semantiche
    """
    
    if not risks and not professions and not topics:
        return 0, False
    
    # Prima verifica che il documento esista
    check_result = session.run("""
        MATCH (doc:TechnicalDocument {doc_id: $doc_id})
        RETURN doc.doc_id as id
    """, doc_id=doc_id)
    
    if not check_result.single():
        return 0, False
    
    total_links = 0
    
    # Link Risks
    if risks:
        result = session.run("""
            MATCH (doc:TechnicalDocument {doc_id: $doc_id})
            UNWIND $risks as risk
            MATCH (r:RiskCategory {name: risk.name})
            MERGE (doc)-[rel:ADDRESSES]->(r)
            ON CREATE SET 
                rel.type = risk.type,
                rel.layer = 3,
                rel.created_at = datetime()
            RETURN count(rel) as cnt
        """, doc_id=doc_id, risks=risks)
        total_links += result.single()['cnt']
    
    # Link Professions
    if professions:
        result = session.run("""
            MATCH (doc:TechnicalDocument {doc_id: $doc_id})
            UNWIND $professions as prof
            MATCH (p:ProfessionalTag {categoria: prof.categoria})
            MERGE (doc)-[rel:TARGETS]->(p)
            ON CREATE SET
                rel.layer = 3,
                rel.created_at = datetime()
            RETURN count(rel) as cnt
        """, doc_id=doc_id, professions=professions)
        total_links += result.single()['cnt']
    
    # Link Topics
    if topics:
        result = session.run("""
            MATCH (doc:TechnicalDocument {doc_id: $doc_id})
            UNWIND $topics as topic
            MATCH (t:Topic {name: topic.name})
            MERGE (doc)-[rel:DISCUSSES]->(t)
            ON CREATE SET
                rel.layer = 3,
                rel.created_at = datetime()
            RETURN count(rel) as cnt
        """, doc_id=doc_id, topics=topics)
        total_links += result.single()['cnt']
    
    return total_links, True


def create_entity_law_bridges(session, law_to_entities):
    """Crea bridge tra entità semantiche e leggi"""
    
    if not law_to_entities:
        return 0, 0
    
    # Prepare mappings
    mappings = []
    for law_id, entities in law_to_entities.items():
        if entities['risks'] or entities['professions']:
            mappings.append({
                'law_id': law_id,
                'risks': list(entities['risks']),
                'professions': list(entities['professions'])
            })
    
    if not mappings:
        return 0, 0
    
    # Risk → Law (GOVERNED_BY)
    try:
        result_risk = session.run("""
            UNWIND $mappings as mapping
            MATCH (law:CanonicalLaw {law_id: mapping.law_id})
            WITH DISTINCT law, mapping
            UNWIND mapping.risks as risk_name
            MATCH (r:RiskCategory {name: risk_name})
            MERGE (r)-[rel:GOVERNED_BY]->(law)
            ON CREATE SET 
                rel.created_at = datetime(), 
                rel.confidence = 0.9,
                rel.layer = 3
            RETURN count(rel) as cnt
        """, mappings=mappings)
        cnt_risk = result_risk.single()['cnt']
    except Exception as e:
        print(f"  ⚠️ Risk bridge error: {str(e)[:100]}")
        cnt_risk = 0
    
    # Profession → Law (REGULATED_BY)
    try:
        result_prof = session.run("""
            UNWIND $mappings as mapping
            MATCH (law:CanonicalLaw {law_id: mapping.law_id})
            WITH DISTINCT law, mapping
            UNWIND mapping.professions as prof_name
            MATCH (p:ProfessionalTag {categoria: prof_name})
            MERGE (p)-[rel:REGULATED_BY]->(law)
            ON CREATE SET 
                rel.created_at = datetime(), 
                rel.confidence = 0.9,
                rel.layer = 3
            RETURN count(rel) as cnt
        """, mappings=mappings)
        cnt_prof = result_prof.single()['cnt']
    except Exception as e:
        print(f"  ⚠️ Prof bridge error: {str(e)[:100]}")
        cnt_prof = 0
    
    return cnt_risk, cnt_prof


# ======================== MAIN PIPELINE ========================

def extract_semantic_layer(source='INAIL'):
    """Main extraction pipeline"""
    
    print(f"\n🚀 LAYER 3 v2: SEMANTIC EXTRACTION ({source})")
    print("="*70)
    
    data_dir = INAIL_ENRICHED if source == 'INAIL' else OSHA_ENRICHED
    
    if not data_dir.exists():
        print(f"❌ Directory not found: {data_dir}")
        return
    
    json_files = sorted(list(data_dir.glob("*.json")))
    
    print(f"\n📂 Found {len(json_files):,} JSON files in {data_dir}")
    
    if len(json_files) == 0:
        print(f"❌ No files found!")
        return
    
    # Extract entities
    print(f"\n📄 Step 1/4: Extracting semantic entities...\n")
    
    extractor = SemanticEntityExtractorV2()
    all_entities = []
    
    for json_file in tqdm(json_files, desc="Parsing"):
        entities = extractor.extract_from_json(json_file, source)
        if entities:
            all_entities.append(entities)
    
    extractor.print_stats()
    
    print(f"\n✅ Extracted entities from {len(all_entities):,} documents")
    print(f"   Unique risk categories: {len(extractor.risk_categories)}")
    print(f"   Unique professions: {len(extractor.professional_tags)}")
    print(f"   Unique topics: {len(extractor.topics)}")
    
    # Show top categories
    print(f"\n📊 Top 10 Risk Categories:")
    for cat, cnt in sorted(extractor.risk_categories.items(), key=lambda x: -x[1])[:10]:
        print(f"   {cat}: {cnt}")
    
    print(f"\n👷 Top 10 Professions:")
    for prof, cnt in sorted(extractor.professional_tags.items(), key=lambda x: -x[1])[:10]:
        print(f"   {prof}: {cnt}")
    
    # Create nodes
    print(f"\n📄 Step 2/4: Creating semantic nodes...\n")
    
    with driver.session() as session:
        # Risks
        risk_entities = []
        for doc_entities in all_entities:
            risk_entities.extend(doc_entities['risk_categories'])
        
        if risk_entities:
            cnt_risk = create_semantic_nodes_batch(session, risk_entities, 'risk')
            print(f"   ✅ RiskCategory: {cnt_risk:,} nodes")
        
        # Professions
        prof_entities = []
        for doc_entities in all_entities:
            prof_entities.extend(doc_entities['professional_tags'])
        
        if prof_entities:
            cnt_prof = create_semantic_nodes_batch(session, prof_entities, 'profession')
            print(f"   ✅ ProfessionalTag: {cnt_prof:,} nodes")
        
        # Topics
        topic_entities = []
        for doc_entities in all_entities:
            topic_entities.extend(doc_entities['topics'])
        
        if topic_entities:
            cnt_topic = create_semantic_nodes_batch(session, topic_entities, 'topic')
            print(f"   ✅ Topic: {cnt_topic:,} nodes")
    
    # Link documents
    print(f"\n📄 Step 3/4: Linking documents to entities...\n")
    
    linked_count = 0
    missing_count = 0
    missing_samples = []
    
    with driver.session() as session:
        for entities in tqdm(all_entities, desc="Linking"):
            cnt, matched = link_document_to_entities(
                session,
                entities['doc_id'],
                entities['risk_categories'],
                entities['professional_tags'],
                entities['topics']
            )
            if matched:
                linked_count += 1
            else:
                missing_count += 1
                if len(missing_samples) < 10:
                    missing_samples.append(entities['doc_id'])
    
    print(f"\n   ✅ Linked {linked_count:,} documents")
    
    if missing_count > 0:
        print(f"   ⚠️  Missing in Neo4j: {missing_count:,} documents")
        print(f"\n   📋 Sample missing doc_ids:")
        for doc_id in missing_samples[:5]:
            print(f"      - {doc_id}")
    
    # Create entity-law bridges
    print(f"\n📄 Step 4/4: Creating entity-law bridges...\n")
    
    with driver.session() as session:
        cnt_risk_law, cnt_prof_law = create_entity_law_bridges(
            session, 
            extractor.law_to_entities
        )
    
    print(f"   ✅ Risk → Law bridges: {cnt_risk_law:,}")
    print(f"   ✅ Profession → Law bridges: {cnt_prof_law:,}")
    
    # Summary
    print(f"\n{'='*70}")
    print(f"✅ LAYER 3 v2 COMPLETE ({source})")
    print(f"{'='*70}\n")


def verify_layer3():
    """Verifica stato Layer 3"""
    
    print("\n🔍 LAYER 3 VERIFICATION")
    print("="*70)
    
    with driver.session() as session:
        # Nodi
        print("\n📊 NODES:")
        
        result = session.run("MATCH (r:RiskCategory) RETURN count(r) as c")
        print(f"   RiskCategory: {result.single()['c']:,}")
        
        result = session.run("MATCH (p:ProfessionalTag) RETURN count(p) as c")
        print(f"   ProfessionalTag: {result.single()['c']:,}")
        
        result = session.run("MATCH (t:Topic) RETURN count(t) as c")
        print(f"   Topic: {result.single()['c']:,}")
        
        # Relazioni
        print("\n🔗 RELATIONSHIPS:")
        
        result = session.run("MATCH ()-[r:ADDRESSES]->() RETURN count(r) as c")
        print(f"   ADDRESSES (Doc→Risk): {result.single()['c']:,}")
        
        result = session.run("MATCH ()-[r:TARGETS]->() RETURN count(r) as c")
        print(f"   TARGETS (Doc→Prof): {result.single()['c']:,}")
        
        result = session.run("MATCH ()-[r:DISCUSSES]->() RETURN count(r) as c")
        print(f"   DISCUSSES (Doc→Topic): {result.single()['c']:,}")
        
        result = session.run("MATCH ()-[r:GOVERNED_BY]->() RETURN count(r) as c")
        print(f"   GOVERNED_BY (Risk→Law): {result.single()['c']:,}")
        
        result = session.run("MATCH ()-[r:REGULATED_BY]->() RETURN count(r) as c")
        print(f"   REGULATED_BY (Prof→Law): {result.single()['c']:,}")
        
        # Top categories
        print("\n🏆 TOP RISK CATEGORIES:")
        result = session.run("""
            MATCH (r:RiskCategory)
            RETURN r.name as name, r.frequency as freq, r.severity_level as severity
            ORDER BY r.frequency DESC
            LIMIT 10
        """)
        for rec in result:
            print(f"   {rec['name']}: {rec['freq']} ({rec['severity']})")
        
        # Top professions
        print("\n👷 TOP PROFESSIONS:")
        result = session.run("""
            MATCH (p:ProfessionalTag)
            RETURN p.categoria as name, p.frequency as freq, p.related_category as cat
            ORDER BY p.frequency DESC
            LIMIT 10
        """)
        for rec in result:
            print(f"   {rec['name']}: {rec['freq']} → {rec['cat']}")
        
        # Sample connections
        print("\n📋 SAMPLE CONNECTIONS:")
        result = session.run("""
            MATCH (doc:TechnicalDocument)-[:ADDRESSES]->(r:RiskCategory)
            RETURN doc.doc_id as doc, r.name as risk
            LIMIT 5
        """)
        for rec in result:
            print(f"   {rec['doc'][:40]}... → {rec['risk']}")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Layer 3 v2: Semantic Bridge")
    parser.add_argument('--source', choices=['INAIL', 'OSHA', 'BOTH'], 
                       default='BOTH', help='Data source')
    parser.add_argument('--clean', action='store_true',
                       help='Clean Layer 3 data before running')
    parser.add_argument('--verify', action='store_true',
                       help='Verify Layer 3 status')
    parser.add_argument('--skip-schema', action='store_true',
                       help='Skip schema setup')
    
    args = parser.parse_args()
    
    # Solo verify
    if args.verify and not args.clean and args.source == 'BOTH':
        verify_layer3()
        driver.close()
        print("\n✅ Done!")
        sys.exit(0)
    
    # Clean se richiesto
    if args.clean:
        clean_layer3_only()
    
    # Schema
    if not args.skip_schema:
        ensure_schema()
    
    # Extract
    if args.source in ['INAIL', 'BOTH']:
        extract_semantic_layer('INAIL')
    
    if args.source in ['OSHA', 'BOTH']:
        extract_semantic_layer('OSHA')
    
    # Verify finale
    verify_layer3()
    
    driver.close()
    print("\n✅ Done!")
