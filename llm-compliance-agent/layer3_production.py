#!/usr/bin/env python3

"""

✅ LAYER 3: SEMANTIC BRIDGE - PRODUCTION VERSION



Features:

1. ✅ Usa doc_id_utils per allineamento perfetto

2. ✅ Entity extraction robusta

3. ✅ Gestione graceful di documenti mancanti

4. ✅ Statistiche dettagliate

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

from doc_id_utils import clean_doc_id, extract_doc_id_from_metadata, validate_doc_id



# ======================== CONFIG ========================



NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "")



INAIL_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/INAILProj/enriched_json")

OSHA_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/OSHAProj/enriched_json")



BATCH_SIZE = 100



driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))





# ======================== SCHEMA ========================



def ensure_schema():

    """Setup schema constraints"""

    

    print("\n🔧 Ensuring schema constraints...")

    

    with driver.session() as session:

        # TechnicalDocument unique

        try:

            session.run("""

                CREATE CONSTRAINT technical_document_doc_id IF NOT EXISTS

                FOR (d:TechnicalDocument)

                REQUIRE d.doc_id IS UNIQUE

            """)

            print("   ✅ TechnicalDocument.doc_id UNIQUE")

        except Exception as e:

            if "already exists" in str(e).lower():

                print("   ⭐ TechnicalDocument.doc_id (exists)")

        

        # Layer 3 constraints

        constraints = [

            "CREATE CONSTRAINT RiskCategory_name IF NOT EXISTS FOR (r:RiskCategory) REQUIRE r.name IS UNIQUE",

            "CREATE CONSTRAINT ProfessionalTag_categoria IF NOT EXISTS FOR (p:ProfessionalTag) REQUIRE p.categoria IS UNIQUE",

            "CREATE CONSTRAINT Topic_name IF NOT EXISTS FOR (t:Topic) REQUIRE t.name IS UNIQUE",

        ]

        

        for constraint in constraints:

            try:

                session.run(constraint)

                name = constraint.split()[2]

                print(f"   ✅ {name}")

            except Exception as e:

                if "already exists" in str(e).lower():

                    name = constraint.split()[2]

                    print(f"   ⭐ {name} (exists)")

    

    print()





# ======================== ENTITY EXTRACTOR ========================



class SemanticEntityExtractor:

    """Estrattore di entità semantiche da enriched JSON"""

    

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

        ✅ FIXED: Usa doc_id_utils per doc_id canonico

        """

        

        try:

            with open(json_path, 'r', encoding='utf-8') as f:

                data = json.load(f)

        except Exception as e:

            self.error_count += 1

            return None

        

        web_meta = data.get('web_metadata', {})

        sem_meta = data.get('semantic_metadata', {})

        

        # ✅ USA FUNZIONE CENTRALIZZATA

        doc_id = extract_doc_id_from_metadata(web_meta, json_path)

        

        # ✅ VALIDAZIONE

        if not validate_doc_id(doc_id):

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

        # 1. RISK CATEGORIES

        # ============================================

        

        main_cat = sem_meta.get('categoria_principale', '')

        if main_cat and len(main_cat) > 3:

            entities['risk_categories'].append({

                'name': main_cat.strip(),

                'type': 'main',

                'severity': self._infer_severity(main_cat)

            })

            self.risk_categories[main_cat.strip()] += 1

        

        for tema in sem_meta.get('temi_principali', [])[:5]:

            if tema and len(tema) > 3:

                entities['risk_categories'].append({

                    'name': tema.strip(),

                    'type': 'theme',

                    'severity': self._infer_severity(tema)

                })

                self.risk_categories[tema.strip()] += 1

        

        # ============================================

        # 2. PROFESSIONAL TAGS

        # ============================================

        

        for prof in sem_meta.get('categorie_professionali', [])[:10]:

            if prof and len(prof) > 3:

                entities['professional_tags'].append({

                    'categoria': prof.strip(),

                    'sector': self._infer_sector(prof),

                    'source': source

                })

                self.professional_tags[prof.strip()] += 1

        

        # ============================================

        # 3. TOPICS

        # ============================================

        

        for kw in sem_meta.get('parole_chiave', [])[:15]:

            if kw and len(kw) > 2:

                entities['topics'].append({

                    'name': kw.lower().strip(),

                    'source': source

                })

                self.topics[kw.lower().strip()] += 1

        

        # ============================================

        # 4. CITED LAWS

        # ============================================

        

        laws_structured = sem_meta.get('laws_structured', [])

        

        for law in laws_structured[:20]:

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

    

    def _infer_severity(self, risk_name: str) -> str:

        """Infer risk severity from name"""

        high_keywords = ['morte', 'grave', 'caduta', 'esplosion', 'incendi', 

                         'elettric', 'crollo', 'infortuni gravi', 'tossic']

        medium_keywords = ['lesion', 'danno', 'esposizion', 'stress', 'disturb']

        

        risk_lower = risk_name.lower()

        

        if any(kw in risk_lower for kw in high_keywords):

            return 'HIGH'

        elif any(kw in risk_lower for kw in medium_keywords):

            return 'MEDIUM'

        return 'LOW'

    

    def _infer_sector(self, profession: str) -> str:

        """Infer professional sector"""

        prof_lower = profession.lower()

        

        sectors = [

            (['edil', 'costruz', 'muratori', 'cantier'], 'Construction'),

            (['industri', 'manifattur', 'operai', 'fabbric'], 'Manufacturing'),

            (['sanitari', 'medic', 'infermier', 'ospedale'], 'Healthcare'),

            (['agricol', 'contadin', 'allevator'], 'Agriculture'),

            (['trasport', 'logistic', 'autisti', 'magazzin'], 'Transportation'),

            (['ufficio', 'impiegat', 'amministrat'], 'Office'),

        ]

        

        for keywords, sector in sectors:

            if any(kw in prof_lower for kw in keywords):

                return sector

        

        return 'General'

    

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

            r.frequency = 1,

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

            p.sector = entity.sector,

            p.source = entity.source,

            p.frequency = 1,

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

            t.created_at = datetime()

        ON MATCH SET

            t.frequency = t.frequency + 1

        RETURN count(t) as cnt

        """

    

    result = session.run(query, entities=entities_batch)

    return result.single()['cnt']





def link_document_to_entities(session, doc_id, risks, professions, topics):

    """

    ✅ FIXED: EXACT doc_id matching + graceful handling

    

    Returns:

        (success_count, is_matched)

    """

    

    if not risks and not professions and not topics:

        return 0, False

    

    # Prima verifica che il documento esista

    check_query = """

    MATCH (doc:TechnicalDocument {doc_id: $doc_id})

    RETURN doc.doc_id as id

    """

    

    result = session.run(check_query, doc_id=doc_id)

    if not result.single():

        return 0, False

    

    # Se esiste, crea i link

    query = """

    MATCH (doc:TechnicalDocument {doc_id: $doc_id})

    

    // Link Risks

    WITH doc

    UNWIND $risks as risk

    MATCH (r:RiskCategory {name: risk.name})

    MERGE (doc)-[rel:ADDRESSES]->(r)

    ON CREATE SET rel.type = risk.type, rel.created_at = datetime()

    

    // Link Professions

    WITH doc

    UNWIND $professions as prof

    MATCH (p:ProfessionalTag {categoria: prof.categoria})

    MERGE (doc)-[:TARGETS]->(p)

    

    // Link Topics

    WITH doc

    UNWIND $topics as topic

    MATCH (t:Topic {name: topic.name})

    MERGE (doc)-[:DISCUSSES]->(t)

    

    RETURN count(*) as cnt

    """

    

    try:

        rec = session.run(

            query,

            doc_id=doc_id,

            risks=risks or [],

            professions=professions or [],

            topics=topics or []

        ).single()

        

        return rec['cnt'] if rec else 0, True

    except Exception as e:

        return 0, False





def create_entity_law_bridges(session, law_to_entities):

    """Crea bridge tra entità semantiche e leggi"""

    

    if not law_to_entities:

        return 0, 0

    

    # Risk → Law

    query_risk = """

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

    """

    

    # Profession → Law

    query_prof = """

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

    """

    

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

    

    # Execute

    try:

        result_risk = session.run(query_risk, mappings=mappings)

        cnt_risk = result_risk.single()['cnt']

    except Exception as e:

        print(f"  ⚠️ Risk bridge error: {str(e)[:100]}")

        cnt_risk = 0

    

    try:

        result_prof = session.run(query_prof, mappings=mappings)

        cnt_prof = result_prof.single()['cnt']

    except Exception as e:

        print(f"  ⚠️ Prof bridge error: {str(e)[:100]}")

        cnt_prof = 0

    

    return cnt_risk, cnt_prof





# ======================== MAIN PIPELINE ========================



def extract_semantic_layer(source='INAIL'):

    """Main extraction pipeline"""

    

    print(f"\n🚀 LAYER 3: SEMANTIC EXTRACTION ({source})")

    print("="*70)

    

    data_dir = INAIL_ENRICHED if source == 'INAIL' else OSHA_ENRICHED

    json_files = sorted(list(data_dir.glob("*.json")))

    

    print(f"\n📂 Found {len(json_files):,} JSON files")

    

    if len(json_files) == 0:

        print(f"❌ No files in {data_dir}")

        return

    

    # Extract entities

    print(f"\n📄 Step 1/4: Extracting semantic entities...\n")

    

    extractor = SemanticEntityExtractor()

    all_entities = []

    

    for json_file in tqdm(json_files, desc="Parsing"):

        entities = extractor.extract_from_json(json_file, source)

        if entities:

            all_entities.append(entities)

    

    extractor.print_stats()

    

    print(f"\n✅ Extracted entities from {len(all_entities):,} documents")

    print(f"   Unique risks: {len(extractor.risk_categories)}")

    print(f"   Unique professions: {len(extractor.professional_tags)}")

    print(f"   Unique topics: {len(extractor.topics)}")

    

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

        print(f"      (doc_id exists in JSON but NOT in Neo4j TechnicalDocument)")

        print(f"\n   📋 Sample missing doc_ids (first 10):")

        for doc_id in missing_samples:

            print(f"      - {doc_id}")

        print(f"\n   💡 Run: python3 backfill_missing_docs.py")

    

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

    print(f"✅ LAYER 3 COMPLETE ({source})")

    print(f"{'='*70}\n")





if __name__ == "__main__":

    

    parser = argparse.ArgumentParser()

    parser.add_argument('--source', choices=['INAIL', 'OSHA', 'BOTH'], 

                       default='INAIL', help='Data source')

    parser.add_argument('--skip-schema', action='store_true',

                       help='Skip schema setup')

    

    args = parser.parse_args()

    

    if not args.skip_schema:

        ensure_schema()

    

    if args.source in ['INAIL', 'BOTH']:

        extract_semantic_layer('INAIL')

    

    if args.source in ['OSHA', 'BOTH']:

        extract_semantic_layer('OSHA')

    

    driver.close()

    print("\n✅ Done!\n")
