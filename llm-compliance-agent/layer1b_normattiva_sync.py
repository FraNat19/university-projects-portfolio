# -*- coding: utf-8 -*-

#!/usr/bin/env python3

"""

LAYER 1B: Normattiva → Neo4j Enrichment - FIXED



✅ Usa law_id come primary key

✅ Arricchisce nodi esistenti da Layer 0b/1

✅ Usa EnhancedLawExtractor per consistenza

✅ Converte date ISO correttamente

"""



from neo4j import GraphDatabase

from pathlib import Path

import xml.etree.ElementTree as ET

from tqdm import tqdm

import argparse

import re

from typing import Dict, List, Optional

import sys

import os



# ✅ IMPORT EXTRACTOR

sys.path.insert(0, os.path.dirname(__file__))

from law_extraction_enhanced import EnhancedLawExtractor



# ======================== CONFIG ========================

NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "thesis2024")

NEO4J_AUTH = (NEO4J_USER, NEO4J_PASSWORD)



NORMATTIVA_DIR = Path("/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/normattiva_filtered")



BATCH_SIZE = 500



# XML Namespaces

NS = {

    'nir': 'http://www.normeinrete.it/nir/2.2/',

    'dsp': 'http://www.normeinrete.it/nir/disposizioni/2.2/',

    'h': 'http://www.w3.org/HTML/1998/html4',

    'xlink': 'http://www.w3.org/1999/xlink'

}



# ======================== CLIENTS ========================

neo4j_driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)





# ======================== HELPER FUNCTIONS ========================



def extract_text_from_element(elem) -> str:

    """Extract all text from XML element recursively"""

    if elem is None:

        return ""



    text_parts = []

    if elem.text:

        text_parts.append(elem.text.strip())



    for child in elem:

        text_parts.append(extract_text_from_element(child))

        if child.tail:

            text_parts.append(child.tail.strip())



    return " ".join([t for t in text_parts if t])





def extract_vigenza(xml_path: Path, root) -> tuple:

    """

    Estrai date di vigenza da filename XML

    Format: YYYY-MM-DD_NNNNUXXXX_VIGENZA_YYYY-MM-DD_VN.xml

    """

    filename = xml_path.name.upper()



    # Extract snapshot date (quella dopo VIGENZA_)

    match = re.search(r'VIGENZA_(\d{4}-\d{2}-\d{2})', filename)

    vigenza_end = match.group(1) if match else None



    # ✅ vigenza_start dal dataDoc norm=YYYYMMDD (convertito in ISO)

    data_doc_elem = root.find('.//nir:dataDoc', NS)

    vigenza_start = None

    if data_doc_elem is not None:

        norm_date = data_doc_elem.get('norm', '')

        if len(norm_date) == 8 and norm_date.isdigit():

            vigenza_start = f"{norm_date[:4]}-{norm_date[4:6]}-{norm_date[6:8]}"



    # Questo dataset è snapshot "in vigenza" se il filename contiene VIGENZA

    is_vigente = "VIGENZA_" in filename



    return vigenza_start, vigenza_end, is_vigente





def parse_normattiva_xml(xml_path: Path, extractor: EnhancedLawExtractor) -> Optional[Dict]:

    """

    Parse completo XML Normattiva

    ✅ FIXED: Usa EnhancedLawExtractor per law_id canonico

    """

    try:

        tree = ET.parse(xml_path)

        root = tree.getroot()



        tipo_doc_elem = root.find('.//nir:tipoDoc', NS)

        num_doc_elem = root.find('.//nir:numDoc', NS)

        data_doc_elem = root.find('.//nir:dataDoc', NS)

        titolo_doc_elem = root.find('.//nir:titoloDoc', NS)

        urn_elem = root.find('.//nir:urn', NS)



        tipo_text = tipo_doc_elem.text if tipo_doc_elem is not None else ""

        numero = num_doc_elem.text if num_doc_elem is not None else ""

        titolo = extract_text_from_element(titolo_doc_elem) if titolo_doc_elem is not None else ""



        # ✅ ESTRAI DATA REALE (norm="YYYYMMDD" → ISO)

        data_norm = data_doc_elem.get('norm', '') if data_doc_elem is not None else ""

        anno = data_norm[:4] if len(data_norm) >= 4 else ""



        # Converti in ISO

        data_iso = None

        if len(data_norm) == 8:

            try:

                data_iso = f"{data_norm[:4]}-{data_norm[4:6]}-{data_norm[6:8]}"

            except:

                pass



        # ✅ USA EXTRACTOR per law_id canonico

        synthetic_citation = f"{tipo_text} {numero}/{anno}"

        laws, conf = extractor.extract_with_confidence(synthetic_citation, max_laws=1)



        if laws:

            law = laws[0]

            law_id = law['law_id']

            tipo_norm = law['tipo']

            tipo_short = law['tipo_short']

        else:

            # Fallback manuale

            folder_path = str(xml_path.parent.parent.name).lower()

            if 'decreti legislativi' in folder_path:

                tipo_short = 'dlgs'

                tipo_norm = 'decreto.legislativo'

            elif 'dpr' in folder_path:

                tipo_short = 'dpr'

                tipo_norm = 'decreto.del.presidente.della.repubblica'

            elif 'leggi' in folder_path:

                tipo_short = 'legge'

                tipo_norm = 'legge'

            else:

                tipo_short = 'altro'

                tipo_norm = 'altro'



            law_id = f"{tipo_short}-{numero}-{anno}"



        # URN

        urn_raw = urn_elem.get('valore', '') if urn_elem is not None else ''

        if urn_raw and urn_raw != 'urn:' and len(urn_raw) > 10:

            urn = urn_raw

        else:

            urn = f"urn:nir:stato:{tipo_norm}:{data_iso or anno+'-01-01'};{numero}"



        vigenza_start, vigenza_end, is_vigente = extract_vigenza(xml_path, root)



        # Extract articles

        articoli = []

        for articolo_elem in root.findall('.//nir:articolo', NS):

            art_id = articolo_elem.get('id', '')

            num_elem = articolo_elem.find('.//nir:num', NS)

            art_num_text = extract_text_from_element(num_elem) if num_elem is not None else ""

            art_num_clean = re.sub(r'[^\d]', '', art_num_text)

            art_text = extract_text_from_element(articolo_elem)



            if art_text and len(art_text) > 50 and art_num_clean:

                articoli.append({

                    'article_id': f"{law_id}_art{art_num_clean}",

                    'article_num': art_num_clean,

                    'article_num_display': art_num_text,

                    'text': art_text[:5000],

                    'full_text': art_text

                })



        # Pubblicazione

        pubblicazione_elem = root.find('.//nir:pubblicazione', NS)

        pub_data_raw = pubblicazione_elem.get('norm', '') if pubblicazione_elem is not None else ''

        pub_num = pubblicazione_elem.get('num', '') if pubblicazione_elem is not None else ''



        # ✅ Converti pub_data in ISO

        pub_data_iso = None

        if len(pub_data_raw) == 8:

            try:

                pub_data_iso = f"{pub_data_raw[:4]}-{pub_data_raw[4:6]}-{pub_data_raw[6:8]}"

            except:

                pass

        safe_year = anno if (isinstance(anno, str) and len(anno) == 4 and anno.isdigit()) else None

        safe_fallback_date = f"{safe_year}-01-01" if safe_year else None  

        return {

            'law_id': law_id,  # ✅ CANONICO

            'tipo': tipo_norm,

            'tipo_short': tipo_short,

            'numero': numero,

            'anno': anno,

            'data': data_iso or safe_fallback_date,

            'titolo': titolo[:1000],

            'urn': urn,

            'vigenza_start': vigenza_start,

            'vigenza_end': vigenza_end,

            'is_vigente': is_vigente,

            'pubblicazione_data': pub_data_iso or pub_data_raw,

            'pubblicazione_num': pub_num,

            'num_articles': len(articoli),

            'articoli': articoli,

            'xml_path': str(xml_path)

        }



    except Exception as e:

        return None





# ======================== NEO4J OPS ========================



def enrich_canonical_law_batch(session, laws_batch):

    """

    ✅ ARRICCHISCE CanonicalLaw usando law_id

    """

    query = """

    UNWIND $laws as law

    WITH law

    WHERE law.law_id IS NOT NULL AND trim(law.law_id) <> ''



    // ✅ Match by law_id (primary key)

    MATCH (cl:CanonicalLaw {law_id: law.law_id})



    // Enrich (non sovrascrive dati già presenti)

    SET

      cl.urn = COALESCE(cl.urn, law.urn),

      cl.tipo = COALESCE(cl.tipo, law.tipo),

      cl.tipo_short = COALESCE(cl.tipo_short, law.tipo_short),

      cl.numero = COALESCE(cl.numero, law.numero),

      cl.anno = COALESCE(cl.anno, law.anno),

      cl.titolo = COALESCE(cl.titolo, law.titolo),

      cl.data = COALESCE(cl.data, CASE

        WHEN law.data IS NOT NULL

        THEN date(law.data)

        ELSE NULL

      END),

      cl.vigenza_start = CASE

        WHEN law.vigenza_start IS NOT NULL

        THEN date(law.vigenza_start)

        ELSE cl.vigenza_start

      END,

      cl.vigenza_end = CASE

        WHEN law.vigenza_end IS NOT NULL

        THEN date(law.vigenza_end)

        ELSE cl.vigenza_end

      END,

      cl.is_vigente = law.is_vigente,

      cl.pubblicazione_data = COALESCE(cl.pubblicazione_data, law.pubblicazione_data),

      cl.pubblicazione_num = COALESCE(cl.pubblicazione_num, law.pubblicazione_num),

      cl.num_articles = COALESCE(cl.num_articles, law.num_articles),

      cl.has_xml = true,

      cl.enriched = true,

      cl.enriched_at = datetime()



    RETURN count(cl) as enriched

    """

    result = session.run(query, laws=laws_batch)

    return result.single()["enriched"]





def enrich_article_nodes_batch(session, articles_batch):

    """

    ✅ ARRICCHISCE Article esistenti (NON crea nuovi!)

    """

    query = """

    UNWIND $articles as art

    WITH art

    WHERE art.parent_law_id IS NOT NULL AND trim(art.parent_law_id) <> ''



    // ✅ Match parent by law_id

    MATCH (law:CanonicalLaw {law_id: art.parent_law_id})



    // ✅ MATCH Article esistente (non MERGE!)

    MATCH (law)-[:HAS_ARTICLE]->(a:Article)

    WHERE a.article_num = art.article_num



    // Arricchisci con testo completo da XML

    SET a.text = COALESCE(art.text, a.text),

        a.full_text = COALESCE(art.full_text, a.full_text),

        a.article_num_display = COALESCE(art.article_num_display, a.article_num_display),

        a.from_xml = true,

        a.enriched_at = datetime()



    RETURN count(a) as enriched

    """

    result = session.run(query, articles=articles_batch)

    return result.single()["enriched"]





def sync_normattiva_to_neo4j():

    """Main sync"""

    print(f"\n📄 Starting Normattiva → Neo4j sync (LAYER 1B - FIXED)...")

    print(f"📂 Source: {NORMATTIVA_DIR}")



    # ✅ CREATE EXTRACTOR

    extractor = EnhancedLawExtractor()



    print(f"\n🔍 Scanning for XML files...\n")

    xml_files = []

    for top_folder in NORMATTIVA_DIR.iterdir():

        if top_folder.is_dir() and not top_folder.name.endswith('.json'):

            for law_folder in top_folder.iterdir():

                if law_folder.is_dir():

                    xml_files.extend(list(law_folder.glob("*.xml")))



    print(f"📄 Found {len(xml_files):,} XML files\n")



    if len(xml_files) == 0:

        print("❌ NO XML FILES FOUND!\n")

        return 0, 0



    print(f"📦 Parsing XMLs...\n")

    laws_data = []

    all_articles = []

    skipped_count = 0



    for xml_file in tqdm(xml_files, desc="Parsing XMLs"):

        parsed = parse_normattiva_xml(xml_file, extractor)

        if parsed:

            laws_data.append(parsed)



            # ✅ Usa law_id come parent_law_id

            for art in parsed['articoli']:

                all_articles.append({

                    'article_id': art['article_id'],

                    'article_num': art['article_num'],

                    'article_num_display': art['article_num_display'],

                    'text': art['text'],

                    'full_text': art['full_text'],

                    'parent_law_id': parsed['law_id'],  # ✅ law_id

                    'law_tipo': parsed['tipo'],

                    'law_numero': parsed['numero'],

                    'law_anno': parsed['anno']

                })

        else:

            skipped_count += 1



    print(f"\n✅ Parsed:")

    print(f"   Laws: {len(laws_data):,}")

    print(f"   Articles: {len(all_articles):,}")

    print(f"   Skipped: {skipped_count:,}\n")



    print(f"⬆️ Enriching CanonicalLaw nodes...\n")

    enriched_count = 0



    with neo4j_driver.session() as session:

        for i in tqdm(range(0, len(laws_data), BATCH_SIZE), desc="Enrich CanonicalLaw"):

            batch = laws_data[i:i+BATCH_SIZE]

            batch_data = []

            for law in batch:

                batch_data.append({

                    'law_id': law['law_id'],  # ✅ CANONICO

                    'tipo': law['tipo'],

                    'tipo_short': law['tipo_short'],

                    'numero': law['numero'],

                    'anno': law['anno'],

                    'data': law['data'],

                    'titolo': law['titolo'],

                    'urn': law['urn'],

                    'vigenza_start': law['vigenza_start'],

                    'vigenza_end': law['vigenza_end'],

                    'is_vigente': law['is_vigente'],

                    'pubblicazione_data': law['pubblicazione_data'],

                    'pubblicazione_num': law['pubblicazione_num'],

                    'num_articles': law['num_articles']

                })



            enriched = enrich_canonical_law_batch(session, batch_data)

            enriched_count += enriched



    print(f"\n✅ Enriched {enriched_count:,} CanonicalLaw nodes\n")



    print(f"⬆️ Enriching Article nodes (existing only)...\n")

    enriched_articles = 0



    with neo4j_driver.session() as session:

        for i in tqdm(range(0, len(all_articles), BATCH_SIZE), desc="Enrich Articles"):

            batch = all_articles[i:i+BATCH_SIZE]

            enriched = enrich_article_nodes_batch(session, batch)

            enriched_articles += enriched



    print(f"\n✅ Enriched {enriched_articles:,} Article nodes\n")



    return enriched_count, enriched_articles





def verify_enrichment():

    """Verifica arricchimento"""

    with neo4j_driver.session() as session:



        print("\n🔍 Verification...\n")



        # CanonicalLaw enriched

        result = session.run("""

            MATCH (cl:CanonicalLaw)

            RETURN

                count(cl) as total,

                count(CASE WHEN cl.enriched = true THEN 1 END) as enriched,

                count(CASE WHEN cl.is_vigente = true THEN 1 END) as vigenti,

                count(CASE WHEN cl.vigenza_start IS NOT NULL THEN 1 END) as with_vigenza

        """)

        stats = result.single()



        print(f"📊 CanonicalLaw:")

        print(f"  Total: {stats['total']:,}")

        print(f"  Enriched: {stats['enriched']:,}")

        print(f"  Vigenti: {stats['vigenti']:,}")

        print(f"  With vigenza dates: {stats['with_vigenza']:,}")



        # Article nodes enriched

        result = session.run("""

            MATCH (a:Article)

            RETURN

                count(a) as total,

                count(CASE WHEN a.from_xml = true THEN 1 END) as enriched_from_xml

        """)

        art_stats = result.single()

        print(f"\n📄 Article nodes:")

        print(f"  Total: {art_stats['total']:,}")

        print(f"  Enriched from XML: {art_stats['enriched_from_xml']:,}")



        # HAS_ARTICLE relationships

        result = session.run("""

            MATCH ()-[r:HAS_ARTICLE]->()

            RETURN count(r) as cnt

        """)

        has_art_count = result.single()['cnt']

        print(f"\n🔗 HAS_ARTICLE relationships: {has_art_count:,}")



        # Sample

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WHERE law.enriched = true

            OPTIONAL MATCH (law)-[:HAS_ARTICLE]->(art:Article)

            RETURN law.law_id, law.titolo, law.is_vigente, count(art) as num_arts

            ORDER BY num_arts DESC

            LIMIT 5

        """)



        print(f"\n📋 Sample Enriched Laws:")

        for record in result:

            vigente = "✅" if record['law.is_vigente'] else "❌"

            print(f"  {record['law.law_id']}: {record['num_arts']} articles {vigente}")

            print(f"    {record['law.titolo'][:70]}...")





def check_layer2_readiness():

    """Verifica readiness Layer 2"""

    with neo4j_driver.session() as session:



        result = session.run("""

            MATCH (doc:TechnicalDocument)

            WITH count(doc) as docs

            MATCH (law:CanonicalLaw)

            WHERE law.enriched = true

            WITH docs, count(law) as laws

            MATCH (art:Article)

            WITH docs, laws, count(art) as arts

            MATCH ()-[c:CITES]->()

            WITH docs, laws, arts, count(c) as cites

            MATCH ()-[ha:HAS_ARTICLE]->()

            RETURN docs, laws, arts, cites, count(ha) as has_arts

        """)



        stats = result.single()



        print("\n" + "="*60)

        print("🎯 LAYER 2 READINESS CHECK")

        print("="*60)



        ready = (

            stats['docs'] > 100 and

            stats['laws'] > 1000 and

            stats['arts'] > 10000

        )



        if ready:

            print(f"\n✅ READY FOR LAYER 2!")

            print(f"   - TechnicalDocument: {stats['docs']:,}")

            print(f"   - CanonicalLaw (enriched): {stats['laws']:,}")

            print(f"   - Article: {stats['arts']:,}")

            print(f"   - CITES: {stats['cites']:,}")

            print(f"   - HAS_ARTICLE: {stats['has_arts']:,}")

            print(f"\nNext step:")

            print(f"  python3 layer2_extract_relations.py --sync")

        else:

            print(f"\n⚠️ NOT READY!")

            print(f"   - TechnicalDocument: {stats['docs']:,}")

            print(f"   - CanonicalLaw: {stats['laws']:,}")

            print(f"   - Article: {stats['arts']:,} (need >10K)")

            print(f"\n💡 Run: python3 layer1b_normattiva_sync.py --sync")



        print("\n" + "="*60)





if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--sync', action='store_true')

    parser.add_argument('--verify', action='store_true')

    parser.add_argument('--check', action='store_true')



    args = parser.parse_args()



    if not any([args.sync, args.verify, args.check]):

        print("\nUsage:")

        print("  --sync    : Sync Normattiva → Neo4j")

        print("  --verify  : Verify enrichment")

        print("  --check   : Check Layer 2 readiness\n")

    else:

        if args.sync:

            enriched, articles = sync_normattiva_to_neo4j()



        if args.verify:

            verify_enrichment()



        if args.check:

            check_layer2_readiness()



        neo4j_driver.close()

        print("\n✅ Done!")
