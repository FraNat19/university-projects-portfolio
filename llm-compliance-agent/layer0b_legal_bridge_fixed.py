#!/usr/bin/env python3



"""

✅ LAYER 0b FIXED: Sync con laws_structured

FIXES:

1. ✅ Sync laws_structured da technical_documents_text

2. ✅ Gestione date reali vs fallback (Fix errore date < 1000)

3. ✅ MERGE (non sovrascrive Normattiva)

"""



from qdrant_client import QdrantClient

from neo4j import GraphDatabase

import re

from datetime import datetime

import argparse

import os



QDRANT_HOST = os.getenv("QDRANT_HOST", "127.0.0.1")

QDRANT_PORT = int(os.getenv("QDRANT_PORT", "6333"))



NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "")



qdrant = QdrantClient(host=QDRANT_HOST, port=QDRANT_PORT, check_compatibility=False)

neo4j_driver = GraphDatabase.driver(

    NEO4J_URI, 

    auth=(NEO4J_USER, NEO4J_PASSWORD),

    max_connection_pool_size=10,

    connection_timeout=30.0

)



def parse_italian_date(date_str: str) -> str:

    """Converte date italiane in ISO"""

    

    if not date_str or not isinstance(date_str, str) or date_str.strip() == '':

        return None

    

    date_str = date_str.strip()

    

    # Già ISO

    if re.match(r'^\d{4}-\d{2}-\d{2}$', date_str):

        try:

            datetime.strptime(date_str, '%Y-%m-%d')

            return date_str

        except ValueError:

            return None

    

    # DD-MM-YYYY o DD/MM/YYYY

    match = re.match(r'^(\d{1,2})[-/](\d{1,2})[-/](\d{4})$', date_str)

    if match:

        giorno = match.group(1).zfill(2)

        mese = match.group(2).zfill(2)

        anno = match.group(3)

        iso_date = f"{anno}-{mese}-{giorno}"

        try:

            datetime.strptime(iso_date, '%Y-%m-%d')

            return iso_date

        except ValueError:

            return None

    

    # Italiano

    mesi_ita = {

        'gennaio': '01', 'febbraio': '02', 'marzo': '03',

        'aprile': '04', 'maggio': '05', 'giugno': '06',

        'luglio': '07', 'agosto': '08', 'settembre': '09',

        'ottobre': '10', 'novembre': '11', 'dicembre': '12'

    }

    

    match = re.match(r'(\d{1,2})\s+(\w+)\s+(\d{4})', date_str)

    if match:

        giorno = match.group(1).zfill(2)

        mese_ita = match.group(2).lower()

        anno = match.group(3)

        

        mese = mesi_ita.get(mese_ita)

        if mese:

            iso_date = f"{anno}-{mese}-{giorno}"

            try:

                datetime.strptime(iso_date, '%Y-%m-%d')

                return iso_date

            except ValueError:

                return None

    

    return None





def sync_legal_context():

    """Sync legal_context → CanonicalLaw (PRIORITÀ)"""

    

    print("\n📚 Syncing legal_context → CanonicalLaw (PRIORITY)...\n")

    

    offset = None

    created = 0

    

    with neo4j_driver.session() as session:

        

        batch_num = 0

        

        while True:

            scroll_result = qdrant.scroll(

                collection_name="legal_context",

                limit=100,

                offset=offset,

                with_payload=True,

                with_vectors=False

            )

            

            points, next_offset = scroll_result

            if not points:

                break

            

            batch_num += 1

            print(f"📦 Batch {batch_num}: {len(points)} points")

            

            batch_data = []

            for point in points:

                p = point.payload

                

                law_id = p.get('law_id', '')

                if not law_id:

                    continue

                

                # ✅ Usa tipo_short da payload

                tipo_short = p.get('tipo_short', '')

                if not tipo_short:

                    # Fallback

                    law_type = p.get('law_type', '')

                    tipo_map = {

                        'Decreto Legislativo': 'dlgs',

                        'decreto.legislativo': 'dlgs',

                        'DPR': 'dpr',

                        'decreto.del.presidente.della.repubblica': 'dpr',

                        'Legge': 'legge',

                        'legge': 'legge',

                        'Decreto Legge': 'dl',

                        'decreto.legge': 'dl',

                        'Codice': 'codice',

                        'codice': 'codice',

                        'Testo Unico': 'tu',

                        'testo.unico': 'tu',

                        'Regio Decreto': 'rd',

                        'regio.decreto': 'rd'

                    }

                    tipo_short = tipo_map.get(law_type, law_type.lower().replace(' ', '_'))

                

                # ✅ Data ISO

                data_iso = parse_italian_date(p.get('data', ''))

                

                # ✅ URN clean

                urn_raw = p.get('urn', '').strip()

                urn_clean = urn_raw if urn_raw and urn_raw != 'urn:' and len(urn_raw) > 10 else None

                

                batch_data.append({

                    'law_id': law_id,

                    'urn': urn_clean,

                    'tipo': tipo_short,

                    'numero': p.get('numero', ''),

                    'anno': p.get('anno', ''),

                    'titolo': p.get('titolo')[:500] if p.get('titolo') else '',

                    'data': data_iso,

                })

            

            if not batch_data:

                offset = next_offset

                continue

            

            # ✅ SET prioritizza Normattiva (non sovrascrive se già esiste)

            result = session.run("""

                UNWIND $batch as law

                

                MERGE (cl:CanonicalLaw {law_id: law.law_id})

                ON CREATE SET

                    cl.created_from = 'qdrant_legal_context',

                    cl.created_at = datetime()

                

                SET

                    cl.urn = COALESCE(cl.urn, law.urn),

                    cl.tipo = law.tipo,

                    cl.numero = law.numero,

                    cl.anno = law.anno,

                    cl.titolo = COALESCE(law.titolo, cl.titolo),

                    cl.data = CASE 

                        WHEN law.data IS NOT NULL 

                        THEN date(law.data)

                        ELSE cl.data 

                    END,

                    cl.has_qdrant_legal = true,

                    cl.synced_at = datetime()

                

                RETURN count(cl) as cnt

            """, batch=batch_data)

            

            created += result.single()['cnt']

            offset = next_offset

            

            if offset is None:

                break

    

    print(f"\n✅ Synced {created:,} from legal_context\n")

    return created





def sync_legal_articles():

    """Sync legal_articles → Article"""

    

    print("📄 Syncing legal_articles → Article...\n")

    

    offset = None

    created = 0

    

    with neo4j_driver.session() as session:

        

        batch_num = 0

        

        while True:

            scroll_result = qdrant.scroll(

                collection_name="legal_articles",

                limit=100,

                offset=offset,

                with_payload=True,

                with_vectors=False

            )

            

            points, next_offset = scroll_result

            if not points:

                break

            

            batch_num += 1

            print(f"📦 Batch {batch_num}: {len(points)} points")

            

            batch_data = []

            for point in points:

                p = point.payload

                

                article_id = p.get('article_id')

                parent_law_id = p.get('parent_law_id')

                

                if not article_id or not parent_law_id:

                    continue

                

                article_num_raw = p.get('article_num', '')

                match = re.search(r'\d+', article_num_raw)

                article_num = match.group(0) if match else article_num_raw

                

                if not article_num:

                    continue

                

                batch_data.append({

                    'article_id': article_id,

                    'parent_law_id': parent_law_id,

                    'article_num': article_num,

                    'article_num_display': article_num_raw,

                    'text': p.get('article_text', '')[:2000],

                    'full_text': p.get('full_text', '')[:5000],

                })

            

            if not batch_data:

                offset = next_offset

                continue

            

            result = session.run("""

                UNWIND $batch as art

                

                MATCH (law:CanonicalLaw {law_id: art.parent_law_id})

                

                MERGE (a:Article {article_id: art.article_id})

                SET 

                    a.article_num = art.article_num,

                    a.article_num_display = art.article_num_display,

                    a.text = art.text,

                    a.full_text = art.full_text,

                    a.parent_law_id = art.parent_law_id,

                    a.from_qdrant = true,

                    a.chunk_type = 'content',

                    a.updated_at = datetime()

                

                MERGE (law)-[r:HAS_ARTICLE {layer: 0}]->(a)

                SET r.source = 'qdrant_legal_articles',

                    r.created_at = datetime()

                

                RETURN count(a) as cnt

            """, batch=batch_data)

            

            created += result.single()['cnt']

            offset = next_offset

            

            if offset is None:

                break

    

    print(f"\n✅ Synced {created:,} Articles\n")

    return created





def sync_technical_laws():

    """

    ✅ NUOVO: Sync laws_structured da technical_documents_text (CON FIX DATE)

    """

    

    print("\n📄 Syncing laws from technical_documents_text (FILL GAPS)...\n")

    

    offset = None

    created = 0

    

    with neo4j_driver.session() as session:

        

        batch_num = 0

        

        while True:

            scroll_result = qdrant.scroll(

                collection_name="technical_documents_text",

                limit=100,

                offset=offset,

                with_payload=True,

                with_vectors=False

            )

            

            points, next_offset = scroll_result

            if not points:

                break

            

            batch_num += 1

            

            batch_data = []

            

            for point in points:

                p = point.payload

                

                # ✅ ESTRAI laws_structured

                laws_structured = p.get('laws_structured', [])

                

                for law in laws_structured:

                    

                    law_id = law.get('law_id')

                    if not law_id:

                        continue

                    

                    # Parse data (potrebbe essere None)

                    data_iso = parse_italian_date(law.get('data', ''))

                    

                    # ✅ FIX 1: Validazione Anno Python-side

                    anno_raw = law.get('anno', '')

                    anno_clean = None

                    try:

                        if anno_raw:

                            y = int(str(anno_raw))

                            # Accettiamo solo anni ragionevoli per evitare "679"

                            if 1000 <= y <= 9999:

                                anno_clean = str(y)

                    except (ValueError, TypeError):

                        pass

                    

                    batch_data.append({

                        'law_id': law_id,

                        'urn': law.get('urn', ''),

                        'tipo': law.get('tipo_short', ''),

                        'numero': law.get('numero', ''),

                        'anno': anno_clean, # Passiamo l'anno validato

                        'titolo': '',  

                        'data': data_iso,

                        'data_available': law.get('data_available', False)

                    })

            

            if not batch_data:

                offset = next_offset

                continue

            

            # ✅ FIX 2: MERGE con controllo Regex e Range in Cypher per la massima sicurezza

            result = session.run("""

                UNWIND $batch as law

                

                MERGE (cl:CanonicalLaw {law_id: law.law_id})

                ON CREATE SET

                    cl.created_from = 'qdrant_technical_docs',

                    cl.created_at = datetime()

                

                SET

                    cl.urn = COALESCE(cl.urn, law.urn),

                    cl.tipo = COALESCE(cl.tipo, law.tipo),

                    cl.numero = COALESCE(cl.numero, law.numero),

                    cl.anno = COALESCE(cl.anno, law.anno),

                    cl.data = CASE 

                        WHEN cl.data IS NOT NULL THEN cl.data

                        WHEN law.data IS NOT NULL THEN date(law.data)

                        WHEN law.anno IS NOT NULL 

                             AND toString(law.anno) =~ '^[0-9]{4}$'

                             AND toInteger(law.anno) >= 1000 

                             AND toInteger(law.anno) <= 9999 

                             THEN date(toString(law.anno) + '-01-01')

                        ELSE cl.data 

                    END,

                    cl.has_qdrant_technical = true,

                    cl.synced_at = datetime()

                

                RETURN count(cl) as cnt

            """, batch=batch_data)

            

            created += result.single()['cnt']

            

            print(f"📦 Batch {batch_num}: processed {len(batch_data)} laws")

            

            offset = next_offset

            

            if offset is None:

                break

    

    print(f"\n✅ Synced {created:,} laws from technical docs\n")

    return created





def verify_sync():

    """Verifica"""

    

    with neo4j_driver.session() as session:

        

        print("\n" + "="*70)

        print("🔍 VERIFICATION")

        print("="*70 + "\n")

        

        result = session.run("""

            MATCH (law:CanonicalLaw)

            RETURN 

                count(law) as total,

                count(CASE WHEN law.has_qdrant_legal = true THEN 1 END) as from_normattiva,

                count(CASE WHEN law.has_qdrant_technical = true THEN 1 END) as from_technical,

                count(CASE WHEN law.data IS NOT NULL THEN 1 END) as with_date,

                count(CASE WHEN law.titolo IS NOT NULL AND law.titolo <> '' THEN 1 END) as with_title

        """)

        

        stats = result.single()

        print(f"📚 CanonicalLaw:")

        print(f"   Total: {stats['total']:,}")

        print(f"   From Normattiva: {stats['from_normattiva']:,}")

        print(f"   From Technical: {stats['from_technical']:,}")

        print(f"   With date: {stats['with_date']:,}")

        print(f"   With title: {stats['with_title']:,}\n")

        

        result = session.run("""

            MATCH (art:Article)

            RETURN count(art) as total

        """)

        

        print(f"📄 Article: {result.single()['total']:,}\n")

        

        result = session.run("""

            MATCH ()-[r:HAS_ARTICLE]->()

            RETURN count(r) as cnt

        """)

        

        print(f"🔗 HAS_ARTICLE: {result.single()['cnt']:,}\n")

        

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WHERE law.data IS NULL OR law.data = date('1900-01-01')

            RETURN count(law) as cnt

        """)

        

        missing_dates = result.single()['cnt']

        print(f"⚠️  Laws with missing/fallback dates: {missing_dates:,}\n")

        

        result = session.run("""

            MATCH (law:CanonicalLaw)-[:HAS_ARTICLE]->(art:Article)

            RETURN law.law_id, count(art) as cnt

            ORDER BY cnt DESC

            LIMIT 5

        """)

        

        print(f"📋 Top laws by articles:")

        for rec in result:

            print(f"   {rec['law.law_id']}: {rec['cnt']} articles")

        

        print()





if __name__ == "__main__":

    

    parser = argparse.ArgumentParser()

    parser.add_argument('--sync', action='store_true')

    parser.add_argument('--verify', action='store_true')

    

    args = parser.parse_args()

    

    if not any([args.sync, args.verify]):

        print("\nUsage:")

        print("  --sync    : Sync ALL collections")

        print("  --verify  : Verify sync\n")

    else:

        if args.sync:

            print("\n" + "="*70)

            print("🔗 LAYER 0b: Complete Sync (Normattiva + Technical)")

            print("="*70)

            

            # ORDINE CRITICO:

            # 1. Normattiva (PRIORITÀ per titoli/date)

            laws_norm = sync_legal_context()

            articles = sync_legal_articles()

            

            # 2. Technical docs (FILL GAPS)

            laws_tech = sync_technical_laws()

            

            print("="*70)

            print(f"✅ DONE!")

            print(f"   Normattiva laws: {laws_norm:,}")

            print(f"   Articles: {articles:,}")

            print(f"   Technical laws: {laws_tech:,}")

            print("="*70 + "\n")

        

        if args.verify:

            verify_sync()

    

    neo4j_driver.close()
