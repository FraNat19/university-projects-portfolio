#!/usr/bin/env python3

"""

✅ LAYER 1 v3: CREA Article al volo se mancanti



FIXES:

1. ✅ Crea Article mancanti durante CITES

2. ✅ Gestione robusta fallback

3. ✅ Logging dettagliato

"""



from qdrant_client import QdrantClient

from neo4j import GraphDatabase

from tqdm import tqdm

import argparse



import os

from datetime import datetime  



QDRANT_HOST = os.getenv("QDRANT_HOST", "127.0.0.1")

QDRANT_PORT = int(os.getenv("QDRANT_PORT", "6333"))

QDRANT_COLLECTION = os.getenv("QDRANT_COLLECTION", "technical_documents_text")



NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "thesis2024")



BATCH_SIZE = int(os.getenv("BATCH_SIZE", "500"))



qdrant = QdrantClient(host=QDRANT_HOST, port=QDRANT_PORT, check_compatibility=False)

neo4j_driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))





def create_technical_document_nodes_batch(session, docs_batch):

    """Crea TechnicalDocument nodes"""

    

    query = """

    UNWIND $docs as doc

    

    MERGE (d:TechnicalDocument {doc_id: doc.doc_id})

    SET d.source = doc.source,

        d.title = doc.title,

        d.chunk_type = doc.chunk_type,

        d.chunk_index = doc.chunk_index,

        d.total_chunks = doc.total_chunks,

        d.data_pubblicazione = doc.data_pubblicazione,

        d.categoria_principale = doc.categoria_principale,

        d.sintesi = doc.sintesi,

        d.pdf_url = doc.pdf_url,

        d.has_law_citations = doc.has_law_citations,

        d.num_law_citations = doc.num_law_citations,

        d.confidence_score = doc.confidence_score,

        d.updated_at = datetime()

    

    RETURN count(d) as created

    """

    

    result = session.run(query, docs=docs_batch)

    return result.single()['created']





def create_canonical_law_stubs_batch(session, laws_batch):

    """Crea stub CanonicalLaw"""

    

    query = """

    UNWIND $laws as law

    

    MERGE (cl:CanonicalLaw {law_id: law.law_id})

    ON CREATE SET

        cl.created_from = 'layer1_stub',

        cl.created_at = datetime()

    

    SET

        cl.qdrant_prefix = law.law_id,

        cl.tipo = law.tipo,

        cl.tipo_short = law.tipo_short,

        cl.numero = law.numero,

        cl.anno = law.anno,

        cl.data = CASE

          WHEN law.data IS NOT NULL AND law.data <> '' THEN date(law.data)

          WHEN cl.data IS NULL

          AND law.anno IS NOT NULL

          AND toString(law.anno) =~ '^[0-9]{4}$'

          AND toInteger(law.anno) >= 1000 AND toInteger(law.anno) <= 9999

          THEN date(toString(law.anno) + '-01-01')

          ELSE cl.data

        END,

        cl.jurisdiction = law.jurisdiction,

        cl.urn = COALESCE(cl.urn, law.urn),

        cl.has_citations = true,

        cl.synced_at = datetime()

    

    RETURN count(cl) as created

    """

    

    result = session.run(query, laws=laws_batch)

    return result.single()['created']





def create_granular_cites_relationships(session, relations_batch):

    """

    ✅ CREA Article al volo se mancante

    """

    

    article_rels = [r for r in relations_batch if r.get('article_num')]

    law_rels = [r for r in relations_batch if not r.get('article_num')]

    

    cnt_article = 0

    cnt_law = 0

    cnt_article_created = 0

    

    # Query 1: CITES → Article (CREA se mancante)

    if article_rels:

        result = session.run("""

            UNWIND $relations as rel

            

            MATCH (doc:TechnicalDocument {doc_id: rel.doc_id})

            MATCH (law:CanonicalLaw {law_id: rel.law_id})

            

            // ✅ CREA Article se non esiste

            MERGE (law)-[:HAS_ARTICLE]->(art:Article {

                parent_law_id: rel.law_id,

                article_num: rel.article_num

            })

            ON CREATE SET

                art.article_id = rel.law_id + '_art' + rel.article_num,

                art.created_from = 'layer1_citation',

                art.created_at = datetime(),

                art.chunk_type = 'content'

            

            // CITES relationship

            MERGE (doc)-[c:CITES {layer: 1, granularity: 'article'}]->(art)

            SET c.chunk_index = rel.chunk_index,

                c.citation_context = rel.citation_context,

                c.cited_as = rel.cited_as,

                c.full_citation = rel.full_citation,

                c.created_at = datetime()

            

            RETURN 

                count(c) as cites_created,

                count(CASE WHEN art.created_from = 'layer1_citation' THEN 1 END) as articles_created

        """, relations=article_rels)

        

        rec = result.single()

        cnt_article = rec['cites_created']

        cnt_article_created = rec['articles_created']

    

    # Query 2: CITES → CanonicalLaw (citazione generica)

    if law_rels:

        result = session.run("""

            UNWIND $relations as rel

            

            MATCH (doc:TechnicalDocument {doc_id: rel.doc_id})

            MATCH (law:CanonicalLaw {law_id: rel.law_id})

            

            // Solo se NON esiste già CITES granulare

            WHERE NOT exists((doc)-[:CITES {granularity: 'article'}]->()-[:HAS_ARTICLE*0..1]-(law))

            

            MERGE (doc)-[c:CITES {layer: 1, granularity: 'law'}]->(law)

            SET c.chunk_index = rel.chunk_index,

                c.citation_context = rel.citation_context,

                c.cited_as = rel.cited_as,

                c.full_citation = rel.full_citation,

                c.created_at = datetime()

            

            RETURN count(c) as cnt

        """, relations=law_rels)

        

        cnt_law = result.single()['cnt']

    

    return cnt_article, cnt_law, cnt_article_created





def sync_qdrant_to_neo4j_granular():

    """Main sync con creazione Article dinamica"""

    

    print(f"\n📄 Starting GRANULAR sync (with Article creation)...")

    print(f"📊 Collection: {QDRANT_COLLECTION}\n")

    

    try:

        collection_info = qdrant.get_collection(QDRANT_COLLECTION)

        total_chunks = collection_info.points_count

        print(f"📂 Total chunks: {total_chunks:,}\n")

    except Exception as e:

        print(f"\n❌ ERROR: Collection not found!")

        return 0, 0, 0, 0

    

    batch_docs = []

    batch_laws = []

    batch_relations = []

    

    created_docs = 0

    created_laws = 0

    created_article_cites = 0

    created_law_cites = 0

    created_articles = 0

    

    with neo4j_driver.session() as session:

        

        offset = None

        processed_chunks = 0

        

        while True:

            

            scroll_result = qdrant.scroll(

                collection_name=QDRANT_COLLECTION,

                limit=1000,

                offset=offset,

                with_payload=True,

                with_vectors=False

            )

            

            points, next_offset = scroll_result

            

            if not points:

                break

            

            processed_chunks += len(points)

            print(f"📦 Batch: {len(points):,} chunks (total: {processed_chunks:,}/{total_chunks:,})")

            

            for point in points:

                

                payload = point.payload

                

                # Solo summary

                if payload.get('chunk_type') != 'summary':

                    continue

                

                doc_id = payload.get('doc_id', '')

                

                doc_data = {

                    'doc_id': doc_id,

                    'source': payload.get('source', ''),

                    'title': payload.get('title', ''),

                    'chunk_type': 'summary',

                    'chunk_index': 0,

                    'total_chunks': payload.get('total_chunks', 1),

                    'data_pubblicazione': payload.get('data_pubblicazione', ''),

                    'categoria_principale': payload.get('categoria_principale', ''),

                    'sintesi': payload.get('sintesi', ''),

                    'pdf_url': payload.get('pdf_url', ''),

                    'has_law_citations': payload.get('has_law_citations', False),

                    'num_law_citations': payload.get('num_law_citations', 0),

                    'confidence_score': payload.get('confidence_score', 0.0),

                }

                

                batch_docs.append(doc_data)

                

                # ✅ USA laws_structured

                laws_structured = payload.get('laws_structured', [])

                

                seen_laws = set()

                

                for law in laws_structured:

                    

                    law_id = law.get('law_id')

                    

                    if not law_id:

                        continue

                    

                    if law_id in seen_laws:

                        continue

                    seen_laws.add(law_id)

                    

                    # ✅ Law data (con data)

                    law_data = {

                        'law_id': law_id,

                        'tipo': law.get('tipo', 'unknown'),

                        'tipo_short': law.get('tipo_short', 'unknown'),

                        'numero': law.get('numero', ''),

                        'anno': law.get('anno', ''),

                        'data': law.get('data', ''),  # Può essere None

                        'urn': law.get('urn', ''),

                        'jurisdiction': law.get('jurisdiction', 'IT')

                    }

                    

                    batch_laws.append(law_data)

                    

                    # ✅ RELATIONS con granularità

                    article_num = law.get('article_num')

                    full_citation = law.get('full_citation', law_id)

                    

                    relation_data = {

                        'doc_id': doc_id,

                        'law_id': law_id,

                        'article_num': article_num,  # Può essere None

                        'chunk_index': 0,

                        'citation_context': 'document-level',

                        'cited_as': law.get('stringa_originale', ''),

                        'full_citation': full_citation

                    }

                    

                    batch_relations.append(relation_data)

            

            # Flush docs

            if len(batch_docs) >= BATCH_SIZE:

                cnt = create_technical_document_nodes_batch(session, batch_docs)

                created_docs += cnt

                batch_docs = []

            

            # Flush laws

            if len(batch_laws) >= BATCH_SIZE:

                cnt = create_canonical_law_stubs_batch(session, batch_laws)

                created_laws += cnt

                batch_laws = []

            

            # Flush relations

            if len(batch_relations) >= BATCH_SIZE:

                cnt_art, cnt_law, cnt_art_new = create_granular_cites_relationships(session, batch_relations)

                created_article_cites += cnt_art

                created_law_cites += cnt_law

                created_articles += cnt_art_new

                batch_relations = []

            

            offset = next_offset

            if offset is None:

                break

        

        # Final flush

        if batch_docs:

            cnt = create_technical_document_nodes_batch(session, batch_docs)

            created_docs += cnt

        

        if batch_laws:

            cnt = create_canonical_law_stubs_batch(session, batch_laws)

            created_laws += cnt

        

        if batch_relations:

            cnt_art, cnt_law, cnt_art_new = create_granular_cites_relationships(session, batch_relations)

            created_article_cites += cnt_art

            created_law_cites += cnt_law

            created_articles += cnt_art_new

    

    print(f"\n✅ Sync complete!")

    print(f"📊 Stats:")

    print(f"   TechnicalDocument nodes: {created_docs:,}")

    print(f"   CanonicalLaw stubs: {created_laws:,}")

    print(f"   CITES → Article: {created_article_cites:,}")

    print(f"   CITES → Law: {created_law_cites:,}")

    print(f"   Articles CREATED: {created_articles:,} ⭐\n")

    

    return created_docs, created_article_cites, created_law_cites, created_articles





def verify_sync():

    """Verifica completa"""

    

    with neo4j_driver.session() as session:

        

        print("\n🔍 Verification...\n")

        

        # CITES by granularity

        result = session.run("""

            MATCH ()-[c:CITES]->()

            RETURN 

                c.granularity as granularity,

                count(c) as cnt

            ORDER BY cnt DESC

        """)

        

        print("📊 CITES by granularity:")

        for rec in result:

            print(f"   {rec['granularity']}: {rec['cnt']:,}")

        

        # Articles created by layer1

        result = session.run("""

            MATCH (art:Article)

            WHERE art.created_from = 'layer1_citation'

            RETURN count(art) as cnt

        """)

        

        layer1_articles = result.single()['cnt']

        print(f"\n⭐ Articles created by layer1: {layer1_articles:,}")

        

        # Articles from Normattiva

        result = session.run("""

            MATCH (art:Article)

            WHERE art.from_qdrant = true

            RETURN count(art) as cnt

        """)

        

        normattiva_articles = result.single()['cnt']

        print(f"📚 Articles from Normattiva: {normattiva_articles:,}")

        

        # Total articles

        result = session.run("""

            MATCH (art:Article)

            RETURN count(art) as total

        """)

        

        total_articles = result.single()['total']

        print(f"📄 Total Articles: {total_articles:,}\n")

        

        # Check duplicati law_id

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WITH law.law_id as lid, count(*) as cnt

            WHERE cnt > 1

            RETURN count(*) as dup_groups

        """)

        

        dup_groups = result.single()['dup_groups']

        

        if dup_groups == 0:

            print(f"✅ No duplicate laws\n")

        else:

            print(f"⚠️  {dup_groups} duplicate law groups!\n")

        

        # Sample CITES

        result = session.run("""

            MATCH (doc:TechnicalDocument)-[c:CITES]->(target)

            RETURN 

                doc.doc_id as doc,

                labels(target)[0] as target_type,

                c.granularity as granularity,

                c.cited_as as citation

            LIMIT 5

        """)

        

        print("📋 Sample CITES:")

        for rec in result:

            print(f"   {rec['doc'][:30]}... → {rec['target_type']} ({rec['granularity']})")

            if rec['citation']:

                print(f"      Citation: {rec['citation'][:60]}...")

        

        print()





if __name__ == "__main__":

    

    parser = argparse.ArgumentParser()

    parser.add_argument('--sync', action='store_true', help='Sync with Article creation')

    parser.add_argument('--verify', action='store_true', help='Verify sync')

    

    args = parser.parse_args()

    

    if not any([args.sync, args.verify]):

        print("\nUsage:")

        print("  --sync    : Sync with automatic Article creation")

        print("  --verify  : Verify granularity and Article creation\n")

    else:

        if args.sync:

            sync_qdrant_to_neo4j_granular()

        

        if args.verify:

            verify_sync()

    

    neo4j_driver.close()

    print("\n✅ Done!")
