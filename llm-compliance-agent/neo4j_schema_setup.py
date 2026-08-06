#!/usr/bin/env python3

"""

Neo4j Schema Setup: Ontologia Legale + Constraints + Indexes

FIXED VERSION: NO APOC dependencies (universal compatibility)

"""



from neo4j import GraphDatabase

import sys



from neo4j import GraphDatabase

import sys

import os



NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "thesis2024")

NEO4J_AUTH = (NEO4J_USER, NEO4J_PASSWORD)



driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)



def reset_schema_only():

    """

    SAFE RESET: Cancella SOLO relazioni Layer 2+ e nodi derivati

    PRESERVA: Nodi TechnicalDocument da Qdrant

    """



    with driver.session() as session:



        print("\n🧹 SAFE RESET: Removing Layer 2+ data only...")



        # 1. Cancella relazioni Layer 2+

        result = session.run("""

            MATCH ()-[r]->()

            WHERE r.layer >= 2 OR r.layer IS NULL

            DELETE r

            RETURN count(r) as deleted

        """)

        deleted = result.single()['deleted']

        print(f"  ✅ Deleted {deleted:,} relationships")



        # 2. Cancella nodi CanonicalLaw (creati da Layer 2)

        result = session.run("""

            MATCH (c:CanonicalLaw)

            DETACH DELETE c

            RETURN count(c) as deleted

        """)

        deleted = result.single()['deleted']

        print(f"  ✅ Deleted {deleted:,} CanonicalLaw nodes")



        # 3. Cancella nodi Article (creati da Layer 2)

        result = session.run("""

            MATCH (a:Article)

            DETACH DELETE a

            RETURN count(a) as deleted

        """)

        deleted = result.single()['deleted']

        print(f"  ✅ Deleted {deleted:,} Article nodes")



        # 4. Verifica nodi TechnicalDocument preservati

        result = session.run("""

            MATCH (d:TechnicalDocument)

            RETURN count(d) as cnt

        """)

        cnt = result.single()['cnt']

        print(f"  ℹ️  Preserved {cnt:,} TechnicalDocument nodes from Qdrant")



        if cnt == 0:

            print("\n  ⚠️  WARNING: NO TechnicalDocument nodes found!")

            print("     Run Qdrant ingestion (02_text_chunks.py) first")



def setup_schema():

    """

    Setup completo:

    1. Constraints (unicità + esistenza)

    2. Indici (performance)

    3. Verifica nodi esistenti

    """



    with driver.session() as session:



        print("\n🔧 Creating Constraints...")



        constraints = [

            # Article: qdrant_id unico

            "CREATE CONSTRAINT Article_qdrant_id IF NOT EXISTS FOR (a:Article) REQUIRE a.qdrant_id IS UNIQUE",



            # CanonicalLaw: URN unico

            "CREATE CONSTRAINT CanonicalLaw_urn IF NOT EXISTS FOR (n:CanonicalLaw) REQUIRE n.urn IS UNIQUE",

            "CREATE CONSTRAINT CanonicalLaw_qdrant_prefix IF NOT EXISTS FOR (n:CanonicalLaw) REQUIRE n.qdrant_prefix IS UNIQUE",



            # TechnicalDocument: id unico

            "CREATE CONSTRAINT TechnicalDocument_id IF NOT EXISTS FOR (d:TechnicalDocument) REQUIRE d.doc_id IS UNIQUE",

        ]



        for constraint in constraints:

            try:

                session.run(constraint)

                constraint_name = constraint.split()[2]

                print(f"  ✅ {constraint_name}")

            except Exception as e:

                if "already exists" in str(e).lower() or "equivalent" in str(e).lower():

                    constraint_name = constraint.split()[2]

                    print(f"  ⏭️  {constraint_name} (già esiste)")

                else:

                    print(f"  ❌ {constraint.split()[2]}: {e}")



        print("\n📊 Creating Indexes...")



        indexes = [

            # Performance indexes

            "CREATE INDEX Article_law_lookup IF NOT EXISTS FOR (a:Article) ON (a.law_type, a.law_number, a.law_year)",

            "CREATE INDEX CanonicalLaw_tipo IF NOT EXISTS FOR (n:CanonicalLaw) ON (n.tipo)",

            "CREATE INDEX CanonicalLaw_year IF NOT EXISTS FOR (n:CanonicalLaw) ON (n.year)",

            "CREATE INDEX CanonicalLaw_vigente IF NOT EXISTS FOR (n:CanonicalLaw) ON (n.is_vigente)",

            "CREATE INDEX CanonicalLaw_vigenza IF NOT EXISTS FOR (n:CanonicalLaw) ON (n.vigenza_start, n.vigenza_end)",

            "CREATE INDEX TechnicalDocument_source IF NOT EXISTS FOR (d:TechnicalDocument) ON (d.source)",

        ]



        for index in indexes:

            try:

                session.run(index)

                index_name = index.split()[2]

                print(f"  ✅ {index_name}")

            except Exception as e:

                if "already exists" in str(e).lower() or "equivalent" in str(e).lower():

                    index_name = index.split()[2]

                    print(f"  ⏭️  {index_name} (già esiste)")

                else:

                    print(f"  ❌ {index.split()[2]}: {e}")



        print("\n🔍 Verifying Schema...")



        # ✅ NO APOC: Query diretta per contare nodi per label

        result = session.run("""

            MATCH (n)

            WITH labels(n)[0] as label, count(n) as count

            WHERE label IS NOT NULL

            RETURN label, count

            ORDER BY count DESC

        """)



        print("\n📈 Node Counts:")

        records = list(result)

        if records:

            for record in records:

                print(f"  {record['label']}: {record['count']:,}")

        else:

            print("  (no nodes yet)")



        # ✅ NO APOC: Query diretta per contare relazioni per tipo

        result = session.run("""

            MATCH ()-[r]->()

            WITH type(r) as relationshipType, count(r) as count

            RETURN relationshipType, count

            ORDER BY count DESC

        """)



        print("\n🔗 Relationship Counts:")

        records = list(result)

        if records:

            for record in records:

                print(f"  {record['relationshipType']}: {record['count']:,}")

        else:

            print("  (no relationships yet)")



        print("\n✅ Schema setup completato!")



def cleanup_isolated_nodes():

    """

    Rimuove nodi isolati (< 2 relazioni) per evitare esplosione grafo

    PRESERVA: CanonicalLaw, Article, TechnicalDocument

    """



    with driver.session() as session:



        print("\n🧹 Cleaning up isolated nodes...")



        # ✅ NO APOC: Query diretta con COUNT{}

        result = session.run("""

            MATCH (n)

            WHERE COUNT { (n)--() } < 2

            AND NOT n:CanonicalLaw

            AND NOT n:Article

            AND NOT n:TechnicalDocument

            RETURN count(n) as cnt

        """)



        isolated_count = result.single()['cnt']

        print(f"  📊 Nodi isolati trovati: {isolated_count:,}")



        if isolated_count == 0:

            print("  ✅ Nessun nodo isolato da rimuovere")

            return



        # ✅ NO APOC: Cancellazione batch manuale

        batch_size = 1000

        total_deleted = 0



        while True:

            result = session.run("""

                MATCH (n)

                WHERE COUNT { (n)--() } < 2

                AND NOT n:CanonicalLaw

                AND NOT n:Article

                AND NOT n:TechnicalDocument

                WITH n LIMIT $batch_size

                DETACH DELETE n

                RETURN count(n) as deleted

            """, batch_size=batch_size)



            deleted = result.single()['deleted']

            total_deleted += deleted



            if deleted == 0:

                break



            print(f"  🔄 Deleted batch: {deleted:,} nodes (total: {total_deleted:,})")



        print(f"  ✅ Rimossi {total_deleted:,} nodi isolati")



def merge_duplicates():

    """

    Identifica e marca duplicati basati su content similarity

    SAFE: Usa relazione DUPLICATE_OF invece di merge

    """



    with driver.session() as session:



        print("\n🔄 Marking duplicate articles...")



        # Trova duplicati (stesso content)

        result = session.run("""

            MATCH (a1:Article), (a2:Article)

            WHERE id(a1) < id(a2)

            AND a1.content = a2.content

            AND a1.qdrant_id STARTS WITH substring(a2.qdrant_id, 0, 10)

            RETURN count(*) as cnt

        """)



        dup_count = result.single()['cnt']

        print(f"  📊 Duplicati trovati: {dup_count:,}")



        if dup_count == 0:

            print("  ✅ Nessun duplicato da marcare")

            return



        # Marca duplicati con relazione

        result = session.run("""

            MATCH (a1:Article), (a2:Article)

            WHERE id(a1) < id(a2)

            AND a1.content = a2.content

            AND a1.qdrant_id STARTS WITH substring(a2.qdrant_id, 0, 10)



            MERGE (a1)-[:DUPLICATE_OF]->(a2)

            SET a2.is_duplicate = true



            RETURN count(*) as cnt

        """)



        dup_marked = result.single()['cnt']

        print(f"  ✅ Marcati {dup_marked:,} duplicati con relazione DUPLICATE_OF")



def print_graph_stats():

    """

    Stampa statistiche complete del grafo

    FIXED: Handle empty graph gracefully

    """



    with driver.session() as session:



        print("\n" + "="*70)

        print("📊 GRAPH STATISTICS")

        print("="*70)



        # Total nodes/relationships

        try:

            result = session.run("""

                MATCH (n)

                WITH count(n) as nodes

                OPTIONAL MATCH ()-[r]->()

                RETURN nodes, count(r) as rels

            """)



            stats = result.single()

            if stats and stats['nodes'] > 0:

                print(f"\nTotal Nodes: {stats['nodes']:,}")

                print(f"Total Relationships: {stats['rels']:,}")

            else:

                print("\n⚠️  No data in graph")

                print("   Run Qdrant ingestion (02_text_chunks.py) first")

                print("\n" + "="*70)

                return

        except Exception as e:

            print(f"\n❌ Error counting nodes: {e}")

            print("\n" + "="*70)

            return



        # Breakdown per label

        try:

            result = session.run("""

                MATCH (n)

                WITH labels(n)[0] as label, count(n) as count

                WHERE label IS NOT NULL

                RETURN label, count

                ORDER BY count DESC

            """)



            print("\n📦 Nodes by Type:")

            for record in result:

                print(f"  {record['label']}: {record['count']:,}")

        except Exception as e:

            print(f"\n⚠️  Node breakdown unavailable: {e}")



        # Leggi vigenti vs abrogate (se esistono)

        try:

            result = session.run("""

                MATCH (law:CanonicalLaw)

                RETURN

                    count(CASE WHEN law.is_vigente = true THEN 1 END) as vigenti,

                    count(CASE WHEN law.is_vigente = false THEN 1 END) as abrogate

            """)



            stats = result.single()

            if stats and (stats['vigenti'] > 0 or stats['abrogate'] > 0):

                print(f"\n📜 Leggi:")

                print(f"  Vigenti: {stats['vigenti']:,}")

                print(f"  Abrogate: {stats['abrogate']:,}")

        except Exception as e:

            pass  # CanonicalLaw non ancora popolato



        # TechnicalDocuments per source

        try:

            result = session.run("""

                MATCH (d:TechnicalDocument)

                WITH d.source as source, count(d) as count

                RETURN source, count

                ORDER BY count DESC

            """)



            records = list(result)

            if records:

                print(f"\n📄 Technical Documents:")

                for record in records:

                    source = record['source'] or 'unknown'

                    print(f"  {source}: {record['count']:,}")

        except Exception as e:

            pass  # TechnicalDocument non ancora popolato



        # Relazioni per confidence (se esistono)

        try:

            result = session.run("""

                MATCH ()-[r]->()

                WHERE r.confidence IS NOT NULL

                WITH

                    count(CASE WHEN r.confidence >= 0.9 THEN 1 END) as high_conf,

                    count(CASE WHEN r.confidence >= 0.7 AND r.confidence < 0.9 THEN 1 END) as med_conf,

                    count(CASE WHEN r.confidence < 0.7 THEN 1 END) as low_conf,

                    count(r) as total

                WHERE total > 0

                RETURN high_conf, med_conf, low_conf

            """)



            stats = result.single()

            if stats:

                print(f"\n🎯 Relazioni per Confidence:")

                print(f"  High (≥0.9): {stats['high_conf']:,}")

                print(f"  Medium (0.7-0.9): {stats['med_conf']:,}")

                print(f"  Low (<0.7): {stats['low_conf']:,}")

        except Exception as e:

            pass  # Confidence non ancora disponibile



        print("\n" + "="*70)



def check_layer2_readiness():

    """

    Verifica se il grafo è pronto per Layer 2

    """



    with driver.session() as session:



        print("\n" + "="*70)

        print("🎯 LAYER 2 READINESS CHECK")

        print("="*70)



        # Check TechnicalDocument nodes

        result = session.run("""

            MATCH (d:TechnicalDocument)

            RETURN count(d) as cnt

        """)

        doc_count = result.single()['cnt']



        # Check constraints

        result = session.run("SHOW CONSTRAINTS")

        constraints = [r['name'] for r in result]



        has_required_constraints = (

            any('TechnicalDocument' in c for c in constraints) and

            any('CanonicalLaw' in c for c in constraints)

        )



        # Check indexes

        result = session.run("SHOW INDEXES")

        indexes = [r['name'] for r in result]

        has_indexes = len(indexes) > 0



        ready = doc_count > 0 and has_required_constraints and has_indexes



        if ready:

            print(f"\n✅ READY FOR LAYER 2!")

            print(f"\n📊 Current State:")

            print(f"   - TechnicalDocument nodes: {doc_count:,}")

            print(f"   - Constraints: {len(constraints)} defined")

            print(f"   - Indexes: {len(indexes)} defined")

            print(f"\n🚀 Next Step:")

            print(f"   python3 layer2_extract_relations.py --sync")

        else:

            print(f"\n❌ NOT READY FOR LAYER 2!")

            print(f"\n📊 Issues:")

            if doc_count == 0:

                print(f"   ❌ No TechnicalDocument nodes found")

                print(f"      → Run: python3 02_text_chunks.py")

            else:

                print(f"   ✅ TechnicalDocument nodes: {doc_count:,}")



            if not has_required_constraints:

                print(f"   ❌ Missing required constraints")

                print(f"      → Run: python3 neo4j_schema_setup.py --setup")

            else:

                print(f"   ✅ Constraints: OK")



            if not has_indexes:

                print(f"   ⚠️  No indexes defined (optional but recommended)")

                print(f"      → Run: python3 neo4j_schema_setup.py --setup")



        print("\n" + "="*70)



if __name__ == "__main__":



    import argparse

    parser = argparse.ArgumentParser(description="Neo4j Schema Setup (NO APOC)")

    parser.add_argument('--reset', action='store_true', help='Safe reset (preserva TechnicalDocument)')

    parser.add_argument('--setup', action='store_true', help='Setup schema + constraints')

    parser.add_argument('--cleanup', action='store_true', help='Cleanup isolated nodes')

    parser.add_argument('--merge-dups', action='store_true', help='Mark duplicates')

    parser.add_argument('--stats', action='store_true', help='Print graph statistics')

    parser.add_argument('--check', action='store_true', help='Check Layer 2 readiness')

    parser.add_argument('--all', action='store_true', help='Run setup + stats + check')



    args = parser.parse_args()



    if not any(vars(args).values()) or args.all:

        # Default: setup + stats + check

        setup_schema()

        print_graph_stats()

        check_layer2_readiness()

    else:

        if args.reset:

            confirm = input("⚠️  Reset Layer 2+ data? (yes/no): ")

            if confirm.lower() == 'yes':

                reset_schema_only()

        if args.setup:

            setup_schema()

        if args.cleanup:

            cleanup_isolated_nodes()

        if args.merge_dups:

            merge_duplicates()

        if args.stats:

            print_graph_stats()

        if args.check:

            check_layer2_readiness()



    driver.close()

    print("\n✅ Done!")

