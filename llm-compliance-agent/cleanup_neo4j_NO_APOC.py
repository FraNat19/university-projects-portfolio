#!/usr/bin/env python3

"""

🗑️ CLEANUP NEO4J - UNIVERSAL VERSION (NO APOC, NO BACKUP)



COMPATIBLE:

- Works without APOC plugin

- No filesystem backup needed

- Pure Cypher queries

"""



from neo4j import GraphDatabase

from datetime import datetime



NEO4J_URI = "bolt://127.0.0.1:7687"

NEO4J_AUTH = ("neo4j", "thesis2024")



driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)





def print_stats():

    """Print current stats"""

    

    with driver.session() as session:

        print("\n📊 CURRENT STATS")

        print("="*70 + "\n")

        

        result = session.run("""

            MATCH (n)

            WITH labels(n)[0] as label, count(n) as cnt

            WHERE label IS NOT NULL

            RETURN label, cnt

            ORDER BY cnt DESC

        """)

        

        for r in result:

            print(f"  {r['label']}: {r['cnt']:,}")

        

        # Duplicates

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WITH law.numero as num, law.anno as anno, count(*) as cnt

            WHERE cnt > 1

            RETURN count(*) as dup_groups, sum(cnt) as total_dups

        """)

        

        stats = result.single()

        if stats and stats['dup_groups'] > 0:

            print(f"\n⚠️  Duplicates: {stats['dup_groups']} groups ({stats['total_dups']} nodes)")

        else:

            print(f"\n✅ No duplicates")

        

        print("\n" + "="*70 + "\n")





def simple_merge_duplicates():

    """

    ✅ SIMPLE MERGE (NO APOC NEEDED)

    

    Strategy:

    1. For each (numero, anno) group

    2. Keep FIRST node with enriched=true OR highest score

    3. DELETE other nodes (relationships CASCADE)

    """

    

    print("🔍 Searching duplicates...\n")

    

    with driver.session() as session:

        

        # Find duplicate groups

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WHERE law.numero IS NOT NULL AND law.anno IS NOT NULL

            WITH law.numero as num, law.anno as anno, collect(elementId(law)) as ids

            WHERE size(ids) > 1

            RETURN num, anno, ids, size(ids) as count

        """)

        

        dup_groups = list(result)

        

        if not dup_groups:

            print("✅ No duplicates!\n")

            return 0

        

        print(f"⚠️  Found {len(dup_groups)} duplicate groups\n")

        

        total_deleted = 0

        

        for group in dup_groups:

            num = group['num']

            anno = group['anno']

            ids = group['ids']

            

            print(f"📋 Processing {num}/{anno} ({len(ids)} duplicates)")

            

            # Find keeper (best node)

            result = session.run("""

                UNWIND $ids as node_id

                MATCH (law:CanonicalLaw)

                WHERE elementId(law) = node_id

                

                WITH law, elementId(law) as eid,

                     CASE WHEN law.enriched = true THEN 100 ELSE 0 END +

                     CASE WHEN law.urn IS NOT NULL THEN 50 ELSE 0 END +

                     CASE WHEN law.titolo IS NOT NULL THEN 30 ELSE 0 END +

                     CASE WHEN law.has_qdrant_legal = true THEN 20 ELSE 0 END as score

                

                ORDER BY score DESC, eid

                LIMIT 1

                

                RETURN eid as keeper_id, law.law_id as law_id, 

                       law.tipo as tipo, score

            """, ids=ids)

            

            keeper = result.single()

            keeper_id = keeper['keeper_id']

            keeper_law_id = keeper['law_id']

            tipo = keeper['tipo']

            

            # Build canonical law_id if missing

            if not keeper_law_id:

                tipo_map = {

                    'decreto.legislativo': 'dlgs',

                    'decreto.del.presidente.della.repubblica': 'dpr',

                    'decreto.legge': 'dl',

                    'legge': 'legge',

                    'codice': 'codice',

                    'testo.unico': 'tu',

                    'regio.decreto': 'rd'

                }

                tipo_short = tipo_map.get(tipo, 'unknown')

                keeper_law_id = f"{tipo_short}-{num}-{anno}"

            

            print(f"   ✅ Keeper: {keeper_law_id}")

            

            # Update keeper

            session.run("""

                MATCH (law:CanonicalLaw)

                WHERE elementId(law) = $keeper_id

                SET law.law_id = $law_id,

                    law.qdrant_prefix = $law_id,

                    law.cleaned = true,

                    law.cleaned_at = datetime()

            """, keeper_id=keeper_id, law_id=keeper_law_id)

            

            # ✅ SIMPLE DELETE (relationships cascade automatically)

            result = session.run("""

                UNWIND $ids as node_id

                MATCH (law:CanonicalLaw)

                WHERE elementId(law) = node_id

                  AND elementId(law) <> $keeper_id

                

                DETACH DELETE law

                RETURN count(law) as deleted

            """, ids=ids, keeper_id=keeper_id)

            

            deleted = result.single()['deleted']

            print(f"   🗑️  Deleted {deleted} duplicates")

            

            total_deleted += deleted

        

        print(f"\n✅ Total removed: {total_deleted}\n")

        return total_deleted





def cleanup_missing_law_id():

    """Fix nodes with missing law_id"""

    

    print("🔧 Fixing missing law_id...\n")

    

    with driver.session() as session:

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WHERE law.law_id IS NULL

              AND law.numero IS NOT NULL

              AND law.anno IS NOT NULL

            

            WITH law, law.numero as num, law.anno as anno, law.tipo as tipo

            

            SET law.law_id = CASE

                WHEN tipo = 'decreto.legislativo' THEN 'dlgs-' + num + '-' + anno

                WHEN tipo = 'decreto.del.presidente.della.repubblica' THEN 'dpr-' + num + '-' + anno

                WHEN tipo = 'decreto.legge' THEN 'dl-' + num + '-' + anno

                WHEN tipo = 'legge' THEN 'legge-' + num + '-' + anno

                WHEN tipo = 'codice' THEN 'codice-' + num + '-' + anno

                ELSE 'unknown-' + num + '-' + anno

            END,

            law.qdrant_prefix = CASE

                WHEN tipo = 'decreto.legislativo' THEN 'dlgs-' + num + '-' + anno

                WHEN tipo = 'decreto.del.presidente.della.repubblica' THEN 'dpr-' + num + '-' + anno

                WHEN tipo = 'decreto.legge' THEN 'dl-' + num + '-' + anno

                WHEN tipo = 'legge' THEN 'legge-' + num + '-' + anno

                WHEN tipo = 'codice' THEN 'codice-' + num + '-' + anno

                ELSE 'unknown-' + num + '-' + anno

            END

            

            RETURN count(law) as fixed

        """)

        

        fixed = result.single()['fixed']

        print(f"✅ Fixed {fixed:,} nodes with missing law_id\n")

        return fixed





def reset_constraints():

    """Reset constraints (law_id as primary)"""

    

    print("🔧 Resetting constraints...\n")

    

    with driver.session() as session:

        try:

            session.run("DROP CONSTRAINT CanonicalLaw_urn IF EXISTS")

        except:

            pass

        

        try:

            session.run("""

                CREATE CONSTRAINT CanonicalLaw_law_id IF NOT EXISTS

                FOR (law:CanonicalLaw) REQUIRE law.law_id IS UNIQUE

            """)

            print("✅ Constraint created: CanonicalLaw_law_id UNIQUE\n")

        except Exception as e:

            if "already exists" in str(e).lower():

                print("✅ Constraint already exists\n")





def verify():

    """Final verification"""

    

    print("\n" + "="*70)

    print("🔍 VERIFICATION")

    print("="*70 + "\n")

    

    with driver.session() as session:

        

        # Duplicates check

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WITH law.numero as num, law.anno as anno, count(*) as cnt

            WHERE cnt > 1

            RETURN count(*) as dup_groups

        """)

        dup = result.single()['dup_groups']

        

        print(f"Duplicates: {dup} groups" + (" ⚠️" if dup > 0 else " ✅"))

        

        # Stats

        result = session.run("""

            MATCH (law:CanonicalLaw)

            RETURN 

                count(law) as total,

                count(DISTINCT law.law_id) as unique,

                count(CASE WHEN law.law_id IS NULL THEN 1 END) as missing

        """)

        stats = result.single()

        

        print(f"Total CanonicalLaw: {stats['total']:,}")

        print(f"Unique law_id: {stats['unique']:,}")

        print(f"Missing law_id: {stats['missing']:,}")

        

        # Format check

        result = session.run("""

            MATCH (law:CanonicalLaw)

            WHERE law.law_id =~ '[a-z]+-\\d+-\\d{4}'

            RETURN count(law) as correct

        """)

        correct = result.single()['correct']

        print(f"Correct format: {correct:,}")

        

        # Sample

        result = session.run("""

            MATCH (law:CanonicalLaw)

            RETURN law.law_id

            LIMIT 5

        """)

        

        print(f"\nSample:")

        for r in result:

            print(f"  {r['law.law_id']}")

        

        print("\n" + "="*70 + "\n")





if __name__ == "__main__":

    

    print("\n" + "="*70)

    print("🗑️  NEO4J CLEANUP (UNIVERSAL - NO APOC)")

    print("="*70 + "\n")

    

    # Before stats

    print_stats()

    

    # Confirm

    confirm = input("Continue with cleanup? (yes/no): ")

    if confirm.lower() != 'yes':

        print("\n❌ Aborted\n")

        driver.close()

        exit(0)

    

    # Execute

    print("")

    deleted = simple_merge_duplicates()

    fixed = cleanup_missing_law_id()

    reset_constraints()

    

    # After stats

    verify()

    

    print("="*70)

    print("✅ CLEANUP COMPLETE!")

    print("="*70)

    print(f"   Duplicates removed: {deleted}")

    print(f"   Missing law_id fixed: {fixed}")

    print("="*70 + "\n")

    

    driver.close()
