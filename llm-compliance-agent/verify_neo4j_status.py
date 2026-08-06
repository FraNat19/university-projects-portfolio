import os
#!/usr/bin/env python3
"""
Verifica stato completo Neo4j dopo i layer
"""
from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = os.environ.get("NEO4J_PASSWORD", "")

def main():
    d = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    
    with d.session() as s:
        print("\n" + "="*70)
        print("📊 STATO COMPLETO NEO4J")
        print("="*70)
        
        # ====== NODI ======
        print("\n🔷 NODI:")
        
        r = s.run("MATCH (n:CanonicalLaw) RETURN count(n) as c")
        print(f"   CanonicalLaw: {r.single()['c']:,}")
        
        r = s.run("MATCH (n:Article) RETURN count(n) as c")
        print(f"   Article: {r.single()['c']:,}")
        
        r = s.run("MATCH (n:TechnicalDocument) RETURN count(n) as c")
        print(f"   TechnicalDocument: {r.single()['c']:,}")
        
        # ====== RELAZIONI ======
        print("\n🔗 RELAZIONI:")
        
        r = s.run("MATCH ()-[r:HAS_ARTICLE]->() RETURN count(r) as c")
        print(f"   HAS_ARTICLE: {r.single()['c']:,}")
        
        r = s.run("MATCH ()-[r:CITES]->() RETURN count(r) as c")
        print(f"   CITES (totale): {r.single()['c']:,}")
        
        r = s.run("MATCH ()-[r:CITES {layer: 1}]->() RETURN count(r) as c")
        print(f"   CITES layer=1: {r.single()['c']:,}")
        
        r = s.run("MATCH ()-[r:CITES]->() WHERE r.granularity = 'article' RETURN count(r) as c")
        print(f"   CITES granularity=article: {r.single()['c']:,}")
        
        r = s.run("MATCH ()-[r:CITES]->() WHERE r.granularity = 'law' RETURN count(r) as c")
        print(f"   CITES granularity=law: {r.single()['c']:,}")
        
        # ====== PROVENIENZA CanonicalLaw ======
        print("\n📚 CanonicalLaw PER PROVENIENZA:")
        
        r = s.run("""
            MATCH (cl:CanonicalLaw)
            RETURN 
                count(CASE WHEN cl.has_qdrant_legal = true THEN 1 END) as from_normattiva,
                count(CASE WHEN cl.has_qdrant_technical = true THEN 1 END) as from_technical,
                count(CASE WHEN cl.has_qdrant_legal = true AND cl.has_qdrant_technical = true THEN 1 END) as both,
                count(CASE WHEN cl.enriched = true THEN 1 END) as enriched
        """)
        rec = r.single()
        print(f"   Da Normattiva (legal_context): {rec['from_normattiva']:,}")
        print(f"   Da Technical docs: {rec['from_technical']:,}")
        print(f"   Da entrambi: {rec['both']:,}")
        print(f"   Enriched (layer1b): {rec['enriched']:,}")
        
        # ====== PROVENIENZA Article ======
        print("\n📄 Article PER PROVENIENZA:")
        
        r = s.run("""
            MATCH (a:Article)
            RETURN 
                count(CASE WHEN a.from_qdrant = true THEN 1 END) as from_qdrant,
                count(CASE WHEN a.from_xml = true THEN 1 END) as from_xml,
                count(CASE WHEN a.created_from = 'layer1_citation' THEN 1 END) as from_layer1
        """)
        rec = r.single()
        print(f"   Da Qdrant (legal_articles): {rec['from_qdrant']:,}")
        print(f"   Da XML (layer1b): {rec['from_xml']:,}")
        print(f"   Creati da Layer1 (CITES): {rec['from_layer1']:,}")
        
        # ====== TOP LEGGI CITATE ======
        print("\n🏆 TOP 10 LEGGI PIÙ CITATE:")
        
        r = s.run("""
            MATCH (doc:TechnicalDocument)-[c:CITES]->(target)
            WITH CASE 
                WHEN 'Article' IN labels(target) THEN target.parent_law_id
                ELSE target.law_id
            END as law_id, count(c) as citations
            RETURN law_id, citations
            ORDER BY citations DESC
            LIMIT 10
        """)
        for rec in r:
            print(f"   {rec['law_id']}: {rec['citations']} citazioni")
        
        # ====== SAMPLE CITES CON GRANULARITÀ ======
        print("\n📋 SAMPLE CITES (Article-level):")
        
        r = s.run("""
            MATCH (doc:TechnicalDocument)-[c:CITES {granularity: 'article'}]->(art:Article)
            RETURN doc.doc_id as doc, art.parent_law_id as law, art.article_num as art_num, c.full_citation as full_cit
            LIMIT 5
        """)
        for rec in r:
            print(f"   {rec['doc'][:35]}... → {rec['law']} art.{rec['art_num']}")
            print(f"      full_citation: {rec['full_cit']}")
        
        print("\n" + "="*70)
    
    d.close()
    print("Done!")

if __name__ == "__main__":
    main()
