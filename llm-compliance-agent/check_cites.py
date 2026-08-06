#!/usr/bin/env python3
"""
Verifica e cancella CITES layer 1
"""
import argparse
from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = "thesis2024"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true', help='Solo conta, non cancella')
    parser.add_argument('--delete', action='store_true', help='Cancella CITES layer 1')
    args = parser.parse_args()
    
    if not args.check and not args.delete:
        print("Uso: python check_cites.py --check    (solo conta)")
        print("     python check_cites.py --delete   (cancella)")
        return
    
    d = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    
    with d.session() as s:
        # Conta CITES layer 1
        r = s.run("MATCH ()-[c:CITES {layer: 1}]->() RETURN count(c) as total")
        cites_layer1 = r.single()['total']
        print(f"CITES con layer=1: {cites_layer1}")
        
        # Conta altre CITES
        r = s.run("MATCH ()-[c:CITES]->() WHERE c.layer IS NULL OR c.layer <> 1 RETURN count(c) as total")
        cites_other = r.single()['total']
        print(f"CITES altre (layer!=1): {cites_other}")
        
        # Conta HAS_ARTICLE
        r = s.run("MATCH ()-[h:HAS_ARTICLE]->() RETURN count(h) as total")
        has_art = r.single()['total']
        print(f"HAS_ARTICLE: {has_art} (NON verranno toccate)")
        
        if args.delete:
            if cites_layer1 > 0:
                print(f"\nCancello {cites_layer1} CITES layer=1...")
                r = s.run("MATCH ()-[c:CITES {layer: 1}]->() DELETE c RETURN count(c) as deleted")
                deleted = r.single()['deleted']
                print(f"✅ Cancellate: {deleted} relazioni CITES")
            else:
                print("\n✅ Nessuna CITES layer=1 da cancellare")
    
    d.close()
    print("\nDone!")

if __name__ == "__main__":
    main()
