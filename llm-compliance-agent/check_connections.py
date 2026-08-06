#!/usr/bin/env python3
"""
🔍 CHECK CONNECTIONS - Verifica Qdrant + Neo4j

Verifica che tunnel e servizi siano accessibili prima di testare
"""

from qdrant_client import QdrantClient
from neo4j import GraphDatabase
import sys

def check_qdrant():
    """Verifica connessione Qdrant"""
    print("\n📦 Checking Qdrant...")
    print("   Host: 127.0.0.1:16333")
    
    try:
        client = QdrantClient(host="127.0.0.1", port=16333)
        collections = client.get_collections()
        
        print(f"   ✅ Connected!")
        print(f"   Collections found: {len(collections.collections)}")
        
        for coll in collections.collections:
            info = client.get_collection(coll.name)
            print(f"      • {coll.name}: {info.points_count:,} points")
        
        return True
        
    except Exception as e:
        print(f"   ❌ Failed: {str(e)}")
        return False

def check_neo4j():
    """Verifica connessione Neo4j"""
    print("\n🔗 Checking Neo4j...")
    print("   URI: bolt://127.0.0.1:17687")
    
    try:
        driver = GraphDatabase.driver(
            "bolt://127.0.0.1:17687",
            auth=("neo4j", "thesis2024")
        )
        
        with driver.session() as session:
            # Count nodes
            result = session.run("""
                MATCH (n)
                RETURN labels(n)[0] as label, count(n) as count
                ORDER BY count DESC
            """)
            
            print(f"   ✅ Connected!")
            print(f"   Node counts:")
            
            for record in result:
                print(f"      • {record['label']}: {record['count']:,}")
            
            # Count relations
            result = session.run("""
                MATCH ()-[r]->()
                RETURN type(r) as rel_type, count(r) as count
                ORDER BY count DESC
                LIMIT 10
            """)
            
            print(f"   Relationship counts (top 10):")
            
            for record in result:
                print(f"      • {record['rel_type']}: {record['count']:,}")
        
        driver.close()
        return True
        
    except Exception as e:
        print(f"   ❌ Failed: {str(e)}")
        return False

def check_tunnel():
    """Verifica che tunnel sia attivo"""
    import subprocess
    
    print("\n🔌 Checking SSH Tunnel...")
    
    try:
        result = subprocess.run(
            ['lsof', '-iTCP:16333', '-sTCP:LISTEN', '-nP'],
            capture_output=True,
            text=True
        )
        
        if result.stdout:
            print("   ✅ Tunnel appears to be running")
            print("   Port 16333 (Qdrant): LISTENING")
        else:
            print("   ⚠️ Port 16333 not listening")
            print("   Run: ./tunnel_stack.sh start <JOBID>")
            return False
        
        return True
        
    except Exception as e:
        print(f"   ⚠️ Could not check: {str(e)}")
        return True  # Don't fail on this

def main():
    print("\n" + "="*80)
    print("🔍 CONNECTION CHECK - Qdrant + Neo4j")
    print("="*80)
    
    # Check tunnel
    tunnel_ok = check_tunnel()
    
    # Check services
    qdrant_ok = check_qdrant()
    neo4j_ok = check_neo4j()
    
    # Summary
    print("\n" + "="*80)
    print("📊 SUMMARY")
    print("="*80)
    
    print(f"   Tunnel: {'✅ OK' if tunnel_ok else '❌ NOT OK'}")
    print(f"   Qdrant: {'✅ OK' if qdrant_ok else '❌ NOT OK'}")
    print(f"   Neo4j:  {'✅ OK' if neo4j_ok else '❌ NOT OK'}")
    
    if qdrant_ok and neo4j_ok:
        print("\n✅ All systems GO! Ready to test retriever.")
        print("\nRun: python3 test_retriever.py")
        sys.exit(0)
    else:
        print("\n❌ Some systems are not ready.")
        print("\n🔧 Troubleshooting:")
        
        if not tunnel_ok:
            print("   1. Start tunnel:")
            print("      sbatch slurm_stack_qdrant_neo4j.sbatch")
            print("      ./tunnel_stack.sh start <JOBID>")
        
        if not qdrant_ok:
            print("   2. Check Qdrant is running on compute node")
        
        if not neo4j_ok:
            print("   3. Check Neo4j is running on compute node")
        
        sys.exit(1)

if __name__ == "__main__":
    main()
