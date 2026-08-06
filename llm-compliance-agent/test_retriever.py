#!/usr/bin/env python3
"""
🧪 TEST SUITE - Hybrid Retriever GraphRAG

Test queries per workplace safety system
"""

from hybrid_retriever import HybridRetriever
import sys

def test_query(retriever, query_name, query_text, top_k=5):
    """Run a single test query"""
    
    print("\n" + "="*80)
    print(f"🧪 TEST: {query_name}")
    print("="*80)
    print(f"Query: {query_text}\n")
    
    try:
        results = retriever.retrieve(query_text, top_k=top_k)
        retriever.print_results(results, max_results=top_k)
        return True
    except Exception as e:
        print(f"❌ ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    print("\n" + "="*80)
    print("🚀 HYBRID RETRIEVER - TEST SUITE")
    print("="*80)
    
    # Initialize retriever
    try:
        retriever = HybridRetriever()
    except Exception as e:
        print(f"\n❌ Failed to initialize retriever: {str(e)}")
        print("\n⚠️ Make sure:")
        print("   1. Tunnel is running (./tunnel_stack.sh status)")
        print("   2. Qdrant and Neo4j are accessible")
        sys.exit(1)
    
    # Test queries
    tests = [
        {
            'name': 'TEST 1: Rischio Chimico - Amianto',
            'query': 'rischio chimico amianto esposizione lavoratori protezione',
            'top_k': 5
        },
        {
            'name': 'TEST 2: DPI Cantiere Edile',
            'query': 'dispositivi protezione individuale DPI obbligatori cantiere edile',
            'top_k': 5
        },
        {
            'name': 'TEST 3: D.Lgs 81/2008 Testo Unico',
            'query': 'decreto legislativo 81 2008 testo unico sicurezza lavoro',
            'top_k': 5
        },
        {
            'name': 'TEST 4: Rischio Caduta dall\'Alto',
            'query': 'rischio caduta dall\'alto lavori in quota ponteggi',
            'top_k': 5
        },
        {
            'name': 'TEST 5: Sorveglianza Sanitaria',
            'query': 'sorveglianza sanitaria lavoratori medico competente visite',
            'top_k': 5
        }
    ]
    
    # Run tests
    results = []
    for test in tests:
        success = test_query(
            retriever,
            test['name'],
            test['query'],
            test['top_k']
        )
        results.append((test['name'], success))
        
        input("\n⏸️  Press ENTER to continue to next test...")
    
    # Summary
    print("\n" + "="*80)
    print("📊 TEST SUMMARY")
    print("="*80)
    
    for test_name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{status} - {test_name}")
    
    passed = sum(1 for _, s in results if s)
    total = len(results)
    print(f"\nTotal: {passed}/{total} tests passed")
    
    # Cleanup
    retriever.close()
    
    print("\n✅ Test suite complete!\n")

if __name__ == "__main__":
    main()
