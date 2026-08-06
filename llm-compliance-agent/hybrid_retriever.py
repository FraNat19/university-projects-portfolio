#!/usr/bin/env python3
"""
🎯 HYBRID RETRIEVER - GraphRAG per Workplace Safety

Integrazione Qdrant + Neo4j con tunnel support
"""

from qdrant_client import QdrantClient
from neo4j import GraphDatabase
from sentence_transformers import SentenceTransformer
from typing import List, Dict, Any
import numpy as np

class HybridRetriever:
    """
    Sistema di retrieval ibrido:
    1. Vector search in Qdrant (primary)
    2. Graph enrichment via Neo4j
    3. Filtering (leggi abrogate/modificate)
    4. Reranking finale
    """
    
    def __init__(
        self,
        qdrant_host: str = "127.0.0.1",
        qdrant_port: int = 16333,  # ✅ Porta tunnel corretta
        neo4j_uri: str = "bolt://127.0.0.1:17687",  # ✅ Porta tunnel corretta
        neo4j_user: str = "neo4j",
        neo4j_password: str = "thesis2024",
        model_name: str = "intfloat/multilingual-e5-large"
    ):
        print("🔧 Initializing Hybrid Retriever...")
        
        # Qdrant client
        print(f"   Connecting to Qdrant: {qdrant_host}:{qdrant_port}")
        self.qdrant = QdrantClient(host=qdrant_host, port=qdrant_port)
        
        # Neo4j driver
        print(f"   Connecting to Neo4j: {neo4j_uri}")
        self.neo4j = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_password))
        
        # Test connections
        self._test_connections()
        
        # Embedding model
        print(f"   Loading model: {model_name}")
        self.encoder = SentenceTransformer(model_name)
        
        print("✅ Hybrid Retriever ready!\n")
    
    def _test_connections(self):
        """Test Qdrant and Neo4j connections"""
        try:
            # Test Qdrant
            collections = self.qdrant.get_collections()
            print(f"      ✅ Qdrant OK ({len(collections.collections)} collections)")
            
            # Test Neo4j
            with self.neo4j.session() as session:
                result = session.run("RETURN 1 as test")
                result.single()
                print(f"      ✅ Neo4j OK")
        except Exception as e:
            print(f"      ❌ Connection failed: {str(e)}")
            raise
    
    def retrieve(
        self,
        query: str,
        top_k: int = 10,
        collections: List[str] = None,
        include_graph_context: bool = True,
        filter_abrogated: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Main retrieval method
        
        Args:
            query: User query string
            top_k: Number of results to return
            collections: List of Qdrant collections to search (None = all)
            include_graph_context: Whether to enrich with Neo4j graph
            filter_abrogated: Whether to filter out abrogated laws
        
        Returns:
            List of results with scores and metadata
        """
        print(f"\n{'='*80}")
        print(f"🔍 QUERY: {query}")
        print(f"{'='*80}\n")
        
        # 1. Encode query
        print("1️⃣ Encoding query...")
        query_embedding = self.encoder.encode(query, normalize_embeddings=True)
        
        # 2. Vector search in Qdrant
        print("2️⃣ Searching Qdrant...")
        qdrant_results = self._qdrant_search(query_embedding, top_k, collections)
        
        # 3. Extract law IDs
        print("3️⃣ Extracting law references...")
        law_ids = self._extract_law_ids(qdrant_results)
        print(f"   Found {len(law_ids)} unique law references")
        
        # 4. Enrich with graph context
        graph_context = {}
        if include_graph_context and law_ids:
            print("4️⃣ Enriching with Neo4j graph...")
            graph_context = self._neo4j_enrich(law_ids)
            print(f"   Retrieved context for {len(graph_context)} laws")
        
        # 5. Merge and filter results
        print("5️⃣ Processing results...")
        processed_results = self._process_results(
            qdrant_results,
            graph_context,
            filter_abrogated
        )
        
        # 6. Rerank
        print("6️⃣ Reranking...")
        final_results = self._rerank(processed_results, query)[:top_k]
        
        print(f"\n✅ Returning top {len(final_results)} results\n")
        
        return final_results
    
    def _qdrant_search(
        self,
        embedding: np.ndarray,
        top_k: int,
        collections: List[str] = None
    ) -> Dict[str, List]:
        """Search across Qdrant collections"""
        
        if collections is None:
            collections = [
                'technical_documents_text',
                'legal_articles',
                'legal_context'
            ]
        
        results = {}
        
        for collection in collections:
            try:
                hits = self.qdrant.search(
                    collection_name=collection,
                    query_vector=embedding.tolist(),
                    limit=top_k
                )
                results[collection] = hits
                print(f"   {collection}: {len(hits)} results")
            except Exception as e:
                print(f"   ⚠️ {collection}: {str(e)[:50]}")
                results[collection] = []
        
        return results
    
    def _extract_law_ids(self, qdrant_results: Dict[str, List]) -> List[str]:
        """Extract all law_id references from Qdrant results"""
        
        law_ids = set()
        
        for collection, hits in qdrant_results.items():
            for hit in hits:
                payload = hit.payload
                
                # From technical documents
                if 'laws_cited' in payload:
                    law_ids.update(payload['laws_cited'])
                
                if 'laws_structured' in payload:
                    for law in payload['laws_structured']:
                        if 'law_id' in law:
                            law_ids.add(law['law_id'])
                
                # From legal articles
                if 'parent_law_id' in payload:
                    law_ids.add(payload['parent_law_id'])
                
                # From legal context
                if 'law_id' in payload:
                    law_ids.add(payload['law_id'])
        
        return list(law_ids)
    
    def _neo4j_enrich(self, law_ids: List[str]) -> Dict[str, Dict]:
        """
        Get graph context for laws from Neo4j:
        - Abrogation status
        - Modifications
        - Related laws (RICHIAMA)
        - Semantic context (risks, professions)
        """
        
        with self.neo4j.session() as session:
            result = session.run("""
                UNWIND $law_ids as law_id
                
                // Get law + basic info
                MATCH (law:CanonicalLaw {law_id: law_id})
                
                // Check if abrogated
                OPTIONAL MATCH (abroga_law)-[abroga:ABROGA]->(law)
                
                // Check if modified
                OPTIONAL MATCH (modif_law)-[modif:MODIFICA]->(law)
                
                // Get related laws that this law references
                OPTIONAL MATCH (law)-[richiama:RICHIAMA]->(related)
                WHERE related:CanonicalLaw OR related:Article
                
                // Get semantic context
                OPTIONAL MATCH (law)<-[:GOVERNED_BY]-(risk:RiskCategory)
                OPTIONAL MATCH (law)<-[:REGULATED_BY]-(prof:ProfessionalTag)
                
                // Get technical documents citing this law
                OPTIONAL MATCH (doc:TechnicalDocument)-[:CITES]->(law)
                
                RETURN 
                    law.law_id as law_id,
                    law.titolo as titolo,
                    law.data as data,
                    law.tipo as tipo,
                    law.numero as numero,
                    law.anno as anno,
                    
                    // Abrogation info
                    count(DISTINCT abroga) > 0 as is_abrogated,
                    collect(DISTINCT abroga_law.law_id) as abrogated_by,
                    
                    // Modification info
                    count(DISTINCT modif) > 0 as is_modified,
                    collect(DISTINCT modif_law.law_id) as modified_by,
                    
                    // Related laws
                    collect(DISTINCT CASE 
                        WHEN related:CanonicalLaw THEN related.law_id 
                        WHEN related:Article THEN related.article_id
                    END)[0..5] as related_refs,
                    
                    // Semantic context
                    collect(DISTINCT risk.name)[0..5] as risks,
                    collect(DISTINCT prof.categoria)[0..5] as professions,
                    
                    // Citation count
                    count(DISTINCT doc) as cited_by_docs
            """, law_ids=law_ids)
            
            context = {}
            for row in result:
                law_id = row['law_id']
                context[law_id] = {
                    'law_id': law_id,
                    'titolo': row['titolo'],
                    'data': row['data'],
                    'tipo': row['tipo'],
                    'numero': row['numero'],
                    'anno': row['anno'],
                    'is_abrogated': row['is_abrogated'],
                    'abrogated_by': [x for x in row['abrogated_by'] if x],
                    'is_modified': row['is_modified'],
                    'modified_by': [x for x in row['modified_by'] if x],
                    'related_refs': [x for x in row['related_refs'] if x],
                    'risks': [x for x in row['risks'] if x],
                    'professions': [x for x in row['professions'] if x],
                    'cited_by_docs': row['cited_by_docs']
                }
            
            return context
    
    def _process_results(
        self,
        qdrant_results: Dict[str, List],
        graph_context: Dict[str, Dict],
        filter_abrogated: bool
    ) -> List[Dict]:
        """Process and filter results"""
        
        all_results = []
        
        for collection, hits in qdrant_results.items():
            for hit in hits:
                payload = hit.payload
                
                # Extract relevant law_ids for this result
                result_law_ids = set()
                
                if 'laws_cited' in payload:
                    result_law_ids.update(payload['laws_cited'])
                if 'parent_law_id' in payload:
                    result_law_ids.add(payload['parent_law_id'])
                if 'law_id' in payload:
                    result_law_ids.add(payload['law_id'])
                
                # Build result object
                result = {
                    'score': hit.score,
                    'collection': collection,
                    'payload': payload,
                    'law_ids': list(result_law_ids),
                    'warnings': [],
                    'graph_context': {}
                }
                
                # Add graph context
                for law_id in result_law_ids:
                    if law_id in graph_context:
                        context = graph_context[law_id]
                        result['graph_context'][law_id] = context
                        
                        # Add warnings
                        if context['is_abrogated']:
                            result['warnings'].append(
                                f"⚠️ {law_id} abrogata da: {', '.join(context['abrogated_by'])}"
                            )
                        if context['is_modified']:
                            result['warnings'].append(
                                f"ℹ️ {law_id} modificata da: {', '.join(context['modified_by'])}"
                            )
                
                # Filter if abrogated and filter_abrogated=True
                if filter_abrogated:
                    is_all_abrogated = all(
                        graph_context.get(lid, {}).get('is_abrogated', False)
                        for lid in result_law_ids
                        if lid in graph_context
                    )
                    if is_all_abrogated and result_law_ids:
                        continue  # Skip fully abrogated results
                
                all_results.append(result)
        
        return all_results
    
    def _rerank(self, results: List[Dict], query: str) -> List[Dict]:
        """Rerank results based on multiple factors"""
        
        for result in results:
            score = result['score']
            
            # Boost if no warnings
            if not result['warnings']:
                score *= 1.1
            
            # Boost if has graph context
            if result['graph_context']:
                score *= 1.05
            
            # Boost if cited by many docs
            for context in result['graph_context'].values():
                if context.get('cited_by_docs', 0) > 5:
                    score *= 1.1
            
            # Boost technical docs
            if result['collection'] == 'technical_documents_text':
                score *= 1.2
            
            # Penalize if has warnings
            if result['warnings']:
                score *= 0.95
            
            result['final_score'] = score
        
        # Sort by final score
        results.sort(key=lambda x: x['final_score'], reverse=True)
        
        return results
    
    def print_results(self, results: List[Dict], max_results: int = 5):
        """Pretty print results"""
        
        print(f"\n{'='*80}")
        print(f"📊 TOP {min(len(results), max_results)} RESULTS")
        print(f"{'='*80}\n")
        
        for i, result in enumerate(results[:max_results], 1):
            print(f"{i}. Score: {result['final_score']:.4f} (original: {result['score']:.4f})")
            print(f"   Collection: {result['collection']}")
            
            # Title/description
            payload = result['payload']
            if 'titolo' in payload:
                print(f"   📄 {payload['titolo'][:100]}")
            elif 'title' in payload:
                print(f"   📄 {payload['title'][:100]}")
            elif 'doc_id' in payload:
                print(f"   📄 Doc: {payload['doc_id']}")
            
            # Warnings
            if result['warnings']:
                for warning in result['warnings']:
                    print(f"   {warning}")
            
            # Graph context
            if result['graph_context']:
                print(f"   🔗 Referenced laws ({len(result['graph_context'])}):")
                for law_id, context in list(result['graph_context'].items())[:3]:
                    status = ""
                    if context['is_abrogated']:
                        status = " [ABROGATA]"
                    elif context['is_modified']:
                        status = " [MODIFICATA]"
                    
                    print(f"      • {law_id}{status}")
                    if context.get('titolo'):
                        print(f"        {context['titolo'][:80]}")
                    if context.get('risks'):
                        print(f"        Rischi: {', '.join(context['risks'][:3])}")
            
            print()
    
    def close(self):
        """Close connections"""
        self.neo4j.close()
        print("✅ Connections closed")


if __name__ == "__main__":
    # Test retriever
    retriever = HybridRetriever()
    
    # Test query
    results = retriever.retrieve(
        "rischio chimico amianto esposizione lavoratori",
        top_k=5
    )
    
    retriever.print_results(results)
    retriever.close()
