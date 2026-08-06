import os
#!/usr/bin/env python3
"""
🎯 HYBRID RETRIEVER v3 - GraphRAG per Workplace Safety
IMPROVED VERSION with:
- Better reranking (keyword matching)
- Category filtering (remove garbage)
- Balanced boosts
"""

import requests
import numpy as np
import re
from neo4j import GraphDatabase
from sentence_transformers import SentenceTransformer
from typing import List, Dict, Any, Optional, Set

# ============================================
# DOMAIN CONFIGURATION
# ============================================

# Categorie valide (dal tuo dominio)
VALID_CATEGORIES = {
    "Rischio chimico e cancerogeno", "Rischio biologico", 
    "Rischio fisico (rumore, vibrazioni, radiazioni)",
    "Prevenzione incendi ed esplosioni", "Sicurezza macchine e attrezzature", 
    "Edilizia e cantieri", "Agricoltura e silvicoltura", "Trasporti e logistica", 
    "Settore sanitario", "DPI e dispositivi di protezione",
    "Formazione e addestramento", "Sorveglianza sanitaria", 
    "Normativa e adempimenti", "Valutazione e gestione del rischio", 
    "Infortuni e malattie professionali"
}

# Pattern per identificare categorie garbage
GARBAGE_PATTERNS = [
    r'^[A-Za-z]{5,}\s+e\s+[a-z]{5,}$',  # "Aafkag e arialmt"
    r'(.+)\1',  # Duplicati come "campionamento campionamento"
    r'^.{1,3}$',  # Troppo corte
    r'^altro$',  # Troppo generico
    r'[^a-zA-ZàèéìòùÀÈÉÌÒÙ\s\-]',  # Caratteri strani
]

# Keywords per boosting per categoria
CATEGORY_KEYWORDS = {
    "amianto": ["amianto", "asbesto", "fibre", "mesotelioma", "bonifica amianto"],
    "chimico": ["chimico", "cancerogeno", "mutageno", "sostanze", "agenti chimici", "esposizione"],
    "cantiere": ["cantiere", "edilizia", "costruzioni", "ponteggi", "scavi", "caduta"],
    "dpi": ["dpi", "dispositivi", "protezione individuale", "respiratori", "maschere", "guanti"],
    "rumore": ["rumore", "vibrazioni", "radiazioni", "acustico", "ipoacusia"],
    "incendio": ["incendio", "antincendio", "emergenza", "evacuazione", "estintori"],
    "biologico": ["biologico", "virus", "batteri", "infezione", "contaminazione"],
    "ergonomico": ["ergonomico", "postura", "movimentazione", "carichi", "schermo"],
    "macchine": ["macchine", "attrezzature", "impianti", "manutenzione"],
    "sanitario": ["sanitario", "ospedale", "infermieri", "medici", "pazienti"],
}


def is_garbage_category(category: str) -> bool:
    """Check if a category is garbage/noise"""
    if not category or not isinstance(category, str):
        return True
    
    category_lower = category.lower().strip()
    
    # Check against valid categories
    if category in VALID_CATEGORIES:
        return False
    
    # Check garbage patterns
    for pattern in GARBAGE_PATTERNS:
        if re.match(pattern, category_lower):
            return True
    
    # Check for nonsense (no vowels, too many consonants in a row)
    if not re.search(r'[aeiouàèéìòù]', category_lower):
        return True
    
    if re.search(r'[bcdfghjklmnpqrstvwxz]{5,}', category_lower):
        return True
    
    return False


def extract_query_keywords(query: str) -> Set[str]:
    """Extract meaningful keywords from query"""
    # Normalize
    query_lower = query.lower()
    
    # Remove stop words
    stop_words = {'il', 'lo', 'la', 'i', 'gli', 'le', 'un', 'una', 'uno', 
                  'di', 'a', 'da', 'in', 'con', 'su', 'per', 'tra', 'fra',
                  'e', 'o', 'ma', 'che', 'come', 'quando', 'dove', 'quale',
                  'del', 'della', 'dei', 'delle', 'al', 'alla', 'ai', 'alle'}
    
    words = set(re.findall(r'\b\w+\b', query_lower))
    keywords = words - stop_words
    
    # Add compound terms
    for category, terms in CATEGORY_KEYWORDS.items():
        for term in terms:
            if term in query_lower:
                keywords.add(category)
    
    return keywords


def calculate_keyword_relevance(text: str, query_keywords: Set[str]) -> float:
    """Calculate how relevant a text is to query keywords"""
    if not text or not query_keywords:
        return 0.0
    
    text_lower = text.lower()
    
    matches = 0
    total_weight = 0
    
    for keyword in query_keywords:
        # Exact match
        if keyword in text_lower:
            matches += 1
            # Weight by keyword length (longer = more specific)
            weight = min(len(keyword) / 5, 2.0)
            total_weight += weight
        
        # Partial match (for compound words)
        elif any(keyword in word for word in text_lower.split()):
            matches += 0.5
            total_weight += 0.5
    
    if len(query_keywords) == 0:
        return 0.0
    
    # Normalize by number of keywords
    relevance = total_weight / len(query_keywords)
    
    return min(relevance, 2.0)  # Cap at 2.0


class HybridRetrieverV3:
    """
    Sistema di retrieval ibrido v3:
    1. Vector search in Qdrant
    2. Graph enrichment via Neo4j
    3. Keyword-based reranking
    4. Category filtering
    """
    
    def __init__(
        self,
        qdrant_host: str = "127.0.0.1",
        qdrant_port: int = None,  # Auto-detect: 6333 local, 16333 tunnel
        neo4j_uri: str = None,    # Auto-detect: 7687 local, 17687 tunnel
        neo4j_user: str = "neo4j",
        neo4j_password: str = os.getenv("NEO4J_PASSWORD", ""),
        model_name: str = "intfloat/multilingual-e5-large"
    ):
        import os
        
        # Auto-detect ports from environment or use defaults
        if qdrant_port is None:
            qdrant_port = int(os.getenv("QDRANT_PORT", "16333"))
        
        if neo4j_uri is None:
            neo4j_port = os.getenv("NEO4J_PORT", "17687")
            neo4j_uri = os.getenv("NEO4J_URI", f"bolt://127.0.0.1:{neo4j_port}")
        print("🔧 Initializing Hybrid Retriever v3...")
        
        self.qdrant_url = f"http://{qdrant_host}:{qdrant_port}"
        print(f"   Connecting to Qdrant: {self.qdrant_url}")
        
        print(f"   Connecting to Neo4j: {neo4j_uri}")
        self.neo4j = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_password))
        
        self._test_connections()
        
        print(f"   Loading model: {model_name}")
        self.encoder = SentenceTransformer(model_name)
        
        print("✅ Hybrid Retriever v3 ready!\n")
    
    def _test_connections(self):
        """Test connections"""
        try:
            response = requests.get(f"{self.qdrant_url}/collections", timeout=5)
            if response.status_code == 200:
                data = response.json()
                collections = data.get('result', {}).get('collections', [])
                print(f"      ✅ Qdrant OK ({len(collections)} collections)")
            
            with self.neo4j.session() as session:
                result = session.run("""
                    MATCH (n) 
                    WITH labels(n)[0] as label, count(n) as cnt 
                    RETURN label, cnt ORDER BY cnt DESC LIMIT 5
                """)
                print(f"      ✅ Neo4j OK")
                for row in result:
                    print(f"         - {row['label']}: {row['cnt']:,}")
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
        """Main retrieval with improved reranking"""
        
        print(f"\n{'='*80}")
        print(f"🔍 QUERY: {query}")
        print(f"{'='*80}\n")
        
        # Extract query keywords for reranking
        query_keywords = extract_query_keywords(query)
        print(f"📌 Query keywords: {query_keywords}\n")
        
        # 1. Encode
        print("1️⃣ Encoding query...")
        query_embedding = self.encoder.encode(query, normalize_embeddings=True)
        
        # 2. Vector search
        print("2️⃣ Searching Qdrant...")
        qdrant_results = self._qdrant_search(query_embedding, top_k * 3, collections)
        
        # 3. Extract law IDs
        print("3️⃣ Extracting law references...")
        law_ids = self._extract_law_ids(qdrant_results)
        print(f"   Found {len(law_ids)} unique law references")
        
        # 4. Graph enrichment
        graph_context = {}
        if include_graph_context and law_ids:
            print("4️⃣ Enriching with Neo4j graph...")
            graph_context = self._neo4j_enrich(law_ids)
            print(f"   Retrieved context for {len(graph_context)} laws")
        
        # 5. Process
        print("5️⃣ Processing results...")
        processed_results = self._process_results(
            qdrant_results, graph_context, filter_abrogated
        )
        
        # 6. IMPROVED Reranking
        print("6️⃣ Reranking with keyword matching...")
        final_results = self._rerank_v3(processed_results, query, query_keywords)[:top_k]
        
        print(f"\n✅ Returning top {len(final_results)} results\n")
        
        return final_results
    
    def _qdrant_search(self, embedding: np.ndarray, top_k: int, collections: List[str] = None) -> Dict[str, List]:
        """Search Qdrant via REST"""
        
        if collections is None:
            collections = ['technical_documents_text', 'legal_articles', 'legal_context']
        
        results = {}
        
        for collection in collections:
            try:
                url = f"{self.qdrant_url}/collections/{collection}/points/search"
                payload = {
                    "vector": embedding.tolist(),
                    "limit": top_k,
                    "with_payload": True,
                    "with_vector": False
                }
                
                response = requests.post(url, json=payload, timeout=30)
                
                if response.status_code == 200:
                    data = response.json()
                    hits = data.get('result', [])
                    results[collection] = hits
                    print(f"   {collection}: {len(hits)} results")
                else:
                    results[collection] = []
                    
            except Exception as e:
                print(f"   ⚠️ {collection}: {str(e)[:50]}")
                results[collection] = []
        
        return results
    
    def _extract_law_ids(self, qdrant_results: Dict[str, List]) -> List[str]:
        """Extract law IDs from results"""
        
        law_ids = set()
        
        for collection, hits in qdrant_results.items():
            for hit in hits:
                payload = hit.get('payload', {})
                if not payload:
                    continue
                
                if 'neo4j_law_ids' in payload:
                    law_ids.update(payload['neo4j_law_ids'])
                if 'laws_structured' in payload:
                    for law in payload['laws_structured']:
                        if isinstance(law, dict) and 'law_id' in law:
                            law_ids.add(law['law_id'])
                if 'parent_law_id' in payload:
                    law_ids.add(payload['parent_law_id'])
                if 'law_id' in payload:
                    law_ids.add(payload['law_id'])
        
        return [lid for lid in law_ids if lid]
    
    def _neo4j_enrich(self, law_ids: List[str]) -> Dict[str, Dict]:
        """Get graph context from Neo4j"""
        
        if not law_ids:
            return {}
        
        with self.neo4j.session() as session:
            result = session.run("""
                UNWIND $law_ids as law_id
                MATCH (law:CanonicalLaw {law_id: law_id})
                
                OPTIONAL MATCH (abroga_src)-[abroga:ABROGA]->(law)
                OPTIONAL MATCH (abroga_law:CanonicalLaw)-[:HAS_ARTICLE]->(abroga_src)
                WHERE abroga_src:Article
                
                OPTIONAL MATCH (modif_src)-[modif:MODIFICA]->(law)
                OPTIONAL MATCH (modif_law:CanonicalLaw)-[:HAS_ARTICLE]->(modif_src)
                WHERE modif_src:Article
                
                OPTIONAL MATCH (law)<-[:GOVERNED_BY]-(risk:RiskCategory)
                OPTIONAL MATCH (law)<-[:REGULATED_BY]-(prof:ProfessionalTag)
                OPTIONAL MATCH (doc:TechnicalDocument)-[:CITES]->(law)
                
                RETURN 
                    law.law_id as law_id,
                    law.titolo as titolo,
                    law.tipo as tipo,
                    law.numero as numero,
                    law.anno as anno,
                    law.is_vigente as is_vigente,
                    count(DISTINCT abroga) > 0 as is_abrogated,
                    collect(DISTINCT COALESCE(abroga_law.law_id, abroga_src.law_id))[0..3] as abrogated_by,
                    count(DISTINCT modif) > 0 as is_modified,
                    collect(DISTINCT COALESCE(modif_law.law_id, modif_src.law_id))[0..3] as modified_by,
                    collect(DISTINCT risk.name)[0..5] as risks,
                    collect(DISTINCT prof.categoria)[0..5] as professions,
                    count(DISTINCT doc) as cited_by_docs
            """, law_ids=law_ids)
            
            context = {}
            for row in result:
                law_id = row['law_id']
                if law_id:
                    # ✅ FILTER GARBAGE CATEGORIES
                    risks = [r for r in row['risks'] if r and not is_garbage_category(r)]
                    professions = [p for p in row['professions'] if p and not is_garbage_category(p)]
                    
                    context[law_id] = {
                        'law_id': law_id,
                        'titolo': row['titolo'],
                        'tipo': row['tipo'],
                        'numero': row['numero'],
                        'anno': row['anno'],
                        'is_vigente': row['is_vigente'],
                        'is_abrogated': row['is_abrogated'],
                        'abrogated_by': [x for x in row['abrogated_by'] if x],
                        'is_modified': row['is_modified'],
                        'modified_by': [x for x in row['modified_by'] if x],
                        'risks': risks,
                        'professions': professions,
                        'cited_by_docs': row['cited_by_docs']
                    }
            
            return context
    
    def _process_results(self, qdrant_results, graph_context, filter_abrogated) -> List[Dict]:
        """Process results"""
        
        all_results = []
        
        for collection, hits in qdrant_results.items():
            for hit in hits:
                payload = hit.get('payload', {})
                score = hit.get('score', 0.0)
                
                if not payload:
                    continue
                
                # Extract law_ids
                result_law_ids = set()
                if 'neo4j_law_ids' in payload:
                    result_law_ids.update(payload['neo4j_law_ids'])
                if 'laws_structured' in payload:
                    for law in payload['laws_structured']:
                        if isinstance(law, dict) and 'law_id' in law:
                            result_law_ids.add(law['law_id'])
                if 'parent_law_id' in payload:
                    result_law_ids.add(payload['parent_law_id'])
                if 'law_id' in payload:
                    result_law_ids.add(payload['law_id'])
                
                result_law_ids = {lid for lid in result_law_ids if lid}
                
                # Build result
                result = {
                    'score': score,
                    'collection': collection,
                    'payload': payload,
                    'law_ids': list(result_law_ids),
                    'warnings': [],
                    'graph_context': {}
                }
                
                # Add graph context
                for law_id in result_law_ids:
                    if law_id in graph_context:
                        ctx = graph_context[law_id]
                        result['graph_context'][law_id] = ctx
                        
                        if ctx['is_abrogated']:
                            abrogators = ', '.join(ctx['abrogated_by'][:3]) or 'unknown'
                            result['warnings'].append(f"⚠️ {law_id} ABROGATA da: {abrogators}")
                        elif ctx['is_modified']:
                            modifiers = ', '.join(ctx['modified_by'][:3]) or 'unknown'
                            result['warnings'].append(f"ℹ️ {law_id} modificata da: {modifiers}")
                
                # Filter abrogated
                if filter_abrogated and result_law_ids:
                    abrogated_count = sum(
                        1 for lid in result_law_ids
                        if graph_context.get(lid, {}).get('is_abrogated', False)
                    )
                    if abrogated_count == len(result_law_ids) and abrogated_count > 0:
                        continue
                
                all_results.append(result)
        
        return all_results
    
    def _rerank_v3(self, results: List[Dict], query: str, query_keywords: Set[str]) -> List[Dict]:
        """
        IMPROVED Reranking v3:
        - Heavy weight on keyword matching
        - Moderate boost for graph context
        - Light boost for citation count
        """
        
        for result in results:
            base_score = result['score']
            
            # ============================================
            # 1. KEYWORD RELEVANCE (most important!)
            # ============================================
            
            # Check title
            title = result['payload'].get('titolo') or result['payload'].get('title', '')
            title_relevance = calculate_keyword_relevance(str(title), query_keywords)
            
            # Check text/content
            text = result['payload'].get('text', result['payload'].get('article_text', ''))
            text_relevance = calculate_keyword_relevance(str(text)[:2000], query_keywords)
            
            # Check category
            category = result['payload'].get('categoria_principale', '')
            category_relevance = calculate_keyword_relevance(str(category), query_keywords)
            
            # Combined keyword score (title is most important)
            keyword_score = (title_relevance * 3 + text_relevance * 1 + category_relevance * 2) / 6
            
            # ============================================
            # 2. GRAPH CONTEXT BOOST (moderate)
            # ============================================
            
            graph_boost = 1.0
            
            # Small boost if has graph context
            if result['graph_context']:
                graph_boost *= 1.02
            
            # Small boost for vigente laws
            vigente_count = sum(
                1 for ctx in result['graph_context'].values()
                if ctx.get('is_vigente')
            )
            if vigente_count > 0:
                graph_boost *= 1.02
            
            # Check if risks match query keywords
            for ctx in result['graph_context'].values():
                for risk in ctx.get('risks', []):
                    if any(kw in risk.lower() for kw in query_keywords):
                        graph_boost *= 1.05
                        break
            
            # ============================================
            # 3. PENALTIES
            # ============================================
            
            penalty = 1.0
            
            # Penalty for warnings (abrogated/modified)
            warning_count = len(result['warnings'])
            if warning_count > 0:
                penalty *= (0.92 ** warning_count)
            
            # ============================================
            # 4. COLLECTION BOOST (light)
            # ============================================
            
            collection_boost = 1.0
            if result['collection'] == 'technical_documents_text':
                collection_boost = 1.05
            elif result['collection'] == 'legal_articles':
                collection_boost = 1.02
            
            # ============================================
            # FINAL SCORE
            # ============================================
            
            # Formula: base * (1 + keyword_score) * graph_boost * penalty * collection_boost
            final_score = base_score * (1 + keyword_score * 0.5) * graph_boost * penalty * collection_boost
            
            result['final_score'] = final_score
            result['debug_scores'] = {
                'base': base_score,
                'keyword_score': keyword_score,
                'title_relevance': title_relevance,
                'text_relevance': text_relevance,
                'graph_boost': graph_boost,
                'penalty': penalty,
                'collection_boost': collection_boost
            }
        
        # Sort by final score
        results.sort(key=lambda x: x['final_score'], reverse=True)
        
        return results
    
    def print_results(self, results: List[Dict], max_results: int = 5, show_debug: bool = False):
        """Pretty print results"""
        
        print(f"\n{'='*80}")
        print(f"📊 TOP {min(len(results), max_results)} RESULTS")
        print(f"{'='*80}\n")
        
        if not results:
            print("   ❌ No results found!\n")
            return
        
        for i, result in enumerate(results[:max_results], 1):
            print(f"{'─'*80}")
            print(f"#{i} | Score: {result['final_score']:.4f} (vector: {result['score']:.4f})")
            
            if show_debug and 'debug_scores' in result:
                ds = result['debug_scores']
                print(f"   📊 Debug: keyword={ds['keyword_score']:.2f}, title_rel={ds['title_relevance']:.2f}, graph={ds['graph_boost']:.2f}")
            
            print(f"   Collection: {result['collection']}")
            
            payload = result['payload']
            
            # Title
            title = payload.get('titolo') or payload.get('title') or payload.get('doc_id', 'N/A')
            print(f"   📄 {str(title)[:100]}")
            
            # Source & Category
            if 'source' in payload:
                print(f"   📁 Source: {payload['source']}")
            
            category = payload.get('categoria_principale', '')
            if category and not is_garbage_category(category):
                print(f"   🏷️ Categoria: {category}")
            
            # Warnings
            if result['warnings']:
                print(f"   ⚠️ WARNINGS:")
                for warning in result['warnings'][:2]:
                    print(f"      {warning}")
            
            # Graph context (filtered)
            if result['graph_context']:
                # Show only relevant laws (with matching risks)
                relevant_laws = []
                for law_id, ctx in result['graph_context'].items():
                    if ctx.get('risks') or ctx.get('titolo'):
                        relevant_laws.append((law_id, ctx))
                
                if relevant_laws:
                    print(f"   🔗 Key laws ({len(relevant_laws)}):")
                    for law_id, ctx in relevant_laws[:3]:
                        status = ""
                        if ctx['is_abrogated']:
                            status = " [❌]"
                        elif ctx.get('is_vigente'):
                            status = " [✅]"
                        
                        print(f"      • {law_id}{status}")
                        if ctx.get('titolo'):
                            print(f"        {str(ctx['titolo'])[:60]}...")
                        if ctx.get('risks'):
                            # Filter garbage from risks too
                            clean_risks = [r for r in ctx['risks'] if not is_garbage_category(r)]
                            if clean_risks:
                                print(f"        🎯 {', '.join(clean_risks[:3])}")
            
            # Text preview
            text = payload.get('text', payload.get('article_text', ''))
            if text:
                preview = str(text)[:150].replace('\n', ' ')
                print(f"   📝 {preview}...")
            
            print()
    
    def get_law_details(self, law_id: str) -> Optional[Dict]:
        """Get details for a specific law"""
        
        with self.neo4j.session() as session:
            result = session.run("""
                MATCH (law:CanonicalLaw {law_id: $law_id})
                OPTIONAL MATCH (law)-[:HAS_ARTICLE]->(art:Article)
                OPTIONAL MATCH (doc:TechnicalDocument)-[:CITES]->(law)
                RETURN 
                    law,
                    collect(DISTINCT art.article_num)[0..20] as articles,
                    collect(DISTINCT doc.doc_id)[0..10] as citing_docs
            """, law_id=law_id)
            
            row = result.single()
            if row:
                law = row['law']
                return {
                    'law_id': law['law_id'],
                    'titolo': law.get('titolo'),
                    'tipo': law.get('tipo'),
                    'numero': law.get('numero'),
                    'anno': law.get('anno'),
                    'is_vigente': law.get('is_vigente'),
                    'articles': [a for a in row['articles'] if a],
                    'citing_docs': [d for d in row['citing_docs'] if d]
                }
            return None
    
    def close(self):
        """Close connections"""
        self.neo4j.close()
        print("✅ Connections closed")


if __name__ == "__main__":
    retriever = HybridRetrieverV3()
    
    results = retriever.retrieve(
        "rischio amianto cantiere bonifica",
        top_k=5
    )
    
    retriever.print_results(results, show_debug=True)
    retriever.close()
