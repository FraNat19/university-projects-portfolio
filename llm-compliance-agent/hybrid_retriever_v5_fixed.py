#!/usr/bin/env python3

"""

🎯 HYBRID RETRIEVER v5 FIXED - GraphRAG per Workplace Safety

FIX: Estrae leggi dalla QUERY e le cerca SEMPRE nel grafo!



Struttura Layer 2:

- ABROGA: Article → CanonicalLaw (o Article)

- MODIFICA: Article → Article

- RICHIAMA: Article → CanonicalLaw (o Article)

"""



import requests

import numpy as np

import re

from neo4j import GraphDatabase

from sentence_transformers import SentenceTransformer

from typing import List, Dict, Any, Optional, Set, Tuple

import os

import time



# ============================================

# CONFIG

# ============================================



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



CATEGORY_KEYWORDS = {

    "amianto": ["amianto", "asbesto", "fibre", "mesotelioma", "bonifica"],

    "chimico": ["chimico", "cancerogeno", "mutageno", "sostanze", "agenti chimici"],

    "cantiere": ["cantiere", "edilizia", "costruzioni", "ponteggi", "scavi", "caduta"],

    "dpi": ["dpi", "dispositivi", "protezione individuale", "respiratori", "maschere"],

    "incendio": ["incendio", "antincendio", "emergenza", "evacuazione", "minicodice"],

    "biologico": ["biologico", "virus", "batteri", "covid", "infezione"],

}





def is_garbage_category(category: str) -> bool:

    if not category or not isinstance(category, str):

        return True

    if category in VALID_CATEGORIES:

        return False

    if len(category) < 4:

        return True

    if not re.search(r'[aeiouàèéìòù]', category.lower()):

        return True

    return False





def extract_query_keywords(query: str) -> Set[str]:

    query_lower = query.lower()

    stop_words = {'il', 'lo', 'la', 'i', 'gli', 'le', 'un', 'una', 'di', 'a', 'da', 

                  'in', 'con', 'su', 'per', 'e', 'o', 'che', 'del', 'della', 'dei'}

    words = set(re.findall(r'\b\w{3,}\b', query_lower))

    keywords = words - stop_words

    for category, terms in CATEGORY_KEYWORDS.items():

        for term in terms:

            if term in query_lower:

                keywords.add(category)

    return keywords





def calculate_keyword_relevance(text: str, query_keywords: Set[str]) -> float:

    if not text or not query_keywords:

        return 0.0

    text_lower = text.lower()

    matches = sum(1 for kw in query_keywords if kw in text_lower)

    return min(matches / len(query_keywords), 2.0)





# ============================================

# 🆕 FUNZIONE CRITICA: Estrai leggi dalla QUERY

# ============================================



def extract_laws_from_query(query: str) -> List[Dict]:

    """

    Estrae riferimenti a leggi dalla query dell'utente.

    Ritorna lista di dict con law_id e info.

    """

    laws = []

    query_lower = query.lower()

    

    # Pattern per D.Lgs., D.L., DPR, Legge, DM

    patterns = [

        # D.Lgs. 81/2008, d.lgs 81/08, dlgs 81-2008

        (r'd\.?lgs\.?\s*n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'dlgs'),

        (r'decreto\s+legislativo\s+n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'dlgs'),

        

        # DPR 177/2011

        (r'd\.?p\.?r\.?\s*n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'dpr'),

        (r'decreto\s+(?:del\s+)?presidente\s+(?:della\s+)?repubblica\s+n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'dpr'),

        

        # Legge 257/1992

        (r'legge\s+n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'legge'),

        (r'l\.\s*n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'legge'),

        

        # D.L. 78/2010

        (r'd\.?l\.?\s+n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'dl'),

        (r'decreto\s+legge\s+n?\.?\s*(\d+)\s*[/\-]\s*(\d{2,4})', 'dl'),

        

        # D.M. 10/03/1998

        (r'd\.?m\.?\s+(?:del\s+)?(\d{1,2})\s*[/\-]\s*(\d{1,2})\s*[/\-]\s*(\d{2,4})', 'dm_date'),

    ]

    

    seen = set()

    

    for pattern, tipo in patterns:

        for match in re.finditer(pattern, query_lower, re.IGNORECASE):

            if tipo == 'dm_date':

                # DM con data

                day, month, year = match.groups()

                year = year if len(year) == 4 else ('19' + year if int(year) > 50 else '20' + year)

                law_id = f"dm-{day}-{month}-{year}"

            else:

                numero, anno = match.groups()

                anno = anno if len(anno) == 4 else ('19' + anno if int(anno) > 50 else '20' + anno)

                law_id = f"{tipo}-{numero}-{anno}"

            

            if law_id not in seen:

                seen.add(law_id)

                laws.append({

                    'law_id': law_id,

                    'tipo': tipo,

                    'original': match.group(0),

                    'from_query': True

                })

    

    return laws





class HybridRetrieverV5:

    

    def __init__(

        self,

        qdrant_host: str = "127.0.0.1",

        qdrant_port: int = None,

        neo4j_uri: str = None,

        neo4j_user: str = "neo4j",

        neo4j_password: str = os.getenv("NEO4J_PASSWORD", ""),

        model_name: str = "intfloat/multilingual-e5-large",

        debug: bool = False

    ):

        self.debug = debug

        

        if qdrant_port is None:

            qdrant_port = int(os.getenv("QDRANT_PORT", "6333"))

        if neo4j_uri is None:

            neo4j_port = os.getenv("NEO4J_PORT", "7687")

            neo4j_uri = f"bolt://127.0.0.1:{neo4j_port}"

        

        print("🔧 Initializing Hybrid Retriever v5 (FIXED)...")

        

        self.qdrant_url = f"http://{qdrant_host}:{qdrant_port}"

        print(f"   Qdrant: {self.qdrant_url}")

        

        print(f"   Neo4j: {neo4j_uri}")

        self.neo4j = GraphDatabase.driver(neo4j_uri, auth=(neo4j_user, neo4j_password))

        

        self._test_connections()

        

        print(f"   Loading model: {model_name}")

        self.encoder = SentenceTransformer(model_name)

        

        print("✅ Hybrid Retriever v5 FIXED ready!\n")

    

    def _test_connections(self):

        try:

            response = requests.get(f"{self.qdrant_url}/collections", timeout=5)

            if response.status_code == 200:

                data = response.json()

                collections = data.get('result', {}).get('collections', [])

                print(f"      ✅ Qdrant OK ({len(collections)} collections)")

            

            with self.neo4j.session() as session:

                result = session.run("MATCH (n:CanonicalLaw) RETURN count(n) as c")

                cnt = result.single()['c']

                print(f"      ✅ Neo4j OK ({cnt:,} laws)")

        except Exception as e:

            print(f"      ❌ Connection failed: {str(e)}")

            raise

    

    def retrieve(

        self,

        query: str,

        top_k: int = 10,

        collections: List[str] = None,

        include_graph_context: bool = True

    ) -> List[Dict[str, Any]]:

        

        print(f"\n{'='*70}")

        print(f"🔍 QUERY: {query}")

        print(f"{'='*70}\n")

        

        query_keywords = extract_query_keywords(query)

        print(f"📌 Keywords: {query_keywords}\n")

        

        # 🆕 STEP 0: Estrai leggi dalla QUERY

        query_laws = extract_laws_from_query(query)

        query_law_ids = {law['law_id'] for law in query_laws}

        if query_laws:

            print(f"⚖️ Leggi nella query: {[l['law_id'] for l in query_laws]}")

        

        # 1. Encode

        print("1️⃣ Encoding query...")

        query_embedding = self.encoder.encode(query, normalize_embeddings=True)

        

        # 2. Vector search

        print("2️⃣ Searching Qdrant...")

        qdrant_results = self._qdrant_search(query_embedding, top_k * 3, collections)

        

        # 3. Extract law IDs from documents

        print("3️⃣ Extracting law references from documents...")

        doc_law_ids = self._extract_law_ids(qdrant_results)

        print(f"   Found {len(doc_law_ids)} law_ids in documents")

        

        # 🆕 MERGE: leggi dalla query + leggi dai documenti

        all_law_ids = doc_law_ids | query_law_ids

        print(f"   Total unique law_ids to check: {len(all_law_ids)}")

        

        # 4. Graph enrichment

        graph_context = {}

        query_law_context = {}  # Contesto specifico per leggi nella query

        

        if include_graph_context and all_law_ids:

            print("4️⃣ Enriching with Neo4j graph...")

            start = time.time()

            graph_context = self._neo4j_enrich_fast(list(all_law_ids)[:50])

            print(f"   Got context for {len(graph_context)} laws ({time.time()-start:.2f}s)")

            

            # 🆕 Estrai contesto specifico per leggi nella query

            for law_id in query_law_ids:

                if law_id in graph_context:

                    query_law_context[law_id] = graph_context[law_id]

                    ctx = graph_context[law_id]

                    

                    # Log importante!

                    status = "✅ VIGENTE" if ctx.get('is_vigente') else ""

                    if ctx.get('is_abrogated'):

                        status = f"⛔ ABROGATA da {ctx.get('abrogated_by', ['?'])}"

                    elif ctx.get('is_modified'):

                        status = f"📝 MODIFICATA da {ctx.get('modified_by', ['?'])}"

                    

                    print(f"   📋 {law_id}: {status}")

        

        # 5. Process & Rerank

        print("5️⃣ Processing and reranking...")

        results = self._process_and_rerank(

            qdrant_results, 

            graph_context, 

            query_keywords,

            query_law_context  # 🆕 Passa anche il contesto delle leggi nella query

        )

        

        final = results[:top_k]

        print(f"\n✅ Returning top {len(final)} results\n")

        

        # 🆕 Se la query contiene leggi, aggiungi info in cima ai risultati

        if query_law_context:

            # Inietta warning per leggi abrogate/modificate dalla query

            for result in final:

                for law_id, ctx in query_law_context.items():

                    if ctx.get('is_abrogated'):

                        abr = ', '.join(ctx.get('abrogated_by', ['?'])[:2])

                        warning = f"⛔ ATTENZIONE: {law_id} è ABROGATA da {abr}!"

                        if warning not in result['warnings']:

                            result['warnings'].insert(0, warning)

                    elif ctx.get('is_modified'):

                        mod = ', '.join(ctx.get('modified_by', ['?'])[:2])

                        warning = f"📝 NOTA: {law_id} è stata MODIFICATA da {mod}"

                        if warning not in result['warnings']:

                            result['warnings'].insert(0, warning)

        

        return final

    

    def _qdrant_search(self, embedding: np.ndarray, top_k: int, collections: List[str] = None) -> Dict[str, List]:

        if collections is None:

            collections = ['technical_documents_text', 'legal_context', 'legal_articles']

        

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

                    hits = response.json().get('result', [])

                    results[collection] = hits

                    print(f"   {collection}: {len(hits)} results")

                else:

                    results[collection] = []

            except Exception as e:

                print(f"   ⚠️ {collection}: {str(e)[:50]}")

                results[collection] = []

        

        return results

    

    def _extract_law_ids(self, qdrant_results: Dict[str, List]) -> Set[str]:

        law_ids = set()

        for collection, hits in qdrant_results.items():

            for hit in hits:

                payload = hit.get('payload', {})

                # From laws_structured

                for law in payload.get('laws_structured', []):

                    if isinstance(law, dict) and 'law_id' in law:

                        law_ids.add(law['law_id'])

                # Direct fields

                if 'law_id' in payload:

                    law_ids.add(payload['law_id'])

                if 'parent_law_id' in payload:

                    law_ids.add(payload['parent_law_id'])

        return {lid for lid in law_ids if lid}

    

    def _neo4j_enrich_fast(self, law_ids: List[str]) -> Dict[str, Dict]:

        """

        ✅ QUERY OTTIMIZZATA - Esegue 3 query separate

        """

        if not law_ids:

            return {}

        

        context = {}

        

        with self.neo4j.session() as session:

            # QUERY 1: Info base delle leggi

            result = session.run("""

                UNWIND $law_ids as law_id

                MATCH (law:CanonicalLaw {law_id: law_id})

                RETURN 

                    law.law_id as law_id,

                    law.titolo as titolo,

                    law.is_vigente as is_vigente,

                    law.abrogato_da as abrogato_da

            """, law_ids=law_ids)

            

            for row in result:

                context[row['law_id']] = {

                    'law_id': row['law_id'],

                    'titolo': row['titolo'],

                    'is_vigente': row['is_vigente'],

                    'is_abrogated': row['is_vigente'] == False or row['abrogato_da'] is not None,

                    'abrogated_by': [row['abrogato_da']] if row['abrogato_da'] else [],

                    'is_modified': False,

                    'modified_by': []

                }

            

            # QUERY 2: Check ABROGA (target = CanonicalLaw)

            result = session.run("""

                UNWIND $law_ids as law_id

                MATCH (law:CanonicalLaw {law_id: law_id})

                MATCH (src)-[:ABROGA]->(law)

                WITH law.law_id as law_id, collect(DISTINCT

                    CASE 

                        WHEN src:Article THEN src.parent_law_id

                        WHEN src:CanonicalLaw THEN src.law_id

                        ELSE null

                    END

                )[0..3] as abrogators

                RETURN law_id, abrogators

            """, law_ids=law_ids)

            

            for row in result:

                if row['law_id'] in context:

                    abr = [a for a in row['abrogators'] if a]

                    if abr:

                        context[row['law_id']]['is_abrogated'] = True

                        context[row['law_id']]['abrogated_by'] = abr

            

            # QUERY 3: Check MODIFICA (target = CanonicalLaw o Article della law)

            result = session.run("""

                UNWIND $law_ids as law_id

                MATCH (law:CanonicalLaw {law_id: law_id})

                OPTIONAL MATCH (src)-[:MODIFICA]->(law)

                OPTIONAL MATCH (src2)-[:MODIFICA]->(art:Article)<-[:HAS_ARTICLE]-(law)

                WITH law.law_id as law_id, 

                     collect(DISTINCT CASE WHEN src:Article THEN src.parent_law_id ELSE src.law_id END) +

                     collect(DISTINCT CASE WHEN src2:Article THEN src2.parent_law_id ELSE src2.law_id END) as all_modifiers

                UNWIND all_modifiers as modifier

                WITH law_id, collect(DISTINCT modifier)[0..3] as modifiers

                WHERE size(modifiers) > 0

                RETURN law_id, modifiers

            """, law_ids=law_ids)

            

            for row in result:

                if row['law_id'] in context:

                    mods = [m for m in row['modifiers'] if m]

                    if mods:

                        context[row['law_id']]['is_modified'] = True

                        context[row['law_id']]['modified_by'] = mods

        

        return context

    

    def _process_and_rerank(

        self, 

        qdrant_results: Dict[str, List], 

        graph_context: Dict[str, Dict],

        query_keywords: Set[str],

        query_law_context: Dict[str, Dict] = None  # 🆕

    ) -> List[Dict]:

        

        all_results = []

        seen_docs = set()

        

        for collection, hits in qdrant_results.items():

            for hit in hits:

                payload = hit.get('payload', {})

                score = hit.get('score', 0.0)

                

                if not payload:

                    continue

                

                # Deduplicate

                doc_id = payload.get('doc_id', payload.get('law_id', str(hit.get('id', ''))))

                if collection == 'technical_documents_text' and doc_id in seen_docs:

                    continue

                seen_docs.add(doc_id)

                

                # Extract law_ids for this result

                result_law_ids = set()

                for law in payload.get('laws_structured', []):

                    if isinstance(law, dict) and 'law_id' in law:

                        result_law_ids.add(law['law_id'])

                if 'law_id' in payload:

                    result_law_ids.add(payload['law_id'])

                

                # Build result

                result = {

                    'score': score,

                    'collection': collection,

                    'payload': payload,

                    'law_ids': list(result_law_ids),

                    'warnings': [],

                    'graph_context': {}

                }

                

                # Add graph context and warnings from documents

                for law_id in result_law_ids:

                    if law_id in graph_context:

                        ctx = graph_context[law_id]

                        result['graph_context'][law_id] = ctx

                        

                        if ctx['is_abrogated']:

                            abr = ', '.join(ctx['abrogated_by'][:2]) if ctx['abrogated_by'] else '?'

                            result['warnings'].append(f"⚠️ {law_id} ABROGATA da: {abr}")

                        elif ctx['is_modified']:

                            mod = ', '.join(ctx['modified_by'][:2]) if ctx['modified_by'] else '?'

                            result['warnings'].append(f"📝 {law_id} modificata da: {mod}")

                

                # Calculate final score

                base_score = score

                

                # Keyword boost

                title = str(payload.get('titolo') or payload.get('title') or '')

                text = str(payload.get('text') or payload.get('sintesi') or '')[:1000]

                keyword_score = (

                    calculate_keyword_relevance(title, query_keywords) * 2 +

                    calculate_keyword_relevance(text, query_keywords)

                ) / 3

                

                # Graph boost

                graph_boost = 1.0

                if result['graph_context']:

                    vigente = sum(1 for c in result['graph_context'].values() if c.get('is_vigente'))

                    graph_boost = 1.0 + (vigente * 0.02)

                

                # Penalty for abrogated

                penalty = 1.0

                abrogated = sum(1 for c in result['graph_context'].values() if c.get('is_abrogated'))

                if abrogated > 0:

                    penalty = 0.9 ** abrogated

                

                # Collection boost

                coll_boost = 1.05 if collection == 'technical_documents_text' else 1.0

                

                result['final_score'] = base_score * (1 + keyword_score * 0.3) * graph_boost * penalty * coll_boost

                

                all_results.append(result)

        

        all_results.sort(key=lambda x: x['final_score'], reverse=True)

        return all_results

    

    def print_results(self, results: List[Dict], max_results: int = 5, show_debug: bool = False):

        print(f"\n{'='*70}")

        print(f"📊 TOP {min(len(results), max_results)} RESULTS")

        print(f"{'='*70}\n")

        

        if not results:

            print("   ❌ No results found!\n")

            return

        

        for i, result in enumerate(results[:max_results], 1):

            print(f"{'─'*70}")

            print(f"#{i} | Score: {result['final_score']:.4f} (vector: {result['score']:.4f})")

            print(f"   📁 Collection: {result['collection']}")

            

            payload = result['payload']

            title = payload.get('titolo') or payload.get('title') or payload.get('doc_id', 'N/A')

            print(f"   📄 {str(title)[:80]}")

            

            source = payload.get('source', '')

            if source:

                print(f"   🏢 Source: {source}")

            

            category = payload.get('categoria_principale', '')

            if category and not is_garbage_category(category):

                print(f"   🏷️ Categoria: {category}")

            

            if result['law_ids']:

                print(f"   ⚖️ Laws ({len(result['law_ids'])}): {', '.join(result['law_ids'][:5])}")

            

            # Warnings

            if result['warnings']:

                for w in result['warnings'][:3]:

                    print(f"   {w}")

            

            # Graph context summary

            if result['graph_context']:

                vigenti = sum(1 for c in result['graph_context'].values() if c.get('is_vigente'))

                abrogati = sum(1 for c in result['graph_context'].values() if c.get('is_abrogated'))

                modificati = sum(1 for c in result['graph_context'].values() if c.get('is_modified'))

                print(f"   🔗 Graph: {vigenti} vigenti, {modificati} modificate, {abrogati} abrogate")

            

            # Text preview

            text = payload.get('text') or payload.get('sintesi') or ''

            if text:

                preview = str(text)[:150].replace('\n', ' ')

                print(f"   📝 {preview}...")

            

            print()

    

    def close(self):

        self.neo4j.close()

        print("✅ Connections closed")





if __name__ == "__main__":

    import argparse

    

    parser = argparse.ArgumentParser()

    parser.add_argument('--query', '-q', type=str, default="Il D.Lgs. 626/1994 è ancora in vigore?")

    parser.add_argument('--top-k', '-k', type=int, default=5)

    parser.add_argument('--no-graph', action='store_true', help='Disable graph enrichment')

    parser.add_argument('--debug', '-d', action='store_true')

    

    args = parser.parse_args()

    

    retriever = HybridRetrieverV5(debug=args.debug)

    

    results = retriever.retrieve(

        args.query, 

        top_k=args.top_k,

        include_graph_context=not args.no_graph

    )

    

    retriever.print_results(results, max_results=args.top_k, show_debug=args.debug)

    retriever.close()
