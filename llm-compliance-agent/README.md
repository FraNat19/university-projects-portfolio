# GraphRAG-based LLM Agent for Workplace Safety Compliance

> Master's Thesis Project | Sapienza University of Rome | 2024-2025

A hybrid Retrieval-Augmented Generation (RAG) system that combines semantic vector search with a legal knowledge graph to provide accurate, legally-validated answers on Italian workplace safety regulations.

## 🎯 Problem Statement

Workplace safety professionals (RSPP, consultants, inspectors) face a critical challenge: navigating thousands of technical documents and legal articles while ensuring cited regulations are still in force. Traditional search systems fail to detect when a law has been **repealed** or **modified**, potentially leading to advice based on obsolete regulations.

**Example:** Many documents still reference *D.Lgs. 626/1994*, which was repealed by *D.Lgs. 81/2008* (Article 304). A standard RAG system would return these documents without warning—our system explicitly flags this.

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                            USER QUERY                               │
│            "Is D.Lgs. 626/1994 still in force?"                     │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                       HYBRID RETRIEVER                              │
│  ┌───────────────────────┐      ┌───────────────────────┐          │
│  │    QDRANT             │      │      NEO4J            │          │
│  │    Vector Search      │      │   Knowledge Graph     │          │
│  │  ┌─────────────────┐  │      │  ┌─────────────────┐  │          │
│  │  │ technical_docs  │  │      │  │ CanonicalLaw    │  │          │
│  │  │ legal_articles  │  │      │  │ Article         │  │          │
│  │  │ legal_context   │  │      │  │ ABROGA/MODIFICA │  │          │
│  │  └─────────────────┘  │      │  └─────────────────┘  │          │
│  └───────────────────────┘      └───────────────────────┘          │
│              │                            │                         │
│              └──────────┬─────────────────┘                         │
│                         ▼                                           │
│              ┌─────────────────────┐                                │
│              │  Context Enrichment │                                │
│              │  + Legal Warnings   │                                │
│              └─────────────────────┘                                │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│                        LLM (Qwen 2.5 72B)                           │
│         Generates grounded response with regulatory status          │
└─────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────┐
│  "No, D.Lgs. 626/1994 was REPEALED by Article 304 of D.Lgs.        │
│   81/2008. The current reference for workplace safety is the        │
│   Testo Unico (D.Lgs. 81/2008)."                                    │
└─────────────────────────────────────────────────────────────────────┘
```

## 🔧 Technical Components

### 1. Data Ingestion Pipeline

| Stage | Description | Tools |
|-------|-------------|-------|
| **PDF Extraction** | Layout-aware parsing preserving tables, figures, structure | Docling (IBM) |
| **Citation Extraction** | Regex-based identification and normalization of legal references | Custom Python |
| **Semantic Enrichment** | Risk categories, professional profiles, keywords | LLM + Clustering |
| **Chunking** | Hierarchical splitting with overlap | LangChain-style |
| **Embedding** | Multilingual dense vectors (1024 dim) | E5-large |

### 2. Vector Database (Qdrant)

Three collections for granular retrieval:

- **`technical_documents_text`**: 2,000+ INAIL/EU-OSHA chunked documents
- **`legal_articles`**: Individual law articles for precise matching
- **`legal_context`**: Law-level metadata (title, publication date, URN)

### 3. Knowledge Graph (Neo4j)

**Node Types:**
- `CanonicalLaw`: Represents a complete law (e.g., D.Lgs. 81/2008)
- `Article`: Individual articles within laws
- `TechnicalDocument`: INAIL/OSHA publications

**Relationship Types:**
- `ABROGA`: Law/Article X repeals Law/Article Y
- `MODIFICA`: Law/Article X amends Article Y  
- `RICHIAMA`: Law/Article X references Law/Article Y
- `CITES`: Technical document cites a law/article

**Example Query:**
```cypher
// Find all laws repealed by D.Lgs. 81/2008
MATCH (a:Article {parent_law_id: 'dlgs-81-2008'})-[:ABROGA]->(repealed)
RETURN a.article_num, repealed.law_id
```

### 4. Hybrid Retriever

The `HybridRetrieverV5` implements a 5-stage pipeline:

1. **Query Encoding**: Embed user query with multilingual-e5-large
2. **Law Extraction**: Parse legal references from query text (e.g., "D.Lgs. 626/1994" → `dlgs-626-1994`)
3. **Vector Search**: Retrieve top-k documents from Qdrant collections
4. **Graph Enrichment**: For each cited law, query Neo4j for:
   - Current validity (`is_vigente`)
   - Abrogation relationships
   - Modification history
5. **Reranking**: Score adjustment with validity boost/penalty

**Key Innovation:** The retriever extracts laws from *both* retrieved documents *and* the query itself, ensuring questions like "Is Law X valid?" always trigger a graph lookup.

### 5. RAG Integration

- **LLM**: Qwen 2.5 72B via Ollama
- **Prompt Engineering**: Three specialized modes (Legal, Practical, General) with automatic selection
- **Context Formatting**: Structured presentation with visual warnings for repealed laws

## 📈 Results

Evaluated on 14 benchmark queries covering fire prevention, asbestos, carcinogens, and construction safety:

| Metric | Vector-Only | Hybrid (Vector + Graph) |
|--------|-------------|------------------------|
| Query Time | ~0.1s | ~40s |
| Repealed Laws Detected | ❌ No | ✅ Yes |
| Modification Tracking | ❌ No | ✅ Yes |
| Legal Traceability | Partial | Full |

**Trade-off:** The additional processing time is justified in safety-critical domains where accuracy outweighs speed.

## 🛠️ Tech Stack

| Component | Technology |
|-----------|------------|
| Vector Database | Qdrant |
| Graph Database | Neo4j |
| Embeddings | intfloat/multilingual-e5-large |
| LLM | Qwen 2.5 72B (Ollama) |
| PDF Processing | Docling (IBM) |
| Language | Python 3.10+ |
| Infrastructure | HPC cluster (NVIDIA A100) |


## 📊 Data Sources

| Source | Type | Volume |
|--------|------|--------|
| **INAIL** | Technical publications (risk assessment, guidelines) | 800+ documents |
| **EU-OSHA** | European safety resources | 200+ documents |
| **Normattiva** | Italian legislation (XML format) | 50+ laws, 3,000+ articles |

## 👤 Author

**Francesco Natali**  
Master's Degree in Statistical Methods and Applications (Data Analyst)  
Sapienza University of Rome

*This project was developed as part of my Master's thesis (A.Y. 2024/2025)*
