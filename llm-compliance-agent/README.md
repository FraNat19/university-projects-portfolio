<div align="center">

# ⚖️ GraphRAG for Workplace Safety Compliance

**A retrieval system that knows when the law it just cited has been repealed**

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/Neo4j-4581C3?style=for-the-badge&logo=neo4j&logoColor=white" alt="Neo4j">
<img src="https://img.shields.io/badge/Qdrant-DC244C?style=for-the-badge" alt="Qdrant">
<img src="https://img.shields.io/badge/SLURM-HPC-1a7f37?style=for-the-badge" alt="SLURM">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>Master's thesis · MSc Statistical Methods and Applications · Sapienza University of Rome · 2024/2025</sub>

</div>

---

## 📋 The problem

A safety consultant asks whether a procedure complies with *D.Lgs. 626/1994*. A standard RAG system finds a dozen documents citing that law, summarises them confidently, and is completely wrong — because 626/1994 was repealed by Article 304 of *D.Lgs. 81/2008*.

The documents are not wrong. They were written when the law was in force. **Vector similarity has no way to represent "this was true and now is not."**

Embeddings encode what text means, not whether it still applies. So this system runs two retrievers side by side: a vector store for meaning, and a knowledge graph for legal status.

```
                        USER QUERY
              "Is D.Lgs. 626/1994 still in force?"
                             │
              ┌──────────────┴──────────────┐
              ▼                             ▼
     ┌─────────────────┐          ┌──────────────────┐
     │     QDRANT      │          │      NEO4J       │
     │  what it means  │          │  whether it holds│
     ├─────────────────┤          ├──────────────────┤
     │ technical_docs  │          │ CanonicalLaw     │
     │ legal_articles  │          │ Article          │
     │ legal_context   │          │ ABROGA / MODIFICA│
     │ technical_images│          │ RICHIAMA / CITES │
     └────────┬────────┘          └────────┬─────────┘
              └──────────────┬─────────────┘
                             ▼
                  Context + validity warnings
                             ▼
                    LLM (Qwen 2.5 72B)
                             ▼
     "No. Repealed by Art. 304 of D.Lgs. 81/2008.
      The current reference is the Testo Unico."
```

---

## 🔑 The design decision that matters

Most hybrid retrievers extract entities from the **retrieved documents** and enrich from there. This one also extracts them from the **query itself**.

```python
# hybrid_retriever_v5_fixed.py
# FIX: Estrae leggi dalla QUERY e le cerca SEMPRE nel grafo
```

The difference shows up exactly where it counts. Ask *"is Law X valid?"* and the vector search returns documents about X's subject matter — none of which necessarily mention X's repeal, because a document about scaffolding safety has no reason to discuss which statute superseded which. Parsing the law reference straight out of the question guarantees the graph lookup happens regardless of what the vector search returned.

Retrieval runs in five stages: encode the query, parse legal references from its text, search the Qdrant collections, look up validity and abrogation history in Neo4j for every law found on either path, then rerank with a boost or penalty based on legal status.

---

## 🗂️ The codebase

Built and run on an HPC cluster with SLURM. The pipeline is layered, and each layer is a separate job because each one takes hours.

### Acquisition and parsing

| File | Does |
|---|---|
| [`inail_hpc_main (1).py`](inail_hpc_main%20\(1\).py) | Scrapes INAIL technical publications |
| [`osha_scraper_hpc_new.py`](osha_scraper_hpc_new.py) | Scrapes EU-OSHA resources |
| [`07_legal_normattiva.py`](07_legal_normattiva.py) | Pulls Italian legislation from Normattiva as Akoma Ntoso XML |
| [`inail_docling (1).py`](inail_docling%20\(1\).py), [`osha_docling (4).py`](osha_docling%20\(4\).py) | Layout-aware PDF parsing with Docling, preserving tables and structure |
| [`vlm_only_script (3).py`](vlm_only_script%20\(3\).py) | Captions extracted figures with a vision-language model |

### Indexing

| File | Does |
|---|---|
| [`00_create_collections.py`](00_create_collections.py) | Creates the four Qdrant collections |
| [`02_text_chunks_FIXED (1).py`](02_text_chunks_FIXED%20\(1\).py) | Hierarchical chunking with overlap, reading `laws_structured` before falling back to semantic metadata |
| [`05_images_vlm.py`](05_images_vlm.py) | Indexes figures, filtering out logos, covers and low-information images |
| [`tf_idf_embedgemma_hpc (1).py`](tf_idf_embedgemma_hpc%20\(1\).py) | Embedding experiments |
| [`enrichment_inail (1).py`](enrichment_inail%20\(1\).py), [`Enrichment_OSHA.ipynb`](Enrichment_OSHA.ipynb) | Risk categories, professional profiles, keyword extraction |

### Graph construction

Six ordered stages, each aligned on `law_id` as the primary key:

| Stage | File | Does |
|---|---|---|
| Schema | [`neo4j_schema_setup.py`](neo4j_schema_setup.py) | Legal ontology, constraints and indexes, with no APOC dependency so it runs on any Neo4j |
| 0b | [`layer0b_legal_bridge_fixed.py`](layer0b_legal_bridge_fixed.py) | Syncs `laws_structured` out of the document collection, separating real publication dates from placeholders |
| 1 | [`layer1_qdrant_bridge_FIXED_V2.py`](layer1_qdrant_bridge_FIXED_V2.py) | Builds `CITES` edges, creating missing `Article` nodes on the fly |
| 1b | [`layer1b_normattiva_sync.py`](layer1b_normattiva_sync.py) | Enriches existing nodes with the official Normattiva metadata |
| 2 | [`layer2_fixed.py`](layer2_fixed.py) | Extracts `ABROGA` and `MODIFICA` relations — the ones that carry legal status |
| 3 | [`layer3_v2_aligned.py`](layer3_v2_aligned.py) | Semantic bridge, using the same category taxonomy as the enrichment step |

### Retrieval and generation

| File | Does |
|---|---|
| [`hybrid_retriever_v5_fixed.py`](hybrid_retriever_v5_fixed.py) | The five-stage hybrid retriever |
| [`graphrag_rag_v3.py`](graphrag_rag_v3.py) | Full RAG loop with three prompt modes — legal, practical, general — selected automatically |
| [`test_retriever.py`](test_retriever.py) | Benchmark queries across fire prevention, asbestos, carcinogens and construction |

### Operations

SLURM jobs ([`slurm_graphrag_full.sbatch`](slurm_graphrag_full.sbatch), [`slurm_stack_qdrant_neo4j.sbatch`](slurm_stack_qdrant_neo4j.sbatch), the scraping jobs), stack control ([`start_stack.sh`](start_stack.sh), [`status_stack.sh`](status_stack.sh), [`tunnel_stack.sh`](tunnel_stack.sh)) and health checks ([`check_connections.py`](check_connections.py), [`verify_neo4j_status.py`](verify_neo4j_status.py), [`check_cites.py`](check_cites.py)).

---

## 🕸️ Graph model

**Nodes** — `CanonicalLaw` (a whole statute), `Article` (an article within one), `TechnicalDocument` (an INAIL or EU-OSHA publication).

**Relations** — `ABROGA` (X repeals Y), `MODIFICA` (X amends Y), `RICHIAMA` (X references Y), `CITES` (a document cites a law or article).

```cypher
// Everything repealed by D.Lgs. 81/2008
MATCH (a:Article {parent_law_id: 'dlgs-81-2008'})-[:ABROGA]->(repealed)
RETURN a.article_num, repealed.law_id
```

Deduplication was the hard part. An early version keyed citations by `law_id`, which collapsed every reference to the same statute into one edge and destroyed article-level precision. The current version keys on `full_citation`, so "Art. 71 comma 4 of D.Lgs. 81/2008" stays distinct from "Art. 26 of D.Lgs. 81/2008". [`reenrich_laws_v3_FINAL.py`](reenrich_laws_v3_FINAL.py) exists to repair graphs built under the old assumption.

---

## 📊 Data and stack

| Source | Type | Volume |
|---|---|---|
| INAIL | Technical publications, risk assessment, guidelines | 800+ documents |
| EU-OSHA | European safety resources | 200+ documents |
| Normattiva | Italian legislation, Akoma Ntoso XML | 50+ laws, 3,000+ articles |

| Component | Choice |
|---|---|
| Vector store | Qdrant — `technical_documents_text`, `legal_articles`, `legal_context`, `technical_images` |
| Graph store | Neo4j |
| Embeddings | `intfloat/multilingual-e5-large`, 1024 dimensions |
| Generation | Qwen 2.5 72B via Ollama |
| PDF parsing | Docling |
| Infrastructure | HPC cluster, NVIDIA A100, SLURM |

---

## 📈 Results

Evaluated on 14 benchmark queries covering fire prevention, asbestos, carcinogens and construction safety.

| | Vector only | Hybrid (vector + graph) |
|---|---|---|
| Query time | ~0.1 s | ~40 s |
| Repealed laws detected | ✗ | ✓ |
| Modification tracking | ✗ | ✓ |
| Legal traceability | Partial | Full |

**Four hundred times slower**, and worth it here. A consultant waiting forty seconds for a correct answer is fine; a consultant citing a repealed statute in a safety assessment is a liability. The trade only makes sense because this is a domain where being wrong is expensive — it would be the wrong call for a search box.

---

## ⚙️ Running it

```bash
cp .env.example .env    # fill in HUGGINGFACE_TOKEN and NEO4J_PASSWORD
set -a && source .env && set +a
```

No credentials are stored in the code. The scripts read `HUGGINGFACE_TOKEN` and `NEO4J_PASSWORD` from the environment, and the SLURM jobs refuse to start if `NEO4J_PASSWORD` is unset rather than falling back to a default.

**What will not run as-is.** Paths point at the cluster where this was built (`/mnt/beegfs/proj/dss.dmaia/...`) and need repointing. The corpora are not redistributed: INAIL and EU-OSHA publications are third-party documents, and the Normattiva XML is bulk-downloaded separately. `inail_docling (1).py`, `osha_docling (4).py` and the enrichment notebook are Colab exports, so they carry `!pip install` lines that are not valid Python outside a notebook.

**On the multiple versions.** Several components appear more than once — `hybrid_retriever` at v1, v3 and v5, `layer1_qdrant_bridge` as FIXED and FIXED_V2, `layer3` as production and v2_aligned. The later file is the current one in each case. The earlier ones are kept because the fixes between them are the substance of the work: the deduplication key, the query-side law extraction, the alignment of category taxonomies between enrichment and the semantic bridge. This folder is a working thesis codebase rather than a packaged release, and it is being tidied incrementally.

[`Thesis_Presentation.pdf`](Thesis_Presentation.pdf) covers the architecture and results as presented.

---

## 👤 Author

**Francesco Natali** · MSc in Statistical Methods and Applications, Sapienza University of Rome
