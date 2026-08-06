
#!/usr/bin/env python3

"""

🎯 QDRANT INDEXING - Configuration

"""

from pathlib import Path



BASE_DIR = Path("/mnt/beegfs/proj/dss.dmaia")

THESIS_DIR = BASE_DIR / "thesis_graph_rag"

INAIL_DIR = BASE_DIR / "INAILProj"

OSHA_DIR = BASE_DIR / "OSHAProj"

NORMATTIVA_DIR = THESIS_DIR / "data" / "normattiva_filtered"



INAIL_ENRICHED = INAIL_DIR / "enriched_json"

INAIL_IMAGES = INAIL_DIR / "images"

INAIL_TABLES = INAIL_DIR / "tables"

INAIL_ENRICHED_TABLES = INAIL_DIR / "enriched_tables"



OSHA_ENRICHED = OSHA_DIR / "enriched_json"

OSHA_IMAGES = OSHA_DIR / "images"

OSHA_TABLES = OSHA_DIR / "tables"

OSHA_ENRICHED_TABLES = OSHA_DIR / "enriched_tables"



# ✅ FIXED: 127.0.0.1 invece di master1

QDRANT_HOST = "127.0.0.1"

QDRANT_PORT = 6333



NEO4J_URI = "bolt://127.0.0.1:7687"

NEO4J_USER = "neo4j"

NEO4J_PASSWORD = "thesis2024"



MODEL_NAME = "intfloat/multilingual-e5-large"

EMBEDDING_DIM = 1024

CHUNK_SIZE = 1000

OVERLAP = 200

BATCH_SIZE = 100

EMBEDDING_BATCH_SIZE = 32



COLLECTIONS = {

    'technical_text': 'technical_documents_text',

    'technical_images': 'technical_images',

    'technical_tables': 'technical_tables',

    'legal_articles': 'legal_articles',

    'legal_context': 'legal_context'

}



PAYLOAD_INDEXES = {

    'technical_text': [

        ('source', 'keyword'),

        ('categoria_principale', 'keyword'),

        ('data_pubblicazione', 'keyword'),

    ],

    'technical_images': [

        ('source', 'keyword'),

        ('parent_doc_id', 'keyword'),

    ],

    'technical_tables': [

        ('source', 'keyword'),

        ('parent_doc_id', 'keyword'),

    ],

    'legal_articles': [

        ('law_id', 'keyword'),

        ('law_type', 'keyword'),

    ],

    'legal_context': [

        ('law_id', 'keyword'),

        ('anno', 'keyword'),

    ]

}



if __name__ == "__main__":

    print("✅ Configuration module loaded!")

    print(f"   INAIL enriched: {INAIL_ENRICHED}")

    print(f"   OSHA enriched: {OSHA_ENRICHED}")

    print(f"   Normattiva: {NORMATTIVA_DIR}")

