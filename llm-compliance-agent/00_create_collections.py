
#!/usr/bin/env python3

"""

🎯 CREATE QDRANT COLLECTIONS

"""

import sys

import os

sys.path.insert(0, os.path.dirname(__file__))



import importlib.util

spec = importlib.util.spec_from_file_location("config", "01_config.py")

config = importlib.util.module_from_spec(spec)

spec.loader.exec_module(config)



from qdrant_client import QdrantClient

from qdrant_client.models import Distance, VectorParams



client = QdrantClient(host=config.QDRANT_HOST, port=config.QDRANT_PORT)



print("="*80)

print("🔧 CREATING QDRANT COLLECTIONS")

print("="*80)



for key, collection_name in config.COLLECTIONS.items():

    try:

        client.get_collection(collection_name=collection_name)

        print(f"✅ {collection_name} already exists")

    except Exception:

        print(f"⚙️  Creating {collection_name}...")

        client.create_collection(

            collection_name=collection_name,

            vectors_config=VectorParams(size=config.EMBEDDING_DIM, distance=Distance.COSINE)

        )

        print(f"✅ {collection_name} created!")



print("="*80)

print("✅ ALL COLLECTIONS READY!")

print("="*80)

