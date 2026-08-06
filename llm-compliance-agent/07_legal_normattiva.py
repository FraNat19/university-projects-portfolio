#!/usr/bin/env python3

"""

⚖️ LEGAL NORMATTIVA - FIXED VERSION



FIXES:

✅ USA EnhancedLawExtractor per consistenza

✅ Estrae date REALI da XML

✅ Gestisce correttamente URN da XML

"""



import json

import hashlib

import re

import xml.etree.ElementTree as ET

from pathlib import Path

from typing import List, Dict, Optional

from tqdm import tqdm

import sys

import os

import torch



sys.path.insert(0, os.path.dirname(__file__))

import importlib.util

spec = importlib.util.spec_from_file_location("config", "01_config.py")

config = importlib.util.module_from_spec(spec)

spec.loader.exec_module(config)



from qdrant_client import QdrantClient

from qdrant_client.models import PointStruct, Distance, VectorParams

from sentence_transformers import SentenceTransformer



# ✅ IMPORT EXTRACTOR

from law_extraction_enhanced import EnhancedLawExtractor



NORMATTIVA_FILTERED = Path("/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/normattiva_filtered")



NS = {

    'nir': 'http://www.normeinrete.it/nir/2.2/',

    'dsp': 'http://www.normeinrete.it/nir/disposizioni/2.2/',

    'h': 'http://www.w3.org/HTML/1998/html4',

    'xlink': 'http://www.w3.org/1999/xlink'

}





def extract_text_from_element(elem) -> str:

    """Extract text recursively"""

    if elem is None:

        return ""

    

    text_parts = []

    if elem.text:

        text_parts.append(elem.text.strip())

    

    for child in elem:

        text_parts.append(extract_text_from_element(child))

        if child.tail:

            text_parts.append(child.tail.strip())

    

    return " ".join([t for t in text_parts if t])





def determine_law_type_and_id(xml_path: Path, root, extractor: EnhancedLawExtractor) -> tuple:

    """

    ✅ USA EnhancedLawExtractor + estrae DATA REALE da XML

    """

    

    # Estrai da XML

    tipo_doc = root.find('.//nir:tipoDoc', NS)

    num_doc = root.find('.//nir:numDoc', NS)

    data_doc = root.find('.//nir:dataDoc', NS)

    urn_elem = root.find('.//nir:urn', NS)

    

    tipo_text = tipo_doc.text if tipo_doc is not None else ""

    numero = num_doc.text if num_doc is not None else ""

    

    # ✅ ESTRAI DATA REALE (formato: norm="19480331")

    data_norm = data_doc.get('norm', '') if data_doc is not None else ""

    anno = data_norm[:4] if len(data_norm) >= 4 else ""

    

    # Converti data_norm in ISO (19480331 → 1948-03-31)

    data_iso = None

    if len(data_norm) == 8:

        try:

            data_iso = f"{data_norm[:4]}-{data_norm[4:6]}-{data_norm[6:8]}"

        except:

            pass

    

    # ✅ ESTRAI URN REALE

    urn_real = urn_elem.get('valore', '') if urn_elem is not None else ''

    

    # Costruisci citazione sintetica

    synthetic_citation = f"{tipo_text} {numero}/{anno}"

    

    # ✅ USA EXTRACTOR

    laws, conf = extractor.extract_with_confidence(synthetic_citation, max_laws=1)

    

    if laws:

        law = laws[0]

        

        # ✅ SOVRASCRIVI con data e URN reali

        law['data'] = data_iso or f"{anno}-01-01"

        law['data_available'] = data_iso is not None

        

        # URN: usa quello da XML se valido, altrimenti quello generato

        if urn_real and urn_real != 'urn:' and len(urn_real) > 10:

            law['urn'] = urn_real

        

        return (

            law['tipo'],

            law['law_id'],

            law['tipo_short'],

            numero,

            anno,

            law['data'],

            law['urn']

        )

    

    # Fallback se extractor fallisce

    return determine_law_type_fallback(xml_path, root, data_iso)





def determine_law_type_fallback(xml_path: Path, root, data_iso: str) -> tuple:

    """Fallback manuale"""

    

    parent_folder = xml_path.parent.name

    grandparent_folder = xml_path.parent.parent.name

    

    num_doc = root.find('.//nir:numDoc', NS)

    data_doc = root.find('.//nir:dataDoc', NS)

    numero = num_doc.text if num_doc is not None else ""

    anno = data_doc.get('norm', '')[:4] if data_doc is not None else ""

    

    if 'Decreti Legislativi' in grandparent_folder:

        law_type = "decreto.legislativo"

        law_id = f"dlgs-{numero}-{anno}"

        tipo_short = "dlgs"

    elif 'DPR' in grandparent_folder:

        law_type = "decreto.del.presidente.della.repubblica"

        law_id = f"dpr-{numero}-{anno}"

        tipo_short = "dpr"

    elif 'Leggi' in grandparent_folder:

        law_type = "legge"

        law_id = f"legge-{numero}-{anno}"

        tipo_short = "legge"

    else:

        law_type = "altro"

        law_id = f"altro-{numero}-{anno}"

        tipo_short = "altro"

    

    urn = f"urn:nir:stato:{law_type}:{data_iso or anno+'-01-01'};{numero}"

    

    return law_type, law_id, tipo_short, numero, anno, data_iso or f"{anno}-01-01", urn





def extract_legal_metadata(xml_path: Path, extractor: EnhancedLawExtractor) -> Optional[Dict]:

    """Estrae metadata completi da XML"""

    

    try:

        tree = ET.parse(xml_path)

        root = tree.getroot()

        

        # ✅ USA EXTRACTOR

        law_type, law_id, tipo_short, numero, anno, data, urn = determine_law_type_and_id(

            xml_path, root, extractor

        )

        

        # Intestazione

        titolo_doc = root.find('.//nir:titoloDoc', NS)

        titolo = titolo_doc.text if titolo_doc is not None else ""

        

        # Pubblicazione

        pubblicazione = root.find('.//nir:pubblicazione', NS)

        pub_data = pubblicazione.get('norm', '') if pubblicazione is not None else ''

        pub_num = pubblicazione.get('num', '') if pubblicazione is not None else ''

        

        # Converti pub_data in ISO

        pub_data_iso = None

        if len(pub_data) == 8:

            try:

                pub_data_iso = f"{pub_data[:4]}-{pub_data[4:6]}-{pub_data[6:8]}"

            except:

                pass

        

        # Visti

        visti = []

        articolato = root.find('.//nir:articolato', NS)

        if articolato is not None:

            for p in articolato.findall('.//h:p', NS):

                text = extract_text_from_element(p)

                if any(kw in text for kw in ['Visto', 'Vista', 'Visti']):

                    if len(text) > 20:

                        visti.append(text[:500])

        

        # Modifiche

        modifiche = []

        for p in root.findall('.//h:p', NS):

            text = extract_text_from_element(p)

            if any(kw in text.lower() for kw in ['aggiornamento', 'modificato', 'sostituito', 'abrogato']):

                if len(text) > 20:

                    modifiche.append(text[:500])

        

        # Articoli

        articoli = []

        for articolo in root.findall('.//nir:articolo', NS):

            art_id = articolo.get('id', '')

            num_elem = articolo.find('.//nir:num', NS)

            art_num = extract_text_from_element(num_elem) if num_elem is not None else ""

            art_text = extract_text_from_element(articolo)

            

            if art_text and len(art_text) > 50:

                articoli.append({

                    'article_id': art_id,

                    'article_num': art_num,

                    'text': art_text[:2000],

                    'full_text': art_text

                })

        

        return {

            'law_id': law_id,

            'law_type': law_type,

            'tipo_short': tipo_short,

            'numero': numero,

            'anno': anno,

            'data': data,  # ✅ Data ISO reale

            'titolo': titolo,

            'urn': urn,  # ✅ URN reale da XML

            'pubblicazione_data': pub_data_iso or pub_data,

            'pubblicazione_num': pub_num,

            'visti': visti[:10],

            'modifiche': modifiche[:10],

            'articoli': articoli,

            'xml_path': str(xml_path),

            'num_articles': len(articoli)

        }

    

    except Exception as e:

        print(f"❌ Error: {xml_path.name}: {str(e)[:100]}")

        return None





def create_legal_context_payload(metadata: Dict, model: SentenceTransformer):

    """Payload per legal_context"""

    

    visti_text = "\n".join(metadata['visti'][:5])

    modifiche_text = "\n".join(metadata['modifiche'][:3])

    

    embedding_text = f"""

{metadata['law_type']} n. {metadata['numero']} del {metadata['data']}



Titolo: {metadata['titolo']}



URN: {metadata['urn']}

Pubblicazione: GU n.{metadata['pubblicazione_num']} del {metadata['pubblicazione_data']}



Numero articoli: {metadata['num_articles']}



Riferimenti normativi:

{visti_text}



Modifiche:

{modifiche_text}

"""

    

    embedding = model.encode(embedding_text, normalize_embeddings=True)

    

    payload = {

        'law_id': metadata['law_id'],

        'canonical_id': metadata['law_id'],

        'law_type': metadata['law_type'],

        'tipo_short': metadata['tipo_short'],

        'numero': metadata['numero'],

        'anno': metadata['anno'],

        'data': metadata['data'],  # ✅ Data reale

        'titolo': metadata['titolo'][:500],

        'urn': metadata['urn'],  # ✅ URN reale

        'pubblicazione_data': metadata['pubblicazione_data'],

        'pubblicazione_num': metadata['pubblicazione_num'],

        'visti': metadata['visti'],

        'modifiche': metadata['modifiche'],

        'num_articles': metadata['num_articles'],

        'xml_path': metadata['xml_path'],

        'text': embedding_text

    }

    

    return {

        'id': hashlib.md5(metadata['law_id'].encode()).hexdigest(),

        'vector': embedding.tolist(),

        'payload': payload

    }





def create_legal_article_payloads(metadata: Dict, model: SentenceTransformer):

    """Payloads per legal_articles"""

    

    payloads = []

    

    for articolo in metadata['articoli']:

        article_full_id = f"{metadata['law_id']}_{articolo['article_id']}"

        

        embedding_text = f"""

Articolo da {metadata['law_type']} n. {metadata['numero']}/{metadata['anno']}

{articolo['article_num']}



{articolo['text']}



Legge: {metadata['titolo'][:200]}

"""

        

        embedding = model.encode(embedding_text, normalize_embeddings=True)

        

        payload = {

            'article_id': article_full_id,

            'canonical_article_id': article_full_id,

            'parent_law_id': metadata['law_id'],

            'article_num': articolo['article_num'],

            'article_text': articolo['text'],

            'full_text': articolo['full_text'][:5000],

            'law_type': metadata['law_type'],

            'law_numero': metadata['numero'],

            'law_anno': metadata['anno'],

            'law_titolo': metadata['titolo'][:300],

            'text': embedding_text

        }

        

        payloads.append({

            'id': hashlib.md5(article_full_id.encode()).hexdigest(),

            'vector': embedding.tolist(),

            'payload': payload

        })

    

    return payloads





def index_legal_normattiva(client: QdrantClient, model: SentenceTransformer):

    """Index normattiva con EXTRACTOR"""

    

    print("\n" + "="*80)

    print("⚖️ LEGAL NORMATTIVA - FIXED VERSION")

    print("="*80 + "\n")

    

    # ✅ CREATE EXTRACTOR

    extractor = EnhancedLawExtractor()

    

    # Create collections

    print("🗑️ Creating collections...\n")

    

    for coll in ["legal_context", "legal_articles"]:

        try:

            client.delete_collection(coll)

        except:

            pass

        

        client.create_collection(

            collection_name=coll,

            vectors_config=VectorParams(size=1024, distance=Distance.COSINE)

        )

    

    print("✅ Collections created!\n")

    

    # Find XMLs (recursive)

    print("🔍 Scanning XMLs (recursive)...\n")

    

    xml_files = []

    for top_folder in NORMATTIVA_FILTERED.iterdir():

        if top_folder.is_dir() and not top_folder.name.endswith('.json'):

            for law_folder in top_folder.iterdir():

                if law_folder.is_dir():

                    xml_files.extend(list(law_folder.glob("*.xml")))

    

    print(f"📄 Found {len(xml_files):,} XMLs\n")

    

    if not xml_files:

        print("❌ NO XMLs FOUND!\n")

        return 0, 0

    

    # Process

    context_points = []

    article_points = []

    

    for xml_file in tqdm(xml_files, desc="Processing"):

        metadata = extract_legal_metadata(xml_file, extractor)

        

        if metadata and metadata['num_articles'] > 0:

            context_points.append(create_legal_context_payload(metadata, model))

            article_points.extend(create_legal_article_payloads(metadata, model))

    

    print(f"\n✅ Extracted:")

    print(f"   Contexts: {len(context_points):,}")

    print(f"   Articles: {len(article_points):,}\n")

    

    # Upload

    print("⬆️ Uploading legal_context...\n")

    for i in tqdm(range(0, len(context_points), 100), desc="Upload"):

        batch = context_points[i:i+100]

        points = [PointStruct(**p) for p in batch]

        client.upsert(collection_name="legal_context", points=points)

    

    print("\n⬆️ Uploading legal_articles...\n")

    for i in tqdm(range(0, len(article_points), 100), desc="Upload"):

        batch = article_points[i:i+100]

        points = [PointStruct(**p) for p in batch]

        client.upsert(collection_name="legal_articles", points=points)

    

    print(f"\n✅ DONE!")

    print(f"   Contexts: {len(context_points):,}")

    print(f"   Articles: {len(article_points):,}\n")

    

    return len(context_points), len(article_points)





if __name__ == "__main__":

    device = "cuda" if torch.cuda.is_available() else "cpu"

    

    client = QdrantClient(host=config.QDRANT_HOST, port=config.QDRANT_PORT)

    

    print(f"\n🔧 Model: {config.MODEL_NAME}")

    print(f"🖥️ Device: {device.upper()}")

    

    model = SentenceTransformer(config.MODEL_NAME, device=device)

    print(f"✅ Loaded\n")

    

    n_contexts, n_articles = index_legal_normattiva(client, model)

    

    print(f"\n{'='*80}")

    print(f"✅ COMPLETE!")

    print(f"   Contexts: {n_contexts:,}")

    print(f"   Articles: {n_articles:,}")

    print(f"{'='*80}\n")
