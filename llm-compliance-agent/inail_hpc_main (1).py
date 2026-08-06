#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INAIL Web Scraper - HPC Version
Adattato per UB HPC con ColdFront
"""

import argparse
import sys
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional, List
import unicodedata
import re
import time
import json
import random

# SELENIUM
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from webdriver_manager.chrome import ChromeDriverManager

# WEB SCRAPING
from bs4 import BeautifulSoup
from urllib.parse import urljoin, quote_plus

# Configurazione cache per HPC (evita riempire home directory)
if 'HF_HOME' not in os.environ:
    hf_cache = Path.home() / '.cache' / 'huggingface'
    os.environ['HF_HOME'] = str(hf_cache)
    os.environ['TRANSFORMERS_CACHE'] = str(hf_cache)

BASE_CATALOG_URL = "https://www.inail.it/portale/it/inail-comunica/pubblicazioni/catalogo-generale.html"

# ============================================================================
# CONFIGURAZIONE
# ============================================================================

class INAILConfig:
    """Configurazione per ambiente HPC"""
    
    # Path verranno settati da argomenti command-line
    OUTPUT_DIR = Path('/tmp/inail_data')  # Default temporaneo
    PDF_DIR = OUTPUT_DIR / 'pdfs'
    JSON_DIR = OUTPUT_DIR / 'json'
    IMAGES_DIR = OUTPUT_DIR / 'images'
    TABLES_DIR = OUTPUT_DIR / 'tables'

    # VLM Settings
    VLM_ENABLED = True
    VLM_MODEL = "HuggingFaceTB/SmolVLM-256M-Instruct"
    VLM_USE_BFLOAT16 = True
    VLM_USE_FLASH_ATTENTION = True
    VLM_USE_8BIT_QUANTIZATION = False
    VLM_IMAGE_SIZE = 512
    VLM_MAX_NEW_TOKENS = 300
    VLM_BATCH_SIZE = 8

    VLM_PROMPT = """Describe this technical document image in detail for academic research. Focus on:
- Main content type (diagram, chart, photo, schematic)
- Key visual elements and their relationships
- Any visible text, labels, or numerical data
- Technical details relevant to workplace safety or industrial processes
Be precise and comprehensive."""

    @classmethod
    def setup_directories(cls):
        """Crea tutte le directory necessarie"""
        for dir_path in [cls.OUTPUT_DIR, cls.PDF_DIR, cls.JSON_DIR,
                        cls.IMAGES_DIR, cls.TABLES_DIR]:
            dir_path.mkdir(exist_ok=True, parents=True)
        print(f"[INFO] Directory configurate: {cls.OUTPUT_DIR}")

    @classmethod
    def verify_directories(cls):
        """Verifica che tutte le directory esistano"""
        all_exist = True
        for name, path in [
            ('OUTPUT', cls.OUTPUT_DIR),
            ('PDFs', cls.PDF_DIR),
            ('JSON', cls.JSON_DIR),
            ('Images', cls.IMAGES_DIR),
            ('Tables', cls.TABLES_DIR)
        ]:
            exists = path.exists()
            status = "✅" if exists else "❌"
            print(f"  {status} {name}: {path}")
            all_exist = all_exist and exists
        return all_exist

# ============================================================================
# SELENIUM DRIVER (versione HPC)
# ============================================================================

def create_driver():
    """Crea driver Selenium per HPC"""
    options = Options()
    options.add_argument('--headless')
    options.add_argument('--no-sandbox')
    options.add_argument('--disable-dev-shm-usage')
    options.add_argument('--disable-gpu')
    options.add_argument('--disable-extensions')
    options.add_argument('user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36')
    
    try:
        # Usa ChromeDriverManager per gestire automaticamente il driver
        service = Service(ChromeDriverManager().install())
        driver = webdriver.Chrome(service=service, options=options)
    except Exception as e:
        print(f"[WARNING] ChromeDriverManager fallito: {e}")
        print("[INFO] Tentativo con chromedriver di sistema...")
        driver = webdriver.Chrome(options=options)
    
    driver.set_page_load_timeout(180)
    driver.implicitly_wait(20)
    
    print("[INFO] Driver Selenium inizializzato")
    return driver

# ============================================================================
# DOCLING
# ============================================================================

def initialize_docling_converter():
    """Inizializza Docling per parsing PDF"""
    from docling.datamodel.pipeline_options import PdfPipelineOptions
    from docling.datamodel.base_models import InputFormat
    from docling.document_converter import DocumentConverter, PdfFormatOption

    pipeline_options = PdfPipelineOptions()
    pipeline_options.generate_picture_images = True
    pipeline_options.images_scale = 2.0
    pipeline_options.do_ocr = True
    pipeline_options.do_picture_description = False

    converter = DocumentConverter(
        format_options={
            InputFormat.PDF: PdfFormatOption(pipeline_options=pipeline_options)
        }
    )

    print(f"[INFO] Docling OK | VLM: {'ON' if INAILConfig.VLM_ENABLED else 'OFF'}")
    return converter

# ============================================================================
# URL UTILITIES
# ============================================================================

def normalize_query(q: str) -> str:
    """Normalizza query di ricerca"""
    q = unicodedata.normalize("NFKC", q).strip().lower()
    return " ".join(q.split()).replace("'", "'")

def validate_date(date_str: str) -> bool:
    """Valida formato data gg/mm/aaaa"""
    return bool(re.match(r"^\d{2}/\d{2}/\d{4}$", date_str))

def build_inail_search_url(query: str = None, page: int = 1,
                          start_date: str = None, end_date: str = None) -> str:
    """Costruisce URL di ricerca INAIL"""
    if page < 1:
        raise ValueError("Pagina >= 1")

    params = []
    if query and query.strip():
        params.append(f"text={quote_plus(normalize_query(query))}")
    if start_date and validate_date(start_date):
        params.append(f"startDate={quote_plus(start_date)}")
    if end_date and validate_date(end_date):
        params.append(f"endDate={quote_plus(end_date)}")
    params.append(f"page={page}")

    return f"{BASE_CATALOG_URL}?{'&'.join(params)}"

# ============================================================================
# ESTRAZIONE LINK
# ============================================================================

def extract_cards_with_wait(driver, timeout: int = 20, max_retries: int = 4) -> List[Dict[str, str]]:
    """Estrae card pubblicazioni con retry"""
    for attempt in range(max_retries):
        try:
            print(f"    [Tentativo {attempt + 1}/{max_retries}]", end=" ")

            WebDriverWait(driver, timeout).until(
                EC.presence_of_element_located((By.CSS_SELECTOR, "h3.card-title a[href], body"))
            )
            time.sleep(3)

            results = []
            cards_elements = driver.find_elements(By.CSS_SELECTOR, "h3.card-title a[href]")

            for a in cards_elements:
                try:
                    href = a.get_attribute("href")
                    text = a.text.strip()
                    if href and text:
                        results.append({
                            "titolo": text,
                            "url": urljoin("https://www.inail.it", href)
                        })
                except:
                    continue

            if results:
                print(f"✓ {len(results)} card")
                return results
            else:
                print("✗ Nessuna card")

        except TimeoutException:
            print(f"✗ Timeout")
            if attempt < max_retries - 1:
                wait_time = random.uniform(10, 15)
                print(f"    [Pausa {wait_time:.1f}s]")
                time.sleep(wait_time)
                try:
                    driver.refresh()
                    time.sleep(5)
                except:
                    pass

        except Exception as e:
            print(f"✗ Errore: {str(e)[:50]}")
            if attempt < max_retries - 1:
                time.sleep(8)

    print("    [FALLITO]")
    return []

def get_inail_publication_links_with_range(
    driver,
    query: str = None,
    start_page: int = 1,
    end_page: int = 5,
    start_date: str = None,
    end_date: str = None
) -> List[Dict[str, str]]:
    """Estrae link pubblicazioni con range di pagine personalizzato"""
    
    all_links = []
    seen = set()
    consecutive_failures = 0
    MAX_FAILURES = 2

    print(f"\n{'='*80}")
    print(f"ESTRAZIONE LINK CON RANGE PAGINE")
    print(f"{'='*80}")
    print(f"Query: {query if query else '(tutte)'}")
    print(f"Range: pagina {start_page} → {end_page}")
    print(f"{'='*80}\n")

    for page_num in range(start_page, end_page + 1):
        url = build_inail_search_url(
            query=query,
            page=page_num,
            start_date=start_date,
            end_date=end_date
        )

        print(f"[Pagina {page_num}/{end_page}]")

        try:
            driver.get(url)

            if page_num > start_page:
                delay = random.uniform(6, 10)
                print(f"  [Pausa {delay:.1f}s]")
                time.sleep(delay)

            cards = extract_cards_with_wait(driver, timeout=20)

            if not cards:
                consecutive_failures += 1
                print(f"  Pagina vuota ({consecutive_failures}/{MAX_FAILURES})")

                if consecutive_failures >= MAX_FAILURES:
                    print(f"\n  STOP: Troppe pagine vuote")
                    break

                time.sleep(random.uniform(10, 15))
                continue

            consecutive_failures = 0
            new_count = 0

            for card in cards:
                if card["url"] not in seen:
                    seen.add(card["url"])
                    all_links.append(card)
                    new_count += 1

            print(f"  +{new_count} nuove (totale: {len(all_links)})")
            time.sleep(random.uniform(2, 4))

        except KeyboardInterrupt:
            print(f"\n⚠️ INTERRUZIONE MANUALE")
            break

        except Exception as e:
            consecutive_failures += 1
            error_msg = str(e)
            print(f"  Errore: {error_msg[:100]}")

            if consecutive_failures >= MAX_FAILURES:
                print(f"\n  STOP: Troppi errori consecutivi")
                break

            time.sleep(random.uniform(12, 18))

    print(f"\n{'='*80}")
    print(f"Totale: {len(all_links)} pubblicazioni trovate")
    print(f"{'='*80}\n")

    return all_links

# ============================================================================
# PDF PROCESSING
# ============================================================================

def download_pdf_inail(pdf_url: str, output_path: Path) -> bool:
    """Scarica PDF da URL"""
    try:
        os.system(f'wget -q -O {output_path} "{pdf_url}"')
        return output_path.exists() and output_path.stat().st_size > 0
    except:
        return False

def export_images_from_document(doc, doc_id: str) -> List[Dict[str, Any]]:
    """Esporta immagini da documento Docling"""
    images_dir = INAILConfig.IMAGES_DIR / doc_id
    images_dir.mkdir(parents=True, exist_ok=True)
    exported_images = []

    try:
        from docling_core.types.doc import PictureItem
        if hasattr(doc, 'iterate_items'):
            for element, _level in doc.iterate_items():
                if isinstance(element, PictureItem):
                    picture_counter = len(exported_images) + 1
                    img_filename = f"figure_{picture_counter}.png"
                    img_path = images_dir / img_filename

                    try:
                        pil_image = element.get_image(doc)
                        pil_image.save(str(img_path), "PNG")
                        exported_images.append({
                            'path': str(img_path),
                            'filename': img_filename,
                            'caption': getattr(element, 'caption', 'N/A'),
                            'position': len(exported_images),
                            'vlm_description': None
                        })
                        print(f"  [+] {img_filename}")
                    except Exception as e:
                        print(f"  [!] Errore figura: {e}")

        print(f"[OK] {len(exported_images)} immagini esportate")
    except Exception as e:
        print(f"[ERROR] Export immagini: {e}")

    return exported_images

def export_tables_to_files(tables: List[Dict], doc_id: str) -> List[str]:
    """Esporta tabelle come file .txt"""
    if not tables:
        return []

    tables_dir = INAILConfig.TABLES_DIR / doc_id
    tables_dir.mkdir(parents=True, exist_ok=True)
    exported_tables = []

    for table in tables:
        table_id = table['table_id']
        txt_path = tables_dir / f"{table_id}.txt"

        try:
            with open(txt_path, 'w', encoding='utf-8') as f:
                f.write("="*80 + "\n")
                f.write(f"TABELLA: {table_id}\n")
                f.write("="*80 + "\n\n")
                f.write(f"CAPTION: {table.get('caption', 'N/A')}\n\n")

                if table.get('potential_columns'):
                    f.write("COLONNE IDENTIFICATE:\n")
                    for idx, col in enumerate(table['potential_columns'], 1):
                        f.write(f"  {idx}. {col}\n")
                    f.write("\n")

                f.write("METADATI STRUTTURALI:\n")
                f.write(f"  - Righe: {table.get('num_rows', 'N/A')}\n")
                f.write(f"  - Posizione: {table.get('position', 'N/A')}\n")
                f.write(f"  - Colonne: {len(table.get('potential_columns', []))}\n\n")

                f.write("CONTENUTO:\n" + "-"*80 + "\n")
                f.write(table['text_content'])
                f.write("\n" + "-"*80 + "\n")

            exported_tables.append(str(txt_path))
            print(f"  [+] {txt_path.name}")
        except Exception as e:
            print(f"  [!] Errore tabella: {e}")

    return exported_tables

def analyze_document_structure(doc) -> Dict[str, Any]:
    """Analizza struttura documento"""
    structure = {
        'num_tables': 0,
        'num_figures': 0,
        'num_headings': 0,
        'headings': [],
        'tables': [],
        'figures': []
    }

    try:
        markdown = doc.export_to_markdown()
        lines = markdown.split('\n')

        # Headings
        for i, line in enumerate(lines):
            if line.strip().startswith('#'):
                level = len(line) - len(line.lstrip('#'))
                text = line.lstrip('#').strip()
                structure['num_headings'] += 1
                structure['headings'].append({
                    'type': f'heading_level_{level}',
                    'text': text[:200],
                    'position': i,
                    'level': level
                })

        # Tables
        in_table = False
        current_table_lines = []
        table_start_idx = -1

        for i, line in enumerate(lines):
            if '|' in line and ('-' in line or '─' in line):
                if not in_table:
                    in_table = True
                    table_start_idx = i
                    current_table_lines = []

            if in_table:
                if '|' in line:
                    current_table_lines.append(line)
                else:
                    if current_table_lines:
                        caption = 'N/A'
                        for j in range(max(0, table_start_idx-5), table_start_idx):
                            candidate = lines[j].strip()
                            if candidate and 'Tabella' in candidate:
                                caption = candidate
                                break

                        potential_columns = []
                        if current_table_lines:
                            first_row = current_table_lines[0]
                            cols = [c.strip() for c in first_row.split('|') if c.strip()]
                            potential_columns = [col for col in cols if col and col not in ['-', '─']]

                        structure['num_tables'] += 1
                        structure['tables'].append({
                            'table_id': f'table_{structure["num_tables"]}',
                            'caption': caption,
                            'text_content': '\n'.join(current_table_lines),
                            'num_rows': len(current_table_lines),
                            'num_columns': len(potential_columns),
                            'potential_columns': potential_columns,
                            'position': table_start_idx
                        })

                    in_table = False
                    current_table_lines = []

        # Figures
        for i, line in enumerate(lines):
            if '<!-- image -->' in line.lower() or (line.strip().startswith('![') and '](' in line):
                structure['num_figures'] += 1
                caption = 'N/A'
                for j in range(i+1, min(len(lines), i+5)):
                    candidate = lines[j].strip()
                    if candidate and not candidate.startswith('#'):
                        caption = candidate
                        break

                structure['figures'].append({
                    'figure_id': f'figure_{structure["num_figures"]}',
                    'caption': caption,
                    'position': i,
                    'vlm_description': None
                })

    except Exception as e:
        print(f"[!] Errore analisi: {e}")

    return structure

def add_vlm_descriptions_to_images(
    exported_images: List[Dict],
    doc_id: str
) -> List[Dict]:
    """Versione ottimizzata con batch processing e ottimizzazioni HPC"""
    
    if not exported_images or not INAILConfig.VLM_ENABLED:
        return exported_images

    try:
        from transformers import AutoProcessor, AutoModelForVision2Seq, BitsAndBytesConfig
        from PIL import Image
        import torch

        print(f"\n[VLM] Caricamento ottimizzato di {INAILConfig.VLM_MODEL}...")

        device = "cuda" if torch.cuda.is_available() else "cpu"

        # Inizializza processor
        processor = AutoProcessor.from_pretrained(
            INAILConfig.VLM_MODEL,
            trust_remote_code=True,
            size={"longest_edge": INAILConfig.VLM_IMAGE_SIZE}
        )

        # Configura modello
        model_kwargs = {
            "trust_remote_code": True,
            "device_map": "auto"
        }

        if INAILConfig.VLM_USE_8BIT_QUANTIZATION:
            print("[VLM] Usando quantizzazione 8-bit")
            quantization_config = BitsAndBytesConfig(load_in_8bit=True)
            model_kwargs["quantization_config"] = quantization_config
        else:
            if INAILConfig.VLM_USE_BFLOAT16 and torch.cuda.is_available():
                print("[VLM] Usando bfloat16 precision")
                model_kwargs["torch_dtype"] = torch.bfloat16
            else:
                model_kwargs["torch_dtype"] = torch.float16 if torch.cuda.is_available() else torch.float32

        if INAILConfig.VLM_USE_FLASH_ATTENTION and device == "cuda":
            print("[VLM] Usando Flash Attention 2")
            model_kwargs["_attn_implementation"] = "flash_attention_2"
        else:
            model_kwargs["_attn_implementation"] = "eager"

        model = AutoModelForVision2Seq.from_pretrained(
            INAILConfig.VLM_MODEL,
            **model_kwargs
        )

        print(f"[VLM] Modello caricato su {device}")
        print(f"[VLM] Risoluzione: {INAILConfig.VLM_IMAGE_SIZE}x{INAILConfig.VLM_IMAGE_SIZE}")
        print(f"[VLM] Batch size: {INAILConfig.VLM_BATCH_SIZE}")

        # Batch processing
        total_images = len(exported_images)
        batch_size = INAILConfig.VLM_BATCH_SIZE

        for batch_start in range(0, total_images, batch_size):
            batch_end = min(batch_start + batch_size, total_images)
            batch = exported_images[batch_start:batch_end]

            print(f"\n[VLM] Batch {batch_start//batch_size + 1}/{(total_images + batch_size - 1)//batch_size}")

            images = []
            valid_indices = []

            for idx, img_data in enumerate(batch):
                try:
                    image = Image.open(img_data['path']).convert('RGB')

                    max_size = INAILConfig.VLM_IMAGE_SIZE * 4
                    if max(image.size) > max_size:
                        ratio = max_size / max(image.size)
                        new_size = tuple(int(dim * ratio) for dim in image.size)
                        image = image.resize(new_size, Image.Resampling.LANCZOS)

                    images.append(image)
                    valid_indices.append(batch_start + idx)

                except Exception as e:
                    print(f"  [{batch_start + idx + 1}] Errore caricamento: {e}")
                    exported_images[batch_start + idx]['vlm_description'] = None
                    continue

            if not images:
                continue

            try:
                messages_batch = []
                for _ in images:
                    messages_batch.append([{
                        "role": "user",
                        "content": [
                            {"type": "image"},
                            {"type": "text", "text": INAILConfig.VLM_PROMPT}
                        ]
                    }])

                prompts = [
                    processor.apply_chat_template(msgs, add_generation_prompt=True)
                    for msgs in messages_batch
                ]

                inputs = processor(
                    text=prompts,
                    images=images,
                    return_tensors="pt",
                    padding=True
                )
                inputs = {k: v.to(device) for k, v in inputs.items()}

                with torch.no_grad():
                    outputs = model.generate(
                        **inputs,
                        max_new_tokens=INAILConfig.VLM_MAX_NEW_TOKENS,
                        do_sample=False
                    )

                generated_texts = processor.batch_decode(
                    outputs,
                    skip_special_tokens=True
                )

                for i, (global_idx, text) in enumerate(zip(valid_indices, generated_texts)):
                    if "Assistant:" in text:
                        description = text.split("Assistant:")[-1].strip()
                    else:
                        description = text.strip()

                    exported_images[global_idx]['vlm_description'] = description

                    img_path = Path(exported_images[global_idx]['path'])
                    txt_path = img_path.with_suffix('.txt')

                    with open(txt_path, 'w', encoding='utf-8') as f:
                        f.write("="*80 + "\n")
                        f.write(f"VLM DESCRIPTION - {exported_images[global_idx]['filename']}\n")
                        f.write("="*80 + "\n\n")
                        f.write("IMAGE METADATA:\n")
                        f.write(f"  - Filename: {exported_images[global_idx]['filename']}\n")
                        f.write(f"  - Caption: {exported_images[global_idx].get('caption', 'N/A')}\n")
                        f.write(f"  - Position: {exported_images[global_idx].get('position', 'N/A')}\n")
                        f.write(f"  - Model: {INAILConfig.VLM_MODEL}\n")
                        f.write(f"  - Timestamp: {datetime.now().isoformat()}\n\n")
                        f.write("VLM ANALYSIS:\n")
                        f.write("-"*80 + "\n")
                        f.write(description)
                        f.write("\n" + "-"*80 + "\n")

                    exported_images[global_idx]['vlm_description_file'] = str(txt_path)
                    print(f"  [{global_idx + 1}/{total_images}] ✓ {description[:50]}...")

            except Exception as e:
                print(f"  Errore batch processing: {e}")
                for idx in valid_indices:
                    exported_images[idx]['vlm_description'] = None

        # Cleanup
        del model, processor
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        success = sum(1 for img in exported_images if img.get('vlm_description'))
        print(f"\n[VLM] Completato: {success}/{total_images} immagini")

    except Exception as e:
        print(f"[ERROR] VLM fallito: {e}")
        import traceback
        traceback.print_exc()

    return exported_images

# ============================================================================
# SCRAPING
# ============================================================================

def scrape_inail_publication(driver, publication_url: str, converter,
                             save_pdf: bool = True, enable_vlm: Optional[bool] = None) -> Dict[str, Any]:
    """Scrape singola pubblicazione INAIL"""
    use_vlm = enable_vlm if enable_vlm is not None else INAILConfig.VLM_ENABLED

    print(f"\n{'='*80}")
    print(f"SCRAPING PUBBLICAZIONE")
    print(f"{'='*80}\n")

    scraping_metadata = {
        'url': publication_url,
        'timestamp': datetime.now().isoformat(),
        'scraper_version': '6.0_HPC',
        'docling_used': True,
        'vlm_enabled': use_vlm
    }

    try:
        driver.get(publication_url)
        time.sleep(3)

        page = BeautifulSoup(driver.page_source, "lxml")

        # Metadati
        title_elem = page.find("h2", class_="h1")
        title = title_elem.get_text(strip=True) if title_elem else "N/A"

        descr_blocks = page.find_all("p", class_="text-20")
        abstract = descr_blocks[1].get_text(strip=True) if len(descr_blocks) > 1 else ""

        data_elem = page.find("strong", class_="js-date-value")
        data_pub = data_elem.get_text(strip=True).split(", ")[0] if data_elem else "N/A"

        link_pdf = page.find("ul", class_="list-download")
        pdf_url = None
        if link_pdf:
            a_tag = link_pdf.find("a", href=True)
            if a_tag:
                pdf_url = urljoin("https://www.inail.it", a_tag["href"])

        print(f"[+] {title[:60]}...")
        print(f"[+] Data: {data_pub}")

        # PDF Processing
        pdf_data = None
        pdf_local_path = None

        if pdf_url:
            safe_filename = re.sub(r'[^\w\-_.]', '_', title[:50]) + '.pdf'
            pdf_local_path = INAILConfig.PDF_DIR / safe_filename
            doc_id = safe_filename.replace('.pdf', '')

            if download_pdf_inail(pdf_url, pdf_local_path):
                print(f"[INFO] Processing PDF...")

                result = converter.convert(str(pdf_local_path))
                doc = result.document

                markdown_text = doc.export_to_markdown()
                plain_text = doc.export_to_text()
                structure_info = analyze_document_structure(doc)
                exported_images = export_images_from_document(doc, doc_id)

                if use_vlm and exported_images:
                    exported_images = add_vlm_descriptions_to_images(exported_images, doc_id)
                    for idx, fig_info in enumerate(structure_info.get('figures', [])):
                        if idx < len(exported_images):
                            vlm_desc = exported_images[idx].get('vlm_description')
                            if vlm_desc:
                                fig_info['vlm_description'] = vlm_desc

                exported_tables = export_tables_to_files(structure_info['tables'], doc_id)

                pdf_data = {
                    'markdown_content': markdown_text,
                    'plain_text': plain_text,
                    'num_pages': len(doc.pages) if hasattr(doc, 'pages') else 'N/A',
                    'num_tables': structure_info['num_tables'],
                    'num_figures': structure_info['num_figures'],
                    'num_headings': structure_info['num_headings'],
                    'headings': structure_info['headings'],
                    'tables': structure_info['tables'],
                    'figures': structure_info['figures'],
                    'exported_images': exported_images,
                    'exported_tables': exported_tables
                }

                print(f"[OK] {pdf_data['num_pages']} pag | {pdf_data['num_tables']} tab | {pdf_data['num_figures']} fig")

                if not save_pdf:
                    pdf_local_path.unlink()
                    pdf_local_path = None

        return {
            'scraping_metadata': scraping_metadata,
            'web_metadata': {
                'title': title,
                'abstract': abstract,
                'data_pubblicazione': data_pub,
                'pdf_url': pdf_url,
                'pdf_local_path': str(pdf_local_path) if pdf_local_path else None
            },
            'document_content': pdf_data,
            'status': 'success',
            'has_pdf': pdf_url is not None,
            'pdf_processed': pdf_data is not None
        }

    except Exception as e:
        print(f"[ERROR] {e}")
        return {
            'scraping_metadata': scraping_metadata,
            'status': 'error',
            'error_message': str(e)
        }

def save_to_json(data: Dict[str, Any], output_name: Optional[str] = None) -> Path:
    """Salva dati in JSON con verifica"""
    try:
        if output_name is None:
            title = data.get('web_metadata', {}).get('title', 'unknown')
            safe_title = re.sub(r'[^\w\-_.]', '_', title[:50])
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_name = f"{safe_title}_{timestamp}.json"

        output_path = INAILConfig.JSON_DIR / output_name

        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)

        if output_path.exists():
            size_kb = output_path.stat().st_size / 1024
            print(f"  JSON: {output_name} ({size_kb:.1f} KB)")
            return output_path
        else:
            print(f"  ERRORE: File non creato!")
            return None

    except Exception as e:
        print(f"  ERRORE SALVATAGGIO: {e}")
        return None

def load_already_scraped_urls() -> set:
    """Carica URL già processati dai JSON"""
    scraped_urls = set()

    json_dir = INAILConfig.JSON_DIR
    if not json_dir.exists():
        print(f"  Directory JSON non esiste ancora")
        return scraped_urls

    json_files = list(json_dir.glob("*.json"))
    json_files = [f for f in json_files if not f.name.startswith('vector_db_manifest')]

    for json_file in json_files:
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)

            docs = data if isinstance(data, list) else [data]

            for doc in docs:
                url = doc.get('scraping_metadata', {}).get('url')
                if url:
                    scraped_urls.add(url)
        except Exception as e:
            print(f"  Errore lettura {json_file.name}: {e}")
            continue

    return scraped_urls

def batch_scrape_inail_publications_incremental(
    driver,
    converter,
    query: str = None,
    start_page: int = 1,
    end_page: int = 5,
    max_documents: Optional[int] = None,
    enable_vlm: bool = False,
    save_pdfs: bool = False,
    delay_between_docs: int = 4
) -> List[Dict[str, Any]]:
    """Batch scraping incrementale con range pagine"""
    
    # Carica URL già processati
    already_scraped = load_already_scraped_urls()

    print(f"\n{'='*80}")
    print(f"SCRAPING INCREMENTALE - RANGE PAGINE")
    print(f"{'='*80}")
    print(f"Documenti già nel database: {len(already_scraped)}")
    print(f"Pagine da analizzare: {start_page} → {end_page}")
    print(f"{'='*80}\n")

    # Estrai link con range custom
    publication_links = get_inail_publication_links_with_range(
        driver,
        query=query,
        start_page=start_page,
        end_page=end_page
    )

    if not publication_links:
        print("[WARNING] Nessuna pubblicazione trovata in questo range")
        return []

    # Filtra SOLO quelli nuovi
    new_links = [pub for pub in publication_links if pub['url'] not in already_scraped]

    print(f"\n{'='*80}")
    print(f"FILTRO DOCUMENTI")
    print(f"{'='*80}")
    print(f"Trovati nelle pagine {start_page}-{end_page}: {len(publication_links)}")
    print(f"NUOVI (da processare): {len(new_links)}")
    print(f"GIÀ PROCESSATI (skip): {len(publication_links) - len(new_links)}")
    print(f"{'='*80}\n")

    if not new_links:
        print("Nessun nuovo documento in questo range di pagine!")
        return []

    print("📋 NUOVE PUBBLICAZIONI DA SCARICARE:")
    for i, pub in enumerate(new_links[:10], 1):
        print(f"  {i}. {pub['titolo'][:70]}...")
    if len(new_links) > 10:
        print(f"  ... e altre {len(new_links) - 10}\n")
    else:
        print()

    if max_documents:
        new_links = new_links[:max_documents]
        print(f"[INFO] Limitato a {max_documents} nuovi documenti\n")

    # Scraping
    results = []
    failed = []

    print(f"{'='*80}")
    print(f"SCRAPING {len(new_links)} NUOVI DOCUMENTI")
    print(f"{'='*80}\n")

    for idx, pub in enumerate(new_links, 1):
        print(f"\n[{idx}/{len(new_links)}] {pub['titolo'][:60]}...")

        try:
            result = scrape_inail_publication(
                driver=driver,
                publication_url=pub['url'],
                converter=converter,
                save_pdf=save_pdfs,
                enable_vlm=enable_vlm
            )

            if result['status'] == 'success':
                json_path = save_to_json(result)

                if json_path and json_path.exists():
                    results.append(result)
                    print(f"  ✓ SALVATO")
                else:
                    print(f"  ✗ Processato ma NON salvato")
                    failed.append({
                        'url': pub['url'],
                        'title': pub['titolo'],
                        'error': 'JSON save failed'
                    })
            else:
                failed.append({
                    'url': pub['url'],
                    'title': pub['titolo'],
                    'error': result.get('error_message', 'Unknown')
                })
                print(f"  ✗ FAIL")

        except KeyboardInterrupt:
            print(f"\n⚠️ INTERRUZIONE - Salvati: {len(results)}")
            break

        except Exception as e:
            print(f"  ✗ Eccezione: {str(e)[:80]}")
            failed.append({'url': pub['url'], 'title': pub['titolo'], 'error': str(e)})

        if idx < len(new_links):
            delay = random.uniform(delay_between_docs, delay_between_docs + 2)
            time.sleep(delay)

    # Riepilogo
    print(f"\n{'='*80}")
    print(f"COMPLETATO")
    print(f"{'='*80}")
    print(f"Nuovi salvati: {len(results)}")
    print(f"Falliti: {len(failed)}")
    print(f"Totale database: {len(already_scraped) + len(results)}")
    print(f"Path: {INAILConfig.OUTPUT_DIR}")
    print(f"{'='*80}\n")

    if failed:
        print("\n❌ DOCUMENTI FALLITI:")
        for f in failed[:5]:
            print(f"  • {f['title'][:50]}: {f['error']}")
        if len(failed) > 5:
            print(f"  ... e altri {len(failed)-5}")

    return results

# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def parse_arguments():
    """Parser argomenti command-line per batch script"""
    parser = argparse.ArgumentParser(
        description='INAIL Web Scraper per HPC',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Esempi d'uso:
  
  # Scraping pagine 1-5 senza VLM
  python inail_scraper_hpc.py --start-page 1 --end-page 5 --no-vlm
  
  # Scraping con query specifica
  python inail_scraper_hpc.py --query "sicurezza lavoro" --start-page 1 --end-page 10
  
  # Scraping con VLM e salvataggio PDF
  python inail_scraper_hpc.py --start-page 11 --end-page 20 --enable-vlm --save-pdfs
        """
    )

    parser.add_argument(
        '--query',
        type=str,
        default=None,
        help='Query di ricerca (default: tutte le pubblicazioni)'
    )

    parser.add_argument(
        '--start-page',
        type=int,
        default=1,
        help='Pagina di inizio (default: 1)'
    )

    parser.add_argument(
        '--end-page',
        type=int,
        default=5,
        help='Pagina di fine (default: 5)'
    )

    parser.add_argument(
        '--max-docs',
        type=int,
        default=None,
        help='Massimo numero di documenti da processare per questo run'
    )

    parser.add_argument(
        '--enable-vlm',
        action='store_true',
        default=False,
        help='Abilita descrizioni VLM per immagini (richiede GPU)'
    )

    parser.add_argument(
        '--no-vlm',
        action='store_true',
        default=False,
        help='Disabilita esplicitamente VLM'
    )

    parser.add_argument(
        '--save-pdfs',
        action='store_true',
        default=False,
        help='Salva file PDF (usa molto spazio)'
    )

    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Directory di output (default: /projects/academic/[group]/inail_thesis/data)'
    )

    parser.add_argument(
        '--delay',
        type=int,
        default=4,
        help='Delay in secondi tra documenti (default: 4)'
    )

    return parser.parse_args()

# ============================================================================
# MAIN
# ============================================================================

def main_hpc():
    """Main function per esecuzione HPC"""
    args = parse_arguments()

    print(f"\n{'='*80}")
    print("INAIL SCRAPER - HPC MODE")
    print(f"{'='*80}\n")

    # Determina output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        # Prova a rilevare automaticamente il gruppo
        import pwd
        import grp
        user = pwd.getpwuid(os.getuid()).pw_name
        groups = [grp.getgrgid(g).gr_name for g in os.getgroups()]
        
        # Cerca gruppo academic
        academic_groups = [g for g in groups if 'ccr' in g.lower() or 'academic' in g.lower()]
        
        if academic_groups:
            group = academic_groups[0]
            output_dir = Path(f"/projects/academic/{group}/inail_thesis/data")
            print(f"[INFO] Gruppo rilevato automaticamente: {group}")
        else:
            print("[ERROR] Impossibile rilevare gruppo ColdFront")
            print("Specifica --output-dir manualmente")
            sys.exit(1)

    # Configura path
    INAILConfig.OUTPUT_DIR = output_dir
    INAILConfig.PDF_DIR = output_dir / 'pdfs'
    INAILConfig.JSON_DIR = output_dir / 'json'
    INAILConfig.IMAGES_DIR = output_dir / 'images'
    INAILConfig.TABLES_DIR = output_dir / 'tables'

    # Configura VLM
    if args.no_vlm:
        enable_vlm = False
    else:
        enable_vlm = args.enable_vlm

    INAILConfig.VLM_ENABLED = enable_vlm

    print(f"Output directory: {INAILConfig.OUTPUT_DIR}")
    print(f"Query: {args.query if args.query else '(tutte)'}")
    print(f"Pagine: {args.start_page} → {args.end_page}")
    print(f"VLM: {'ON' if enable_vlm else 'OFF'}")
    print(f"Save PDFs: {'YES' if args.save_pdfs else 'NO'}")
    print(f"Delay: {args.delay}s tra documenti")
    print(f"{'='*80}\n")

    # Setup directories
    try:
        INAILConfig.setup_directories()
        
        if not INAILConfig.verify_directories():
            print("[ERROR] Impossibile creare directory")
            sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Setup directory fallito: {e}")
        sys.exit(1)

    # Inizializza driver e converter
    print("\n[SETUP] Inizializzazione componenti...")
    try:
        driver = create_driver()
        converter = initialize_docling_converter()
    except Exception as e:
        print(f"[ERROR] Inizializzazione fallita: {e}")
        sys.exit(1)

    # Carica documenti già processati
    already_scraped = load_already_scraped_urls()
    print(f"\n[INFO] Documenti già nel database: {len(already_scraped)}\n")

    # Esegui scraping
    try:
        results = batch_scrape_inail_publications_incremental(
            driver=driver,
            converter=converter,
            query=args.query,
            start_page=args.start_page,
            end_page=args.end_page,
            max_documents=args.max_docs,
            enable_vlm=enable_vlm,
            save_pdfs=args.save_pdfs,
            delay_between_docs=args.delay
        )
    except Exception as e:
        print(f"\n[ERROR] Scraping fallito: {e}")
        import traceback
        traceback.print_exc()
        driver.quit()
        sys.exit(1)

    # Cleanup
    driver.quit()

    # Report finale
    print(f"\n{'='*80}")
    print("SCRAPING COMPLETATO")
    print(f"{'='*80}")
    print(f"Nuovi documenti: {len(results)}")
    print(f"Totale database: {len(load_already_scraped_urls())}")
    print(f"Output: {INAILConfig.OUTPUT_DIR}")
    print(f"{'='*80}\n")

    return len(results)

if __name__ == "__main__":
    try:
        num_processed = main_hpc()
        sys.exit(0 if num_processed >= 0 else 1)
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Scraping interrotto dall'utente")
        sys.exit(130)
    except Exception as e:
        print(f"\n[FATAL ERROR] {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)