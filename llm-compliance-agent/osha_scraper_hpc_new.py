# -*- coding: utf-8 -*-
"""
OSHA Publications Scraper - HPC Version
Adattato per ambiente HPC senza Google Drive
"""

import re
import requests
from bs4 import BeautifulSoup
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, Any, List
import tempfile
import os
import base64
from io import BytesIO
from urllib.parse import quote_plus, urljoin
import unicodedata
import time
import random

# =========================================================
# IMPORT DOCLING
# =========================================================
from docling.document_converter import DocumentConverter, PdfFormatOption
from docling.datamodel.base_models import InputFormat
from docling.datamodel.pipeline_options import PdfPipelineOptions
from docling_core.types.doc import ImageRefMode

# =========================================================
# ENHANCED TABLE PARSER
# =========================================================

def extract_table_metadata_enhanced(table_lines: List[str], context_lines: List[str]) -> Dict[str, Any]:
    metadata = {'caption': 'N/A', 'title': 'N/A', 'headers': [], 'rows': [], 'num_rows': 0, 'num_columns': 0, 'content_type': 'unknown'}
    descriptive_line = None
    for line in reversed(context_lines[-10:]):
        clean = line.strip()
        if not clean or clean.startswith('#'): continue
        caption_patterns = [r'(?i)^(?:table|tabella|tab\.|tbl\.?)\s*\d*[\.:]\s*(.+)', r'(?i)^(?:table|tabella)\s+\d+[\s\-—–]+(.+)', r'(?i)^caption[\s:]+(.+)']
        for pattern in caption_patterns:
            match = re.search(pattern, clean)
            if match:
                metadata['caption'] = match.group(1).strip()
                metadata['title'] = clean
                break
        if metadata['caption'] != 'N/A': break
        if len(clean) > 20 and not re.match(r'^[\|\-\s]+$', clean): descriptive_line = clean

    if metadata['caption'] == 'N/A' and descriptive_line: metadata['caption'] = descriptive_line[:150]
    header_row = None
    separator_idx = -1
    for idx, line in enumerate(table_lines):
        if re.match(r'^\s*\|[\s\-\|:]+\|\s*$', line):
            separator_idx = idx
            if idx > 0: header_row = table_lines[idx - 1]
            break
    if header_row:
        headers_raw = [cell.strip() for cell in header_row.split('|') if cell.strip()]
        metadata['headers'] = [h for h in headers_raw if h and not re.match(r'^[\-\s\|]+$', h)]
        metadata['num_columns'] = len(metadata['headers'])
    else:
        if len(table_lines) >= 2:
            first_row = table_lines[0]
            headers_raw = [cell.strip() for cell in first_row.split('|') if cell.strip() and len(cell.strip()) > 0 and not re.match(r'^[\-\s\|]+$', cell.strip())]
            if len(headers_raw) >= 2:
                metadata['headers'] = headers_raw[:2]
                metadata['num_columns'] = len(metadata['headers'])

    if separator_idx >= 0: data_rows = table_lines[separator_idx + 1:]
    else: data_rows = table_lines[1:] if len(table_lines) > 1 and metadata['headers'] else table_lines
    for line in data_rows:
        if not line.strip() or re.match(r'^\s*[\|\-\s]+$', line): continue
        cells_raw = line.split('|')
        cells = [cell.strip() for cell in cells_raw if cell.strip()]
        if len(cells) >= 1:
            while len(cells) < metadata['num_columns']: cells.append('')
            metadata['rows'].append(cells[:metadata['num_columns']])
    metadata['num_rows'] = len(metadata['rows'])

    if metadata['headers']:
        headers_text = ' '.join(metadata['headers']).lower()
        if any(kw in headers_text for kw in ['risk', 'hazard', 'danger', 'injury']): metadata['content_type'] = 'risk_assessment'
        elif any(kw in headers_text for kw in ['statistic', 'number', 'percentage', '%', 'rate']): metadata['content_type'] = 'statistical_data'
        elif any(kw in headers_text for kw in ['checklist', 'requirement', 'compliance']): metadata['content_type'] = 'compliance_checklist'
        elif any(kw in headers_text for kw in ['dx', 'a1', 'a2', 'read out', 'multiple answers']): metadata['content_type'] = 'survey_questionnaire'
    if metadata['content_type'] == 'unknown' and metadata['rows']:
        rows_text = ' '.join(' '.join(row).lower() for row in metadata['rows'])
        if any(kw in rows_text for kw in ['website', 'period', 'data source', 'years', 'covered by']): metadata['content_type'] = 'data_source_metadata'
        elif any(kw in rows_text for kw in ['workload', 'harassment', 'bullying', 'digital technologies', 'automation', 'emotional burden']): metadata['content_type'] = 'psychosocial_risks'
        else: metadata['content_type'] = 'general_data'
    return metadata

def format_table_content_readable(metadata: Dict[str, Any]) -> str:
    lines = []
    if metadata['headers']:
        header_line = ' | '.join([h if h else '' for h in metadata['headers']])
        lines.append(header_line)
        lines.append('-' * len(header_line))

    for row in metadata['rows']:
        if metadata['headers']:
            padded_cells = []
            for i, cell in enumerate(row):
                cell_str = str(cell)
                max_len = 1000 if re.match(r'^https?://', cell_str) else 500
                cell_str = cell_str[:max_len]
                header_len = len(metadata['headers'][i]) if i < len(metadata['headers']) and metadata['headers'][i] else len(cell_str)
                padded_cells.append(cell_str.ljust(header_len))
            lines.append(' | '.join(padded_cells))
        else:
            lines.append(' • ' + ' | '.join(str(c)[:500] for c in row))
    return '\n'.join(lines)

def analyze_document_structure_enhanced(doc) -> Dict[str, Any]:
    structure = {'num_tables': 0, 'num_figures': 0, 'num_headings': 0, 'headings': [], 'tables': [], 'figures': []}
    try:
        markdown = doc.export_to_markdown()
        lines = markdown.split('\n')
        for i, line in enumerate(lines):
            if line.strip().startswith('#'):
                level = len(line) - len(line.lstrip('#'))
                text = line.lstrip('#').strip()
                structure['num_headings'] += 1
                structure['headings'].append({'type': f'heading_level_{level}', 'text': text[:200], 'position': i, 'level': level})

        in_table = False
        current_table_lines = []
        table_start_idx = -1
        for i, line in enumerate(lines):
            if '|' in line and ('-' in line or '─' in line or ':' in line):
                if not in_table:
                    in_table = True
                    table_start_idx = i
                    current_table_lines = []
            if in_table:
                if '|' in line:
                    current_table_lines.append(line)
                else:
                    if current_table_lines and len(current_table_lines) >= 2:
                        context_start = max(0, table_start_idx - 15)
                        context_lines = lines[context_start:table_start_idx]
                        metadata = extract_table_metadata_enhanced(current_table_lines, context_lines)
                        readable_content = format_table_content_readable(metadata)
                        structure['num_tables'] += 1
                        structure['tables'].append({
                            'table_id': f'table_{structure["num_tables"]}',
                            'caption': metadata['caption'],
                            'title': metadata['title'],
                            'text_content': readable_content,
                            'raw_markdown': '\n'.join(current_table_lines),
                            'num_rows': metadata['num_rows'],
                            'num_columns': metadata['num_columns'],
                            'potential_columns': metadata['headers'],
                            'content_type': metadata['content_type'],
                            'position': table_start_idx,
                            'enhanced_parsing': True
                        })
                    in_table = False
                    current_table_lines = []

        for i, line in enumerate(lines):
            if '<image' in line.lower() or (line.strip().startswith('![') and '](' in line):
                structure['num_figures'] += 1
                caption = 'N/A'
                for j in range(i+1, min(len(lines), i+5)):
                    candidate = lines[j].strip()
                    if candidate and not candidate.startswith('#'):
                        caption = candidate
                        break
                structure['figures'].append({'figure_id': f'figure_{structure["num_figures"]}', 'caption': caption, 'position': i, 'vlm_description': None})
    except Exception as e:
        print(f"[!] Errore analisi: {e}")
    return structure

def export_tables_to_files_enhanced(tables: List[Dict], doc_id: str, tables_dir: Path) -> List[str]:
    if not tables: return []
    doc_tables_dir = tables_dir / doc_id
    doc_tables_dir.mkdir(parents=True, exist_ok=True)
    exported_tables = []
    for table in tables:
        table_id = table['table_id']
        txt_path = doc_tables_dir / f"{table_id}.txt"
        try:
            with open(txt_path, 'w', encoding='utf-8') as f:
                f.write("="*80 + "\n")
                f.write(f"TABELLA: {table_id}")
                if table.get('enhanced_parsing'): f.write(" (Enhanced Parsing)")
                f.write("\n" + "="*80 + "\n\n")
                f.write("█ METADATI TABELLA\n")
                f.write("-"*80 + "\n")
                f.write(f"Titolo/Caption: {table.get('title', 'N/A')}\n")
                f.write(f"Caption estratta: {table.get('caption', 'N/A')}\n")
                f.write(f"Tipo contenuto: {table.get('content_type', 'N/A')}\n")
                f.write(f"Righe dati: {table.get('num_rows', 'N/A')}\n")
                f.write(f"Colonne: {table.get('num_columns', len(table.get('potential_columns', [])))}\n")
                f.write(f"Posizione nel documento: {table.get('position', 'N/A')}\n\n")
                if table.get('potential_columns'):
                    f.write("█ COLONNE IDENTIFICATE:\n")
                    for idx, col in enumerate(table['potential_columns'], 1):
                        f.write(f"  {idx}. {col}\n")
                    f.write("\n")
                f.write("█ CONTENUTO (formattato):\n")
                f.write("-"*80 + "\n")
                f.write(table['text_content'])
                f.write("\n" + "-"*80 + "\n\n")
            exported_tables.append(str(txt_path))
            print(f"  [+] {txt_path.name} ({'enhanced' if table.get('enhanced_parsing') else 'standard'})")
        except Exception as e:
            print(f"  [!] Errore export tabella {table_id}: {e}")
    return exported_tables


# =========================================================
# CONFIGURAZIONE HPC
# =========================================================
class OSHAConfig:
    # PATH HPC - MODIFICATO
    OUTPUT_DIR = Path('/mnt/beegfs/home/fnatali/docling/OSHA_data')
    PDF_DIR = OUTPUT_DIR / 'pdfs'
    JSON_DIR = OUTPUT_DIR / 'json'
    IMAGES_DIR = OUTPUT_DIR / 'images'
    TABLES_DIR = OUTPUT_DIR / 'tables'
    KEYWORDS_DIR = OUTPUT_DIR / 'keywords'
    
    # VLM disabilitato di default per HPC (può essere lento)
    VLM_ENABLED = False
    VLM_MODEL = "HuggingFaceTB/SmolVLM-256M-Instruct"
    VLM_PROMPT = """Describe this technical document image in detail for academic research. Focus on:
- Main content type (diagram, chart, photo, schematic) - Key visual elements and their relationships
- Any visible text, labels, or numerical data - Technical details relevant to workplace safety or industrial processes
Be precise and comprehensive."""
    
    USER_AGENT = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    TIMEOUT = 30

    @classmethod
    def setup_directories(cls):
        for dir_path in [cls.OUTPUT_DIR, cls.PDF_DIR, cls.JSON_DIR, cls.IMAGES_DIR, cls.TABLES_DIR, cls.KEYWORDS_DIR]:
            dir_path.mkdir(exist_ok=True, parents=True)
        print(f"[INFO] Directory configurate: {cls.OUTPUT_DIR}")

    @classmethod
    def verify_directories(cls):
        all_exist = True
        for name, path in [('OUTPUT', cls.OUTPUT_DIR), ('PDFs', cls.PDF_DIR), ('JSON', cls.JSON_DIR), 
                          ('Images', cls.IMAGES_DIR), ('Tables', cls.TABLES_DIR), ('Keywords', cls.KEYWORDS_DIR)]:
            exists = path.exists()
            status = "✓" if exists else "✗"
            print(f"  {status} {name}: {path}")
            all_exist = all_exist and exists
        return all_exist

def initialize_docling_converter():
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
    print(f"[INFO] Docling OK | VLM: {'ON' if OSHAConfig.VLM_ENABLED else 'OFF'} | Enhanced Tables: ON\n")
    return converter


# =========================================================
# URL UTILITIES
# =========================================================
def normalize_query(q: str) -> str:
    q = unicodedata.normalize("NFKC", q).strip().lower()
    return " ".join(q.split()).replace("'", "'")

def build_osha_search_url(query: str = "", lang: str = "en", page: int = 0,
                          sort: str = "field_publication_date") -> str:
    if lang not in {"it", "en"}:
        raise ValueError("Il parametro 'lang' deve essere 'it' oppure 'en'.")
    if page < 0:
        raise ValueError("Pagina >= 0")
    base_url = f"https://osha.europa.eu/{lang}/publications"
    if query.strip():
        q = normalize_query(query)
        encoded_query = quote_plus(q)
        url = f"{base_url}?search_api_fulltext={encoded_query}&sort_by={sort}"
    else:
        url = f"{base_url}?sort_bef_combine=field_publication_date"
    if page > 0:
        separator = '&' if '?' in url else '?'
        url += f"{separator}page={page}"
    return url

# =========================================================
# ESTRAZIONE LINK CON REQUESTS
# =========================================================
def get_osha_publication_links_requests(
    start_page: int = 0,
    end_page: int = 4,
) -> List[Dict[str, str]]:
    """Estrae i link delle pubblicazioni OSHA utilizzando requests e BeautifulSoup."""

    publications = []
    seen_urls = set()
    print(f"\n{'='*80}")
    print(f"ESTRAZIONE LINK OSHA CON REQUESTS/BS4")
    print(f"{'='*80}")
    print(f"Range: pagina {start_page} → {end_page}")
    print(f"{'='*80}\n")
    headers = {'User-Agent': OSHAConfig.USER_AGENT}

    for page_num in range(start_page, end_page + 1):
        url = build_osha_search_url(page=page_num)
        print(f"[Pagina {page_num}] Fetching: {url}")

        if page_num > start_page:
            time.sleep(random.uniform(5, 8))

        try:
            response = requests.get(url, headers=headers, timeout=OSHAConfig.TIMEOUT)
            response.raise_for_status()

            soup = BeautifulSoup(response.content, 'lxml')
            page_pubs = []
            rows = soup.find_all('div', class_='revamp-row')

            if rows:
                print(f"  [+] Trovate {len(rows)} righe 'revamp-row'")
                for row in rows:
                    right_col = row.find('div', class_='publications-right-column')
                    h2 = right_col.find('h2') if right_col else row.find('h2')
                    if h2:
                        a_tag = h2.find('a')
                        if a_tag and a_tag.get('href'):
                            href = a_tag['href']
                            full_url = f"https://osha.europa.eu{href}" if href.startswith('/') else href
                            title = a_tag.get_text(strip=True)
                            if full_url not in seen_urls:
                                page_pubs.append({'titolo': title, 'url': full_url})
                                seen_urls.add(full_url)
            else:
                no_results = soup.find('div', class_='view-empty')
                if no_results and 'no results' in no_results.text.lower():
                    print("  [INFO] Raggiunta la fine dei risultati. Interrompo.")
                    break
                print("  [!] Nessun elemento 'revamp-row' trovato.")

            if page_pubs:
                print(f"  [FOUND] {len(page_pubs)} pubblicazioni in pagina {page_num}")
                publications.extend(page_pubs)
            else:
                print(f"  [!] Nessun link trovato in pagina {page_num}.")

        except requests.exceptions.RequestException as e:
            print(f"  [CRITICAL] Errore richiesta HTTP pagina {page_num}: {e}. Salto e continuo.")
            time.sleep(random.uniform(10, 15))
            continue
        except Exception as e:
            print(f"  [!] Errore generico pagina {page_num}: {e}. Salto e continuo.")
            time.sleep(random.uniform(10, 15))
            continue

    print(f"\n[RESULT] Totale: {len(publications)} pubblicazioni trovate")
    return publications

# =========================================================
# PDF PROCESSING E SCRAPING
# =========================================================
def download_pdf_osha(pdf_url: str, output_path: Path) -> bool:
    try:
        headers = {'User-Agent': OSHAConfig.USER_AGENT}
        response = requests.get(pdf_url, headers=headers, timeout=OSHAConfig.TIMEOUT, stream=True)
        response.raise_for_status()

        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return output_path.exists() and output_path.stat().st_size > 0
    except Exception as e:
        print(f"[!] Errore download PDF: {e}")
        return False

def export_images_from_document(doc, doc_id: str) -> List[Dict[str, Any]]:
    images_dir = OSHAConfig.IMAGES_DIR / doc_id
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

def add_vlm_descriptions_to_images(exported_images: List[Dict], doc_id: str) -> List[Dict]:
    if not exported_images or not OSHAConfig.VLM_ENABLED:
        return exported_images
    try:
        from transformers import AutoProcessor, AutoModelForVision2Seq
        from PIL import Image
        import torch
        print(f"\n[VLM] Caricamento {OSHAConfig.VLM_MODEL}...")
        processor = AutoProcessor.from_pretrained(OSHAConfig.VLM_MODEL, trust_remote_code=True)
        model = AutoModelForVision2Seq.from_pretrained(
            OSHAConfig.VLM_MODEL,
            torch_dtype=torch.float16 if torch.cuda.is_available() else torch.float32,
            device_map="auto",
            trust_remote_code=True
        )
        for idx, img_data in enumerate(exported_images, 1):
            try:
                print(f"  [{idx}/{len(exported_images)}] {img_data['filename']}")
                image = Image.open(img_data['path']).convert('RGB')
                messages = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": OSHAConfig.VLM_PROMPT}]}]
                prompt = processor.apply_chat_template(messages, add_generation_prompt=True)
                inputs = processor(text=prompt, images=[image], return_tensors="pt").to(model.device)
                with torch.no_grad():
                    outputs = model.generate(**inputs, max_new_tokens=300, do_sample=False)
                full_output = processor.decode(outputs[0], skip_special_tokens=True)
                description = full_output.split("Assistant:")[-1].strip() if "Assistant:" in full_output else full_output.strip()
                img_data['vlm_description'] = description
                img_path = Path(img_data['path'])
                txt_path = img_path.with_suffix('.txt')
                with open(txt_path, 'w', encoding='utf-8') as f:
                    f.write("="*80 + "\n")
                    f.write(f"VLM DESCRIPTION - {img_data['filename']}\n")
                    f.write("="*80 + "\n\n")
                    f.write("VLM ANALYSIS:\n")
                    f.write("-"*80 + "\n")
                    f.write(description)
                    f.write("\n" + "-"*80 + "\n")
                img_data['vlm_description_file'] = str(txt_path)
                print(f"    ✓ {description[:60]}...")
            except Exception as e:
                print(f"    ✗ {e}")
                img_data['vlm_description'] = None
                img_data['vlm_description_file'] = None

        del model, processor
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        success = sum(1 for img in exported_images if img.get('vlm_description'))
        print(f"[VLM] Completato: {success}/{len(exported_images)}\n")
    except Exception as e:
        print(f"[ERROR] VLM fallito: {e}")
    return exported_images

def scrape_osha_publication_requests(
    publication_url: str,
    converter,
    save_pdf: bool = True,
    enable_vlm: Optional[bool] = None
) -> Dict[str, Any]:
    use_vlm = enable_vlm if enable_vlm is not None else OSHAConfig.VLM_ENABLED
    headers = {'User-Agent': OSHAConfig.USER_AGENT}

    print(f"\n{'='*80}")
    print(f"SCRAPING PUBBLICAZIONE")
    print(f"{'='*80}\n")
    scraping_metadata = {
        'url': publication_url,
        'timestamp': datetime.now().isoformat(),
        'scraper_version': '7.0_HPC',
        'docling_used': True,
        'vlm_enabled': use_vlm
    }

    try:
        print(f"[INFO] Fetching metadata: {publication_url}")
        response = requests.get(publication_url, headers=headers, timeout=OSHAConfig.TIMEOUT)
        response.raise_for_status()
        page = BeautifulSoup(response.content, 'lxml')

        # Estrazione metadati
        title_elem = page.find('h1')
        title = title_elem.text.strip() if title_elem else 'N/A'

        descr_blocks = page.find_all(id=re.compile('^tmgmt-'))
        description = re.sub(r'\s+', ' ', ' '.join(p.get_text(strip=True, separator=' ') for p in descr_blocks).replace('\xa0', ' ')).strip()

        date_elem = page.find('p', class_='datetime-style')
        if not date_elem:
            date_elem = page.find('time')
        data_pubblicazione = date_elem.text.strip() if date_elem else 'N/A'

        pdf_url = None
        link_pdf_elem = page.find('div', class_='download-pdf')
        if link_pdf_elem and link_pdf_elem.find('a'):
            pdf_url = link_pdf_elem.find('a')['href']

        if not pdf_url:
            link_pdf_elem = page.find('a', href=re.compile(r'\.pdf$', re.IGNORECASE))
            if link_pdf_elem:
                pdf_url = link_pdf_elem['href']

        if pdf_url:
            pdf_url = urljoin(publication_url, pdf_url)

        # Estrazione Keywords
        keywords_elem = page.find('ul', class_='field__items')
        keywords = [kw.strip() for kw in keywords_elem.text.strip().split('\n')] if keywords_elem else []

        if keywords:
            safe_filename = re.sub(r'[^\w\-_.]', '_', title[:50])
            keywords_path = OSHAConfig.KEYWORDS_DIR / f"{safe_filename}_keywords.json"
            with open(keywords_path, 'w', encoding='utf-8') as f:
                json.dump({'url': publication_url, 'title': title, 'keywords': keywords}, f, ensure_ascii=False, indent=2)

        print(f"[+] {title[:60]}...")
        print(f"[+] Data: {data_pubblicazione}")

        pdf_data = None
        pdf_local_path = None
        if pdf_url:
            safe_filename = re.sub(r'[^\w\-_.]', '_', title[:50]) + '.pdf'
            pdf_local_path = OSHAConfig.PDF_DIR / safe_filename
            doc_id = safe_filename.replace('.pdf', '')

            if download_pdf_osha(pdf_url, pdf_local_path):
                print(f"[INFO] Processing PDF con enhanced table extraction...")

                result = converter.convert(str(pdf_local_path))
                doc = result.document

                markdown_text = doc.export_to_markdown()
                plain_text = doc.export_to_text()
                structure_info = analyze_document_structure_enhanced(doc)
                exported_images = export_images_from_document(doc, doc_id)

                if use_vlm and exported_images:
                    exported_images = add_vlm_descriptions_to_images(exported_images, doc_id)
                    for idx, fig_info in enumerate(structure_info.get('figures', [])):
                        if idx < len(exported_images) and exported_images[idx].get('vlm_description'):
                            fig_info['vlm_description'] = exported_images[idx]['vlm_description']

                exported_tables = export_tables_to_files_enhanced(structure_info['tables'], doc_id, OSHAConfig.TABLES_DIR)

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
                'abstract': description,
                'data_pubblicazione': data_pubblicazione,
                'pdf_url': pdf_url,
                'pdf_local_path': str(pdf_local_path) if pdf_local_path else None
            },
            'document_content': pdf_data,
            'status': 'success',
            'has_pdf': pdf_url is not None,
            'pdf_processed': pdf_data is not None
        }

    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Errore richiesta HTTP: {e}")
        return {'scraping_metadata': scraping_metadata, 'status': 'error', 'error_message': f'HTTP Request Error: {e}'}
    except Exception as e:
        print(f"[ERROR] Errore generico: {e}")
        return {'scraping_metadata': scraping_metadata, 'status': 'error', 'error_message': str(e)}


def save_to_json(data: Dict[str, Any], output_name: Optional[str] = None) -> Path:
    try:
        if output_name is None:
            title = data.get('web_metadata', {}).get('title', 'unknown')
            safe_title = re.sub(r'[^\w\-_.]', '_', title[:50])
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_name = f"{safe_title}_{timestamp}.json"
        output_path = OSHAConfig.JSON_DIR / output_name
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
    scraped_urls = set()
    json_dir = OSHAConfig.JSON_DIR
    if not json_dir.exists():
        return scraped_urls
    for json_file in list(json_dir.glob("*.json")):
        try:
            with open(json_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            docs = data if isinstance(data, list) else [data]
            for doc in docs:
                url = doc.get('scraping_metadata', {}).get('url')
                if url:
                    scraped_urls.add(url)
        except:
            continue
    return scraped_urls

# =========================================================
# BATCH SCRAPING INCREMENTALE
# =========================================================
def batch_scrape_osha_publications_incremental(
    converter,
    start_page: int = 0,
    end_page: int = 4,
    max_documents: Optional[int] = None,
    enable_vlm: bool = False,
    save_pdfs: bool = True,
    delay_between_docs: int = 4
) -> List[Dict[str, Any]]:

    already_scraped = load_already_scraped_urls()
    print(f"\n{'='*80}")
    print(f"SCRAPING INCREMENTALE OSHA")
    print(f"{'='*80}")
    print(f"Documenti già nel database: {len(already_scraped)}")
    print(f"Pagine: {start_page} → {end_page}")
    print(f"{'='*80}\n")

    publication_links = get_osha_publication_links_requests(
        start_page=start_page,
        end_page=end_page
    )

    if not publication_links:
        print("[WARNING] Nessuna pubblicazione trovata")
        return []

    new_links = [pub for pub in publication_links if pub['url'] not in already_scraped]
    print(f"\n{'='*80}")
    print(f"FILTRO DOCUMENTI")
    print(f"{'='*80}")
    print(f"Trovati: {len(publication_links)}")
    print(f"NUOVI: {len(new_links)}")
    print(f"GIÀ PROCESSATI (skip): {len(publication_links) - len(new_links)}")
    print(f"{'='*80}\n")

    if not new_links:
        print("✓ Nessun nuovo documento in questo range!")
        print(f"  Prova ad aumentare il range (es: {end_page+1}-{end_page+5})")
        return []

    print("📋 NUOVE PUBBLICAZIONI DA SCARICARE:")
    for i, pub in enumerate(new_links[:10], 1):
        print(f"  {i}. {pub['titolo'][:70]}...")

    if len(new_links) > 10:
        print(f"  ... e altre {len(new_links) - 10}\n")
    elif len(new_links) > 0:
        print()

    if max_documents:
        new_links = new_links[:max_documents]
        print(f"[INFO] Limitato a {max_documents} nuovi documenti\n")

    results = []
    failed = []
    print(f"{'='*80}")
    print(f"SCRAPING {len(new_links)} NUOVI DOCUMENTI")
    print(f"{'='*80}\n")

    for idx, pub in enumerate(new_links, 1):
        print(f"\n[{idx}/{len(new_links)}] {pub['titolo'][:60]}...")
        try:
            result = scrape_osha_publication_requests(
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
                    failed.append({'url': pub['url'], 'title': pub['titolo'], 'error': 'JSON save failed'})
            else:
                failed.append({'url': pub['url'], 'title': pub['titolo'], 'error': result.get('error_message', 'Unknown')})
                print(f"  ✗ FAIL")
        except KeyboardInterrupt:
            print(f"\nINTERRUZIONE - Salvati: {len(results)}")
            break
        except Exception as e:
            print(f"  ✗ Eccezione: {str(e)[:80]}")
            failed.append({'url': pub['url'], 'title': pub['titolo'], 'error': str(e)})

        if idx < len(new_links):
            delay = random.uniform(delay_between_docs, delay_between_docs + 2)
            time.sleep(delay)

    print(f"\n{'='*80}")
    print(f"COMPLETATO")
    print(f"{'='*80}")
    print(f"Nuovi salvati: {len(results)}")
    print(f"Falliti: {len(failed)}")
    print(f"Totale database: {len(already_scraped) + len(results)}")
    print(f"Path: {OSHAConfig.OUTPUT_DIR}")
    print(f"{'='*80}\n")
    if failed:
        print("\nDOCUMENTI FALLITI:")
        for f in failed[:5]:
            print(f"  • {f['title'][:50]}: {f['error']}")
        if len(failed) > 5:
            print(f"  ... e altri {len(failed)-5}")
    return results

# UTILITIES
def list_directory_contents():
    json_dir = OSHAConfig.JSON_DIR
    pdf_dir = OSHAConfig.PDF_DIR
    images_dir = OSHAConfig.IMAGES_DIR
    tables_dir = OSHAConfig.TABLES_DIR
    keywords_dir = OSHAConfig.KEYWORDS_DIR
    print(f"\n{'='*80}")
    print("STATISTICHE DIRECTORY OSHA")
    print(f"{'='*80}\n")
    if json_dir.exists():
        json_files = list(json_dir.glob("*.json"))
        total_size = sum(f.stat().st_size for f in json_files) / (1024*1024)
        print(f"JSON: {len(json_files)} file ({total_size:.1f} MB)")
        if json_files:
            recent = sorted(json_files, key=lambda x: x.stat().st_mtime, reverse=True)[:3]
            print("   Ultimi 3:")
            for f in recent:
                print(f"     • {f.name}")
    if pdf_dir.exists():
        pdf_files = list(pdf_dir.glob("*.pdf"))
        total_size = sum(f.stat().st_size for f in pdf_files) / (1024*1024)
        print(f"\nPDF: {len(pdf_files)} file ({total_size:.1f} MB)")
    if images_dir.exists():
        image_folders = [d for d in images_dir.iterdir() if d.is_dir()]
        total_images = sum(len(list(d.glob("*.png"))) for d in image_folders)
        total_txt = sum(len(list(d.glob("*.txt"))) for d in image_folders)
        print(f"\nImmagini: {total_images} PNG + {total_txt} descrizioni VLM")
        print(f"   in {len(image_folders)} documenti")
    if tables_dir.exists():
        table_folders = [d for d in tables_dir.iterdir() if d.is_dir()]
        total_tables = sum(len(list(d.glob("*.txt"))) for d in table_folders)
        print(f"\nTabelle: {total_tables} file")
        print(f"   in {len(table_folders)} documenti")
    if keywords_dir.exists():
        keyword_files = list(keywords_dir.glob("*_keywords.json"))
        print(f"\nKeywords: {len(keyword_files)} file")
    print(f"\n{'='*80}\n")

# =========================================================

# MAIN PER HPC

# =========================================================
def main():
    # ==========================================
    # CONFIGURAZIONE SCRAPING (modifica qui prima di lanciare)
    # ==========================================

    START_PAGE = 75           # Pagina di inizio
    NUM_PAGES = 140            # Numero di pagine da scansionare
    MAX_DOCS = None           # Limite documenti (None = tutti)
    ENABLE_VLM = True         # VLM attivo (True/False)
    SAVE_PDFS = True          # Salvare PDF localmente (True/False)
    DELAY_BETWEEN_DOCS = 4    # Secondi tra un documento e l'altro

    # ==========================================
    END_PAGE = START_PAGE + NUM_PAGES - 1

    print("\n" + "="*80)
    print("OSHA SCRAPER - HPC VERSION (AUTOMATICO)")
    print("="*80)
    print("\nFEATURES:")
    print("   • Adattato per HPC (no Google Drive)")
    print("   • Path: /mnt/beegfs/home/fnatali/docling/OSHA_data")
    print("   • Enhanced Table Parser")
    print("   • Skip documenti già processati")
    print("   • Scraping con Requests/BS4\n")

    print("[STEP 1] Setup directories...")
    OSHAConfig.setup_directories()
    print("\n[VERIFICA] Controllo path...")
    if not OSHAConfig.verify_directories():
        print("Alcune directory non esistono!")
        return None

    print("\n[STEP 2] Caricamento statistiche...")
    already_scraped = load_already_scraped_urls()
    print(f"Documenti già processati: {len(already_scraped)}\n")

    print("[STEP 3] Inizializzazione Docling...")
    converter = initialize_docling_converter()

    print("=" * 80)
    print("CONFIGURAZIONE AUTOMATICA")
    print("=" * 80)
    print(f"Range pagine: {START_PAGE} → {END_PAGE} ({NUM_PAGES} pagine)")
    print(f"Nuovi docs max: {MAX_DOCS if MAX_DOCS else '(tutti trovati)'}")
    print(f"VLM: {'ON' if ENABLE_VLM else 'OFF'}")
    print(f"Salva PDF: {'SI' if SAVE_PDFS else 'NO'}")
    print(f"Già nel database: {len(already_scraped)} documenti")
    print(f"Output: {OSHAConfig.JSON_DIR}")
    print("=" * 80 + "\n")

    print("AVVIO AUTOMATICO IN 3 SECONDI...\n")
    time.sleep(3)

    results = batch_scrape_osha_publications_incremental(
        converter=converter,
        start_page=START_PAGE,
        end_page=END_PAGE,
        max_documents=MAX_DOCS,
        enable_vlm=ENABLE_VLM,
        save_pdfs=SAVE_PDFS,
        delay_between_docs=DELAY_BETWEEN_DOCS
    )

    print("\n" + "="*80)
    print("SCRAPING TERMINATO")
    print("="*80)
    if results:
        print(f"✓ Nuovi documenti salvati: {len(results)}")

    total_docs = len(load_already_scraped_urls())
    print(f"Documenti totali nel database: {total_docs}")
    print(f"Path: {OSHAConfig.OUTPUT_DIR}")
    print("\n" + "="*80)
    print("STATISTICHE FINALI:")
    print("="*80)
    list_directory_contents()
    print("="*80 + "\n")
    print("SUGGERIMENTI:")
    print(f"   • Prossimo range: pagina {END_PAGE + 1} → {END_PAGE + NUM_PAGES}")
    print(f"   • Documenti totali: {total_docs}\n")
    return results


if __name__ == "__main__":

    results = main()
