#!/usr/bin/env python3

"""

✅ LAYER 2 HYBRID: RELATION EXTRACTION - FIXED & ALIGNED



ALLINEATO CON:

- layer0b: usa law_id come primary key

- layer1: CITES relationships esistenti

- layer1b: Article nodes con article_id



PIPELINE:

1. Structured: <ref> tags da AKN XML

2. Pattern Mining: Regex ABROGA/MODIFICA

3. LLM Refinement: Ollama post-processing (opzionale)

"""



import os

import re

from lxml import etree

from neo4j import GraphDatabase

from tqdm import tqdm

from collections import defaultdict

import argparse

import json



# ======================== CONFIG ========================

NORMATTIVA_DIR = "/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/normattiva_akn_filtered"

NEO4J_URI = os.getenv("NEO4J_URI", "bolt://127.0.0.1:7687")

NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")

NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "thesis2024")



BATCH_SIZE = 30

MAX_TEXT_LENGTH = 1_000_000



# LLM Settings

USE_LLM_REFINEMENT = False

OLLAMA_MODEL = "gemma2:27b-instruct-q4_0"

OLLAMA_HOST = os.getenv("OLLAMA_HOST", "127.0.0.1:11434")

OLLAMA_URL = f"http://{OLLAMA_HOST}/api/generate"



driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))



# ======================== NAMESPACES ========================

NS = {

    'akn': 'http://docs.oasis-open.org/legaldocml/ns/akn/3.0',

    'eli': 'http://data.europa.eu/eli/ontology#',

    'nakn': 'http://normattiva.it/akn/vocabulary'

}



# ======================== TIPO MAP ========================

TIPO_MAP = {

    'decreto_legislativo': 'dlgs',

    'decreto_legge': 'dl',

    'legge': 'legge',

    'decreto_del_presidente_della_repubblica': 'dpr',

    'regio_decreto': 'rd',

    'regio_decreto_legislativo': 'rdl'

}





class HybridRelationExtractor:

    

    def __init__(self, job_id='single', use_llm=False):

        self.driver = driver

        self.stats = defaultdict(int)

        self.job_id = job_id

        self.use_llm = use_llm

        

        self.relations = []

        self.llm_candidates = []

    

    def extract_metadata(self, xml_path):

        """

        ✅ Estrae metadata + genera law_id CANONICO

        Allineato con layer0b/1/1b

        """

        dirname = os.path.basename(os.path.dirname(xml_path))

        filename = os.path.basename(xml_path)

        

        # Parse dirname: decreto_legislativo_19480331_242

        dir_pattern = r'(.+?)_(\d{8})_(\d+)'

        dir_match = re.match(dir_pattern, dirname)

        if not dir_match:

            return None

        

        tipo_raw, date_str, num = dir_match.groups()

        tipo_clean = tipo_raw.lower().replace(' ', '_')

        tipo_short = TIPO_MAP.get(tipo_clean, 'legge')

        

        year = date_str[:4]

        month = date_str[4:6]

        day = date_str[6:8]

        

        # ✅ law_id CANONICO (stesso formato layer0b/1/1b)

        law_id = f"{tipo_short}-{num}-{year}"

        

        # Vigenza da filename

        file_pattern = r'(\d{4}-\d{2}-\d{2})_(.+?)_VIGENZA_(\d{4}-\d{2}-\d{2})_V0\.xml'

        file_match = re.match(file_pattern, filename)

        

        vigenza_end = None

        if file_match:

            _, _, vigenza_end = file_match.groups()

        

        return {

            'law_id': law_id,  # ✅ PRIMARY KEY

            'tipo_short': tipo_short,

            'num': num,

            'year': year,

            'date': f"{year}-{month}-{day}",

            'vigenza_start': f"{year}-{month}-{day}",

            'vigenza_end': vigenza_end

        }

    

    # ============ LAYER 2A: STRUCTURED (AKN <ref>) ============

    

    def extract_ref_citations(self, root, source_meta):

        """

        Estrae citazioni strutturate da <ref> tags

        ✅ Usa law_id per match Neo4j

        """

        for ref in root.findall('.//akn:ref', NS):

            href = ref.get('href')

            if not href:

                continue

            

            # Parse href: /akn/it/act/legge/stato/1949-02-28/43/!main#art_10

            match = re.search(

                r'/act/([^/]+)/[^/]+/(\d{4})-\d{2}-\d{2}/(\d+)/[^#]*(?:#art[_-](\d+))?',

                href

            )

            if not match:

                continue

            

            tipo_raw, year, num = match.groups()[:3]

            article = match.groups()[3] if len(match.groups()) > 3 else None

            

            tipo = TIPO_MAP.get(tipo_raw.replace('_', '').replace('.', ''), 'legge')

            target_law_id = f"{tipo}-{num}-{year}"

            

            # Confidence

            parent = ref.getparent()

            parent_tag = etree.QName(parent.tag).localname if parent is not None else ''

            

            if parent_tag in ['citations', 'citation']:

                confidence = 0.95

                method = 'akn_citation'

            else:

                confidence = 0.85

                method = 'akn_ref'

            

            self.relations.append({

                'type': 'RICHIAMA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': article,

                'confidence': confidence,

                'method': method

            })

            

            if article:

                self.stats['RICHIAMA_akn_ref_article'] += 1

            else:

                self.stats['RICHIAMA_akn_ref_law'] += 1

    

    # ============ LAYER 2B: PATTERN MINING ============

    

    def extract_pattern_abrogations(self, root, source_meta):

        """

        PATTERN MINING per ABROGA/SOPPRIME/DECADE

        """

        text_parts = []

        for p in root.findall('.//akn:p', NS):

            text_parts.append(etree.tostring(p, method='text', encoding='unicode'))

        

        full_text = " ".join(text_parts)

        

        if len(full_text) > MAX_TEXT_LENGTH:

            self.stats['text_too_long'] += 1

            return

        

        # Pattern 1: ((ARTICOLO ABROGATO DAL...))

        pattern1 = r'\(\((?:ARTICOLO|COMMA|LETTERA|PERIODO|PROVVEDIMENTO)?\s*(ABROGATO|SOPPRESSO)\s+(?:DAL|DALLA)?\s*(?:L\.|legge|D\.LGS|DLGS|D\.L\.|DL|DPR)\.?\s*(?:n\.?)?\s*(\d+)[,/\s]+(\d{4})\s*\)\)'

        

        for match in re.finditer(pattern1, full_text, re.IGNORECASE):

            rel_type_raw, num, year = match.groups()

            

            context = match.group(0)

            tipo = self._extract_tipo_from_context(context)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'ABROGA',

                'src_law_id': target_law_id,  # Chi abroga

                'tgt_law_id': source_meta['law_id'],  # Chi è abrogato

                'tgt_article': None,

                'confidence': 1.0,

                'method': 'pattern_annotation'

            })

            self.stats['ABROGA_pattern_annotation'] += 1

        

        # Pattern 2: "è abrogata/abrogato"

        pattern2 = r'(?:La|Il)\s+(legge|decreto\s+legislativo|D\.LGS|DPR)[^\d]{1,50}?n?\.?\s*(\d+)[^\d]{1,50}?(\d{4})[^\n]{1,100}?(?:è|e\'|sono)\s+abrogat[aeo]'

        

        for match in re.finditer(pattern2, full_text, re.IGNORECASE):

            tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'ABROGA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': None,

                'confidence': 0.9,

                'method': 'pattern_explicit'

            })

            self.stats['ABROGA_pattern_explicit'] += 1

        

        # Pattern 3: "sopprime l'articolo X"

        pattern3 = r'sopprim[eo][^\n]{1,30}?art(?:icolo)?\.?\s+(\d+)[^\n]{1,60}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})'

        

        for match in re.finditer(pattern3, full_text, re.IGNORECASE):

            article, tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'ABROGA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': article,

                'confidence': 0.92,

                'method': 'pattern_suppress'

            })

            self.stats['ABROGA_pattern_suppress'] += 1

        

        # Pattern 4: "decade/cessa di avere efficacia"

        pattern4 = r'(?:decade|cessa)[^\n]{1,80}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})'

        

        for match in re.finditer(pattern4, full_text, re.IGNORECASE):

            tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'ABROGA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': None,

                'confidence': 0.88,

                'method': 'pattern_lapse'

            })

            self.stats['ABROGA_pattern_lapse'] += 1

        

        # Pattern 5: "Sono abrogati gli articoli X, Y, Z"

        pattern5 = r'(?:Sono|È)\s+abrogat[ie]\s+(?:gli?\s+)?articol[io]\s+([\d,\s\-e]+)[^\n]{1,60}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})'

        

        for match in re.finditer(pattern5, full_text, re.IGNORECASE):

            articles_str, tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            articles = re.findall(r'\d+', articles_str)

            

            for article in articles[:10]:

                self.relations.append({

                    'type': 'ABROGA',

                    'src_law_id': source_meta['law_id'],

                    'tgt_law_id': target_law_id,

                    'tgt_article': article,

                    'confidence': 0.85,

                    'method': 'pattern_multi_repeal'

                })

                self.stats['ABROGA_pattern_multi'] += 1

    

    def extract_pattern_modifications(self, root, source_meta):

        """

        PATTERN MINING per MODIFICA/SOSTITUISCE/INTEGRA

        """

        text_parts = []

        for p in root.findall('.//akn:p', NS):

            text_parts.append(etree.tostring(p, method='text', encoding='unicode'))

        

        full_text = " ".join(text_parts)

        

        if len(full_text) > MAX_TEXT_LENGTH:

            return

        

        # Pattern 1: "art. X è modificato"

        pattern1 = r'art(?:icolo)?\.?\s+(\d+)[^\n]{1,60}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS|DPR)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})[^\n]{1,100}?(?:è|e\'|sono)\s+(?:modificat|sostitu)[aeio]'

        

        for match in re.finditer(pattern1, full_text, re.IGNORECASE):

            article, tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'MODIFICA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': article,

                'confidence': 0.9,

                'method': 'pattern_explicit_mod'

            })

            self.stats['MODIFICA_pattern_explicit'] += 1

        

        # Pattern 2: ((MODIFICATO...))

        pattern2 = r'\(\((?:ARTICOLO|COMMA)?\s*MODIFICATO\s+(?:DAL|DALLA)?\s*(?:L\.|legge|D\.LGS|DLGS|D\.L\.|DL)\.?\s*(?:n\.?)?\s*(\d+)[,/\s]+(\d{4})\s*\)\)'

        

        for match in re.finditer(pattern2, full_text, re.IGNORECASE):

            num, year = match.groups()

            

            context = match.group(0)

            tipo = self._extract_tipo_from_context(context)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'MODIFICA',

                'src_law_id': target_law_id,

                'tgt_law_id': source_meta['law_id'],

                'tgt_article': None,

                'confidence': 1.0,

                'method': 'pattern_annotation_mod'

            })

            self.stats['MODIFICA_pattern_annotation'] += 1

        

        # Pattern 3: "sostituisce l'articolo X"

        pattern3 = r'sostituisc[eo][^\n]{1,30}?art(?:icolo)?\.?\s+(\d+)[^\n]{1,60}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})'

        

        for match in re.finditer(pattern3, full_text, re.IGNORECASE):

            article, tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'MODIFICA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': article,

                'confidence': 0.85,

                'method': 'pattern_substitution'

            })

            self.stats['MODIFICA_pattern_substitution'] += 1

        

        # Pattern 4: "aggiorna/integra l'art. X"

        pattern4 = r'(?:aggiorn|integr|rettific)[ao][^\n]{1,30}?art(?:icolo)?\.?\s+(\d+)[^\n]{1,60}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})'

        

        for match in re.finditer(pattern4, full_text, re.IGNORECASE):

            article, tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'MODIFICA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': article,

                'confidence': 0.85,

                'method': 'pattern_update'

            })

            self.stats['MODIFICA_pattern_update'] += 1

        

        # Pattern 5: "Al comma X...sono aggiunte"

        pattern5 = r'(?:Al|Alla|All\')\s+(?:comma|articolo|art\.?)\s+(\d+)[^\n]{1,80}?(?:del|della)\s+(decreto\s+legislativo|legge|D\.LGS|DLGS)\.?\s*n?\.?\s*(\d+)[/\s,]+(\d{4})[^\n]{1,100}?(?:aggiunt|sostitu|modific)'

        

        for match in re.finditer(pattern5, full_text, re.IGNORECASE):

            article, tipo_raw, num, year = match.groups()

            tipo = self._normalize_tipo(tipo_raw)

            target_law_id = f"{tipo}-{num}-{year}"

            

            self.relations.append({

                'type': 'MODIFICA',

                'src_law_id': source_meta['law_id'],

                'tgt_law_id': target_law_id,

                'tgt_article': article,

                'confidence': 0.88,

                'method': 'pattern_amendment'

            })

            self.stats['MODIFICA_pattern_amendment'] += 1

    

    def _normalize_tipo(self, tipo_raw):

        """Normalizza tipo legge"""

        tipo_clean = tipo_raw.upper().replace('.', '').replace(' ', '')

        

        if 'DLGS' in tipo_clean or 'LEGISLATIVO' in tipo_clean:

            return 'dlgs'

        elif 'DL' in tipo_clean and 'DLGS' not in tipo_clean:

            return 'dl'

        elif 'DPR' in tipo_clean:

            return 'dpr'

        

        return 'legge'

    

    def _extract_tipo_from_context(self, context):

        """Estrai tipo da contesto"""

        context_upper = context.upper()

        

        if 'D.LGS' in context_upper or 'DLGS' in context_upper:

            return 'dlgs'

        elif 'D.L' in context_upper and 'D.LGS' not in context_upper:

            return 'dl'

        elif 'DPR' in context_upper:

            return 'dpr'

        

        return 'legge'

    

    # ============ MAIN PROCESSING ============

    

    def process_file(self, xml_path):

        """

        Processa singolo file XML AKN

        """

        try:

            meta = self.extract_metadata(xml_path)

            if not meta:

                self.stats['metadata_error'] += 1

                return

            

            tree = etree.parse(xml_path)

            root = tree.getroot()

            

            # EXTRACTION PIPELINE

            self.extract_ref_citations(root, meta)

            self.extract_pattern_abrogations(root, meta)

            self.extract_pattern_modifications(root, meta)

            

            self.stats['files_ok'] += 1

        

        except Exception as e:

            self.stats['file_error'] += 1

            if self.stats['file_error'] < 5:

                print(f"\n⚠️ {os.path.basename(xml_path)}: {str(e)[:100]}")

    

    def flush_relations(self):
        """
        ✅ GRANULARITÀ IBRIDA CORRETTA:
        - Article → Article specifico (quando href ha #art_N)
        - Article → Law (quando href NON ha #art_N)  
        - Law → Law (per modifiche intere)
        """
    
        if not self.relations:
            return
    
        try:
            with self.driver.session() as s:
                by_type = defaultdict(list)
                for rel in self.relations:
                    by_type[rel['type']].append(rel)
            
                for rel_type, rels in by_type.items():
                
                    # ===============================================
                    # CASO 1: Relazioni con ARTICLE TARGET (granulari)
                    # ===============================================
                    with_article = [r for r in rels if r.get('tgt_article')]
                
                    if with_article:
                        # Raggruppa per source_law per evitare duplicati
                        by_source = defaultdict(list)
                        for r in with_article:
                            by_source[r['src_law_id']].append(r)
                    
                        for src_law_id, src_rels in by_source.items():
                            # Per ogni legge source, usa PRIMO articolo
                            result = s.run(f"""
                                // ✅ Match SOURCE law + PRIMO articolo
                                MATCH (src_law:CanonicalLaw {{law_id: $src_law_id}})
                                MATCH (src_law)-[:HAS_ARTICLE]->(src_art:Article)
                                WITH src_art
                                ORDER BY src_art.article_num
                                LIMIT 1
                            
                                // ✅ Per ogni target article
                                UNWIND $rels as rel
                            
                                // Match TARGET article specifico
                                MATCH (tgt_law:CanonicalLaw {{law_id: rel.tgt_law_id}})
                                MATCH (tgt_law)-[:HAS_ARTICLE]->(tgt_art:Article)
                                WHERE tgt_art.article_num = rel.tgt_article
                            
                                // Crea relazione Article → Article (1 per citazione)
                                MERGE (src_art)-[r:{rel_type}]->(tgt_art)
                                ON CREATE SET 
                                    r.confidence = rel.confidence,
                                    r.method = rel.method,
                                    r.layer = 2,
                                    r.created_at = datetime()
                                ON MATCH SET
                                    r.confidence = CASE 
                                        WHEN rel.confidence > r.confidence 
                                        THEN rel.confidence 
                                        ELSE r.confidence 
                                    END
                            
                                RETURN count(r) as cnt
                            """, src_law_id=src_law_id, rels=src_rels)
                        
                            rec = result.single()
                            if rec and rec['cnt']:
                                self.stats[f'{rel_type}_exact'] += rec['cnt']
                
                    # ===============================================
                    # CASO 2: Relazioni LAW-LEVEL (senza articolo target)
                    # ===============================================
                    without_article = [r for r in rels if not r.get('tgt_article')]
                
                    if without_article:
                        # Raggruppa per source_law
                        by_source = defaultdict(list)
                        for r in without_article:
                            by_source[r['src_law_id']].append(r)
                    
                        for src_law_id, src_rels in by_source.items():
                            result = s.run(f"""
                                // ✅ Match SOURCE law + PRIMO articolo
                                MATCH (src_law:CanonicalLaw {{law_id: $src_law_id}})
                                MATCH (src_law)-[:HAS_ARTICLE]->(src_art:Article)
                                WITH src_art
                                ORDER BY src_art.article_num
                                LIMIT 1
                            
                                // ✅ Per ogni target law
                                UNWIND $rels as rel
                            
                                // Match TARGET law
                                MATCH (tgt_law:CanonicalLaw {{law_id: rel.tgt_law_id}})
                            
                                // Crea relazione Article → Law (citazione generica)
                                MERGE (src_art)-[r:{rel_type}]->(tgt_law)
                                ON CREATE SET
                                    r.confidence = rel.confidence,
                                    r.method = rel.method,
                                    r.layer = 2,
                                    r.created_at = datetime()
                                ON MATCH SET
                                    r.confidence = CASE 
                                        WHEN rel.confidence > r.confidence 
                                        THEN rel.confidence 
                                        ELSE r.confidence 
                                    END
                            
                                RETURN count(r) as cnt
                            """, src_law_id=src_law_id, rels=src_rels)
                        
                            rec = result.single()
                            if rec and rec['cnt']:
                                self.stats[f'{rel_type}_law'] += rec['cnt']
    
        except Exception as e:
            print(f"\n❌ Relations error: {str(e)[:200]}")
    
        self.relations = []

    

    def run(self, start=0, end=None):

        """

        Main execution

        """

        # Trova tutti gli XML

        xml_files = []

        for cat in os.listdir(NORMATTIVA_DIR):

            cat_path = os.path.join(NORMATTIVA_DIR, cat)

            if not os.path.isdir(cat_path):

                continue

            

            for law_dir in os.listdir(cat_path):

                law_path = os.path.join(cat_path, law_dir)

                if not os.path.isdir(law_path):

                    continue

                

                for f in os.listdir(law_path):

                    if f.endswith('.xml'):

                        xml_files.append(os.path.join(law_path, f))

        

        total = len(xml_files)

        print(f"📂 Total XML files: {total:,}")

        

        if end:

            xml_files = xml_files[start:end]

            print(f"📊 Processing: {start:,} → {end:,} ({len(xml_files):,} files)\n")

        

        # Process

        for i, xml_path in enumerate(tqdm(xml_files, desc=f"Layer2-{self.job_id}")):

            self.process_file(xml_path)

            

            # Flush periodico

            if i > 0 and i % BATCH_SIZE == 0:

                self.flush_relations()

        

        # Final flush

        self.flush_relations()

        

        self._print_stats()

    

    def _print_stats(self):

        """Print final statistics"""

        print("\n" + "="*70)

        print(f"✅ LAYER 2 HYBRID - JOB {self.job_id} DONE")

        print("="*70)

        print(f"Files processed: {self.stats['files_ok']:,}")

        print(f"Files errors: {self.stats['file_error']:,}")

        print("")

        

        abroga_total = sum(v for k, v in self.stats.items() 

                          if 'ABROGA' in k and ('_exact' in k or '_law' in k))

        modifica_total = sum(v for k, v in self.stats.items() 

                            if 'MODIFICA' in k and ('_exact' in k or '_law' in k))

        richiama_total = sum(v for k, v in self.stats.items() 

                            if 'RICHIAMA' in k and ('_exact' in k or '_law' in k))

        

        print(f"🎯 SUMMARY:")

        print(f"  ABROGA total: {abroga_total:,}")

        print(f"  MODIFICA total: {modifica_total:,}")

        print(f"  RICHIAMA total: {richiama_total:,}")

        print(f"\n📊 DETAILS:")

        

        for key in sorted(self.stats.keys()):

            if any(x in key for x in ['_exact', '_law', 'pattern', 'akn_ref']):

                print(f"  {key}: {self.stats[key]:,}")

        

        print("="*70)

    

    def close(self):

        self.driver.close()





if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--start', type=int, default=0)

    parser.add_argument('--end', type=int, default=None)

    parser.add_argument('--job-id', type=str, default='single')

    

    args = parser.parse_args()

    

    print(f"🚀 LAYER 2: HYBRID RELATION EXTRACTION (FIXED)")

    print(f"📊 Range: {args.start} → {args.end if args.end else 'ALL'}")

    print()
    
    extractor = HybridRelationExtractor(job_id=args.job_id, use_llm=USE_LLM_REFINEMENT)

    extractor.run(start=args.start, end=args.end)

    extractor.close()



    print("\n✅ Layer 2 Complete!")

 
