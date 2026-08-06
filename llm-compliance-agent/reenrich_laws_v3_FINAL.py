#!/usr/bin/env python3
"""
🔄 RE-ENRICH LAWS V3 - FINAL CORRECTED VERSION

FIXES CRITICI:
✅ NON deduplica per law_id - tiene TUTTE le citazioni con articoli/commi diversi
✅ Cerca articolo/comma PRIMA E DOPO la legge (es. "art. 4 del D.lgs 81/2008")
✅ Deduplica solo citazioni IDENTICHE (stesso law_id + article + comma + lettera)
✅ Aggiunge contesto per ogni citazione
✅ Salva in enriched.laws_structured E semantic_metadata.laws_structured
✅ Aggiunge position per ordinare le citazioni nel documento

Esempio CORRETTO:
- "art. 4 del D.lgs 81/2008"           → dlgs-81-2008-art4
- "D.lgs 81/2008, art. 37, comma 2"    → dlgs-81-2008-art37-comma2
- "D.lgs 81/2008"                       → dlgs-81-2008
- "art. 4, comma 3, D.lgs 81/2008"     → dlgs-81-2008-art4-comma3

TUTTE queste sono DIVERSE e vengono TUTTE mantenute!
"""

import json
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from tqdm import tqdm
import argparse
from datetime import datetime


# ======================== CONFIG ========================
INAIL_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/INAIL_enriched")
OSHA_ENRICHED = Path("/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/data/OSHA_enriched")

CONTEXT_WINDOW = 200  # Caratteri di contesto prima e dopo la citazione
SEARCH_WINDOW_BEFORE = 80  # Caratteri PRIMA della legge dove cercare articolo
SEARCH_WINDOW_AFTER = 120  # Caratteri DOPO la legge dove cercare articolo


# ======================== ENHANCED EXTRACTOR V3 ========================

class EnhancedLawExtractorV3:
    """
    Estrae citazioni legali SENZA deduplicazione aggressiva.
    TIENE TUTTE le citazioni con articoli/commi diversi!
    
    NOVITÀ V3:
    - Cerca articolo/comma PRIMA e DOPO la legge
    - Gestisce pattern come "art. 4 del D.lgs 81/2008"
    - Migliore estrazione di contesto
    """
    
    # Pattern per tipi di legge
    PATTERNS = {
        'decreto.legislativo': [
            # D.Lgs, d.lgs, dlgs, etc.
            r'[Dd]\.?\s*[Ll]\.?\s*[Gg]\.?\s*[Ss]?\.?\s*(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'[Dd]ecreto\s+[Ll]egislativo\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'[Dd]\.?\s*[Ll]gs\.?\s*(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
        ],
        'decreto.del.presidente.della.repubblica': [
            r'[Dd]\.?\s*[Pp]\.?\s*[Rr]\.?\s*(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'[Dd]ecreto\s+(?:del\s+)?[Pp]residente\s+(?:della\s+)?[Rr]epubblica\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
        ],
        'legge': [
            r'[Ll]egge\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'\b[Ll]\.?\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
        ],
        'decreto.legge': [
            r'[Dd]\.?\s*[Ll]\.?\s*(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'[Dd]ecreto[\s\-][Ll]egge\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
        ],
        'decreto.ministeriale': [
            r'[Dd]\.?\s*[Mm]\.?\s+(?:del\s+)?(\d{1,2})[/\-\.\s]+(\d{1,2})[/\-\.\s]+(\d{2,4})',
            r'[Dd]ecreto\s+[Mm]inisteriale\s+(?:del\s+)?(\d{1,2})[/\-\.\s]+(\d{1,2})[/\-\.\s]+(\d{2,4})',
        ],
        'direttiva.ue': [
            r'[Dd]irettiva\s+(?:UE\s+|CE\s+|\(UE\)\s+|\(CE\)\s+)?(\d{4})\s*[/\-]\s*(\d+)',
            r'[Dd]ir(?:ettiva)?\.?\s+(\d{4})\s*[/\-]\s*(\d+)',
        ],
        'regolamento.ue': [
            r'[Rr]egolamento\s+(?:UE\s+|CE\s+|\(UE\)\s+|\(CE\)\s+)?(?:n\.?\s*)?(\d+)\s*[/\-]\s*(\d{4})',
        ],
        'testo.unico': [
            r'[Tt]esto\s+[Uu]nico\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'[Tt]\.?\s*[Uu]\.?\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
        ],
        'regio.decreto': [
            r'[Rr]\.?\s*[Dd]\.?\s*(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
            r'[Rr]egio\s+[Dd]ecreto\s+(?:n\.?\s*)?(\d+)\s*[/\-°]\s*(\d{2,4})',
        ],
        # US Laws (for OSHA)
        'cfr': [
            r'(\d+)\s*[Cc]\.?\s*[Ff]\.?\s*[Rr]\.?\s*(?:§|[Pp]art)?\s*(\d+)',
            r'[Cc]ode\s+of\s+[Ff]ederal\s+[Rr]egulations.*?(\d+)\s*[Cc][Ff][Rr]\s*(\d+)',
        ],
        'usc': [
            r'(\d+)\s*[Uu]\.?\s*[Ss]\.?\s*[Cc]\.?\s*(?:§)?\s*(\d+)',
        ],
        'osha.standard': [
            r'OSHA\s+[Ss]tandard\s+(\d+)\.(\d+)',
            r'29\s*CFR\s*(?:§)?\s*(\d+)\.(\d+)',
        ],
    }
    
    # Pattern per articolo/comma/lettera - cerca sia PRIMA che DOPO
    ARTICLE_PATTERNS = [
        r'art(?:icolo|\.)?[\s\.]*(\d+(?:[\-_]?(?:bis|ter|quater|quinquies|sexies|septies|octies|novies|decies))?)',
        r'articolo[\s\.]*(\d+)',
        # US style
        r'[Ss]ection[\s\.]*(\d+(?:\.\d+)?)',
        r'§[\s]*(\d+(?:\.\d+)?)',
    ]
    
    COMMA_PATTERNS = [
        r'comma[\s\.]*(\d+(?:[\-_]?(?:bis|ter|quater))?)',
        r'co\.[\s]*(\d+)',
        # US style
        r'paragraph[\s\.]*\(?([a-z\d]+)\)?',
        r'subsection[\s\.]*\(?([a-z\d]+)\)?',
    ]
    
    LETTERA_PATTERNS = [
        r'lett(?:era)?[\s\.]*([a-z](?:[\-_]?(?:bis|ter))?)',
        r'lettera[\s\.]*([a-z])',
    ]
    
    def __init__(self):
        self.compiled_patterns = {}
        for tipo, patterns in self.PATTERNS.items():
            self.compiled_patterns[tipo] = [
                re.compile(p, re.IGNORECASE | re.MULTILINE) for p in patterns
            ]
        
        self.article_patterns = [re.compile(p, re.IGNORECASE) for p in self.ARTICLE_PATTERNS]
        self.comma_patterns = [re.compile(p, re.IGNORECASE) for p in self.COMMA_PATTERNS]
        self.lettera_patterns = [re.compile(p, re.IGNORECASE) for p in self.LETTERA_PATTERNS]
    
    
    def extract_all_citations(self, text: str, max_laws: int = 100) -> List[Dict]:
        """
        Estrae TUTTE le citazioni, incluse varianti con articoli/commi diversi.
        NON deduplica per law_id - deduplica solo per full_citation!
        
        Returns: Lista di citazioni ordinate per posizione nel testo
        """
        if not text or len(text) < 10:
            return []
        
        found_citations = []
        
        for tipo, compiled_patterns in self.compiled_patterns.items():
            for pattern in compiled_patterns:
                for match in pattern.finditer(text):
                    # Contesto PRIMA della legge (per "art. 4 del D.lgs 81/2008")
                    start_before = max(0, match.start() - SEARCH_WINDOW_BEFORE)
                    context_before = text[start_before:match.start()]
                    
                    # Contesto DOPO la legge (per "D.lgs 81/2008, art. 37")
                    end_after = min(len(text), match.end() + SEARCH_WINDOW_AFTER)
                    context_after = text[match.end():end_after]
                    
                    # Contesto completo per output
                    start_full = max(0, match.start() - CONTEXT_WINDOW)
                    end_full = min(len(text), match.end() + CONTEXT_WINDOW)
                    full_context = text[start_full:end_full].strip()
                    # Pulisci il contesto
                    full_context = re.sub(r'\s+', ' ', full_context)
                    
                    citation = self._parse_citation(
                        tipo=tipo,
                        match_groups=match.groups(),
                        context_before=context_before,
                        context_after=context_after,
                        full_match=match.group(0),
                        full_context=full_context,
                        position=match.start()
                    )
                    
                    if citation:
                        found_citations.append(citation)
        
        # ✅ DEDUPLICA SOLO CITAZIONI IDENTICHE (full_citation + position vicine)
        deduplicated = self._deduplicate_smart(found_citations)
        
        # Ordina per posizione nel testo
        deduplicated.sort(key=lambda x: x.get('position', 0))
        
        return deduplicated[:max_laws]
    
    
    def _parse_citation(
        self,
        tipo: str,
        match_groups: tuple,
        context_before: str,
        context_after: str,
        full_match: str,
        full_context: str,
        position: int
    ) -> Optional[Dict]:
        """Parsa singola citazione cercando articolo/comma PRIMA E DOPO"""
        
        try:
            # Parse numero e anno in base al tipo
            if tipo == 'decreto.ministeriale' and len(match_groups) == 3:
                numero = f"{match_groups[0]}-{match_groups[1]}"
                anno = self._normalize_year(match_groups[2])
            elif tipo in ['direttiva.ue'] and len(match_groups) == 2:
                anno = match_groups[0]
                numero = match_groups[1]
            elif tipo == 'regolamento.ue' and len(match_groups) == 2:
                numero = match_groups[0]
                anno = match_groups[1]
            elif tipo in ['cfr', 'usc', 'osha.standard'] and len(match_groups) == 2:
                numero = match_groups[0]
                anno = match_groups[1]  # Per CFR, questo è la parte/sezione
            else:
                if len(match_groups) < 2:
                    return None
                numero = match_groups[0]
                anno = self._normalize_year(match_groups[1])
            
            tipo_short = self._get_tipo_short(tipo)
            
            # Per leggi US, formato diverso
            if tipo in ['cfr', 'usc', 'osha.standard']:
                law_id = f"{tipo_short}-{numero}-{anno}"
            else:
                law_id = f"{tipo_short}-{numero}-{anno}"
            
            # ✅ CERCA ARTICOLO/COMMA SIA PRIMA CHE DOPO LA LEGGE
            article_num = None
            comma_num = None
            lettera = None
            
            # Prima cerca PRIMA (es. "art. 4 del D.lgs")
            article_num = self._find_in_context(self.article_patterns, context_before)
            if not article_num:
                # Poi cerca DOPO (es. "D.lgs, art. 4")
                article_num = self._find_in_context(self.article_patterns, context_after)
            
            # Stessa cosa per comma
            comma_num = self._find_in_context(self.comma_patterns, context_before)
            if not comma_num:
                comma_num = self._find_in_context(self.comma_patterns, context_after)
            
            # E lettera
            lettera = self._find_in_context(self.lettera_patterns, context_before)
            if not lettera:
                lettera = self._find_in_context(self.lettera_patterns, context_after)
            
            # Costruisci full_citation UNICA
            full_citation_parts = [law_id]
            if article_num:
                full_citation_parts.append(f"art{article_num}")
            if comma_num:
                full_citation_parts.append(f"comma{comma_num}")
            if lettera:
                full_citation_parts.append(f"lett{lettera}")
            
            full_citation = "-".join(full_citation_parts)
            
            # URN
            urn = self._generate_urn(tipo, numero, anno)
            
            return {
                'law_id': law_id,
                'tipo': tipo,
                'tipo_short': tipo_short,
                'numero': str(numero),
                'anno': str(anno),
                'stringa_originale': full_match.strip(),
                'article_num': article_num,
                'comma_num': comma_num,
                'lettera': lettera,
                'urn': urn,
                'qdrant_prefix': law_id,
                'full_citation': full_citation,
                'jurisdiction': 'IT' if tipo not in ['cfr', 'usc', 'osha.standard'] else 'US',
                'context': full_context[:300],  # Limita contesto
                'position': position
            }
        
        except Exception as e:
            return None
    
    
    def _find_in_context(self, patterns: List[re.Pattern], context: str) -> Optional[str]:
        """Cerca un pattern nel contesto e restituisce il primo match"""
        for pattern in patterns:
            match = pattern.search(context)
            if match:
                return match.group(1).lower()
        return None
    
    
    def _deduplicate_smart(self, citations: List[Dict]) -> List[Dict]:
        """
        Deduplica intelligente:
        - Rimuove citazioni IDENTICHE (stesso full_citation) se sono vicine (<100 chars)
        - TIENE citazioni con stesso law_id ma articoli/commi diversi!
        - TIENE citazioni identiche se sono lontane nel testo (citazioni ripetute)
        """
        if not citations:
            return []
        
        # Ordina per posizione
        citations.sort(key=lambda x: x.get('position', 0))
        
        result = []
        seen_at_position = {}  # full_citation -> ultima posizione vista
        
        PROXIMITY_THRESHOLD = 200  # Considera "vicine" se entro 200 caratteri
        
        for citation in citations:
            key = citation['full_citation']
            pos = citation.get('position', 0)
            
            if key in seen_at_position:
                last_pos = seen_at_position[key]
                # Se la stessa citazione è lontana, tienila (è una citazione ripetuta nel doc)
                if pos - last_pos > PROXIMITY_THRESHOLD:
                    result.append(citation)
                    seen_at_position[key] = pos
                # Altrimenti skip (è un duplicato vicino, probabilmente stesso paragrafo)
            else:
                result.append(citation)
                seen_at_position[key] = pos
        
        return result
    
    
    def _get_tipo_short(self, tipo: str) -> str:
        """Mappa tipo → tipo_short"""
        mapping = {
            'decreto.legislativo': 'dlgs',
            'decreto.del.presidente.della.repubblica': 'dpr',
            'legge': 'legge',
            'decreto.legge': 'dl',
            'decreto.ministeriale': 'dm',
            'direttiva.ue': 'dir',
            'regolamento.ue': 'reg',
            'testo.unico': 'tu',
            'regio.decreto': 'rd',
            'cfr': 'cfr',
            'usc': 'usc',
            'osha.standard': 'osha',
        }
        return mapping.get(tipo, tipo.replace('.', '_'))
    
    
    def _normalize_year(self, year_str: str) -> str:
        """Normalizza anno a 4 cifre"""
        try:
            year = int(year_str)
            if year < 50:
                return str(2000 + year)
            elif year < 100:
                return str(1900 + year)
            else:
                return str(year)
        except:
            return str(year_str)
    
    
    def _generate_urn(self, tipo: str, numero: str, anno: str) -> str:
        """Genera URN NIR (solo per leggi italiane)"""
        if tipo in ['cfr', 'usc', 'osha.standard']:
            return f"urn:us:{tipo}:{numero}:{anno}"
        
        data_placeholder = f"{anno}-01-01"
        return f"urn:nir:stato:{tipo}:{data_placeholder};{numero}"


# ======================== MAIN PROCESSING ========================

def process_json_file(json_path: Path, extractor: EnhancedLawExtractorV3, dry_run: bool = False) -> Tuple[bool, int, int]:
    """
    Processa un singolo file JSON:
    1. Legge il testo
    2. Estrae TUTTE le citazioni (senza dedup aggressiva)
    3. Salva in enriched.laws_structured E semantic_metadata.laws_structured
    
    Returns: (success, num_citations_new, num_citations_old)
    """
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Conta citazioni vecchie
        old_laws = data.get('semantic_metadata', {}).get('laws_structured', [])
        num_old = len(old_laws)
        
        # Estrai testo
        content = data.get('document_content', {})
        text = content.get('markdown_content') or content.get('plain_text', '')
        
        if not text or len(text) < 100:
            return False, 0, num_old
        
        # ✅ Estrai TUTTE le citazioni con nuovo estrattore
        citations = extractor.extract_all_citations(text, max_laws=100)
        
        if dry_run:
            return True, len(citations), num_old
        
        # ✅ Salva in enriched
        if 'enriched' not in data:
            data['enriched'] = {}
        
        data['enriched']['laws_structured'] = citations
        data['enriched']['laws_count'] = len(citations)
        data['enriched']['unique_laws'] = len(set(c['law_id'] for c in citations))
        data['enriched']['reenriched_at'] = datetime.now().isoformat()
        data['enriched']['extractor_version'] = 'EnhancedLawExtractorV3_KeepAllVariants'
        data['enriched']['laws_confidence'] = _calculate_confidence(citations)
        
        # ✅ Aggiorna ANCHE semantic_metadata per compatibilità
        if 'semantic_metadata' not in data:
            data['semantic_metadata'] = {}
        data['semantic_metadata']['laws_structured'] = citations
        data['semantic_metadata']['laws_confidence'] = data['enriched']['laws_confidence']
        
        # Salva
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        
        return True, len(citations), num_old
    
    except Exception as e:
        print(f"❌ Error processing {json_path.name}: {str(e)[:100]}")
        return False, 0, 0


def _calculate_confidence(citations: List[Dict]) -> Dict:
    """Calcola confidence metrics"""
    if not citations:
        return {'overall': 0.0}
    
    unique_laws = len(set(c['law_id'] for c in citations))
    with_articles = sum(1 for c in citations if c.get('article_num'))
    with_comma = sum(1 for c in citations if c.get('comma_num'))
    with_context = sum(1 for c in citations if c.get('context'))
    
    count_score = min(1.0, unique_laws / 5.0)
    granularity_score = with_articles / len(citations) if citations else 0
    
    overall = 0.5 * count_score + 0.5 * granularity_score
    
    return {
        'overall': round(overall, 3),
        'total_citations': len(citations),
        'unique_laws': unique_laws,
        'with_articles': with_articles,
        'with_comma': with_comma,
        'with_context': with_context,
        'granularity': round(granularity_score, 3)
    }


def main():
    parser = argparse.ArgumentParser(
        description='Re-enrich laws in JSON files (V3 - keeps ALL variants)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reenrich_laws_v3_FINAL.py --all              # Process all files
  python reenrich_laws_v3_FINAL.py --inail            # Only INAIL
  python reenrich_laws_v3_FINAL.py --sample 5 --all   # Test on 5 files
  python reenrich_laws_v3_FINAL.py --dry-run --all    # Show what would happen
        """
    )
    parser.add_argument('--inail', action='store_true', help='Process INAIL files')
    parser.add_argument('--osha', action='store_true', help='Process OSHA files')
    parser.add_argument('--all', action='store_true', help='Process all files')
    parser.add_argument('--sample', type=int, default=0, help='Process only N files (for testing)')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without saving')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show details for each file')
    
    args = parser.parse_args()
    
    if not any([args.inail, args.osha, args.all]):
        parser.print_help()
        return
    
    print("\n" + "="*70)
    print("🔄 RE-ENRICH LAWS V3 - KEEP ALL VARIANTS")
    print("="*70)
    print(f"⏰ Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    if args.dry_run:
        print("⚠️  DRY RUN MODE - No files will be modified")
    print()
    
    extractor = EnhancedLawExtractorV3()
    
    files_to_process = []
    
    if args.all or args.inail:
        if INAIL_ENRICHED.exists():
            inail_files = list(INAIL_ENRICHED.glob('*.json'))
            files_to_process.extend(inail_files)
            print(f"📘 INAIL: {len(inail_files)} files found")
        else:
            print(f"⚠️  INAIL path not found: {INAIL_ENRICHED}")
    
    if args.all or args.osha:
        if OSHA_ENRICHED.exists():
            osha_files = list(OSHA_ENRICHED.glob('*.json'))
            files_to_process.extend(osha_files)
            print(f"📗 OSHA: {len(osha_files)} files found")
        else:
            print(f"⚠️  OSHA path not found: {OSHA_ENRICHED}")
    
    if args.sample > 0:
        files_to_process = files_to_process[:args.sample]
        print(f"📋 Sample mode: processing {len(files_to_process)} files")
    
    print(f"\n🔄 Processing {len(files_to_process)} files...\n")
    
    total_new = 0
    total_old = 0
    success_count = 0
    improved_count = 0
    
    for json_file in tqdm(files_to_process, desc="Re-enriching"):
        success, num_new, num_old = process_json_file(json_file, extractor, dry_run=args.dry_run)
        
        if success:
            success_count += 1
            total_new += num_new
            total_old += num_old
            
            if num_new > num_old:
                improved_count += 1
            
            if args.verbose:
                change = "↑" if num_new > num_old else ("=" if num_new == num_old else "↓")
                tqdm.write(f"  {json_file.name[:50]}: {num_old} → {num_new} {change}")
    
    # Summary
    print(f"\n{'='*70}")
    print(f"✅ RE-ENRICHMENT {'SIMULATION ' if args.dry_run else ''}COMPLETE!")
    print(f"{'='*70}")
    print(f"   Files processed: {success_count}/{len(files_to_process)}")
    print(f"   Total citations OLD: {total_old:,}")
    print(f"   Total citations NEW: {total_new:,}")
    print(f"   Change: {'+' if total_new > total_old else ''}{total_new - total_old:,} ({((total_new/total_old - 1) * 100) if total_old > 0 else 0:.1f}%)")
    print(f"   Files improved: {improved_count}")
    if success_count > 0:
        print(f"   Avg citations/file: {total_new/success_count:.1f}")
    print(f"{'='*70}\n")
    
    if args.dry_run:
        print("💡 Run without --dry-run to apply changes\n")


if __name__ == "__main__":
    main()
