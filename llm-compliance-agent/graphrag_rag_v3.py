#!/usr/bin/env python3

"""

🎯 GRAPHRAG SAFETY SYSTEM v2 - Complete RAG with Improved Prompts



USA: hybrid_retriever_v5.py (ORIGINALE)



Integrates:

- Hybrid Retriever v5 (Qdrant + Neo4j)

- Qwen 2.5 72B via Ollama

- Prompt engineering ottimizzato per risposte naturali



Usage:

1. Tunnel: ssh -N -L 16333:localhost:6333 -L 17687:localhost:7687 dgx003 &

2. Ollama: ./ollama serve > /tmp/ollama.log 2>&1 &

3. Run: python graphrag_rag_v2.py -i

"""



import requests

import json

import sys

import os

from typing import List, Dict, Optional



# USA IL RETRIEVER v5 FIXED (estrae leggi dalla query!)

from hybrid_retriever_v5_fixed import HybridRetrieverV5



# ============================================

# CONFIGURATION

# ============================================



OLLAMA_HOST = "http://127.0.0.1:11434"

MODEL_NAME = "qwen2.5:72b"



# ============================================

# SYSTEM PROMPTS MIGLIORATI

# ============================================



SYSTEM_PROMPT_IT = """Sei un esperto consulente di sicurezza sul lavoro italiano con 20 anni di esperienza.



REGOLE FONDAMENTALI:

1. Rispondi SEMPRE in italiano, in modo chiaro e professionale

2. Basa le risposte ESCLUSIVAMENTE sui documenti forniti nel contesto

3. Quando citi una norma, specifica SEMPRE se è VIGENTE, MODIFICATA o ABROGATA

4. Usa un tono professionale ma accessibile, evita burocratese inutile

5. Struttura le risposte in paragrafi brevi e leggibili

6. NON inventare informazioni - se non trovi qualcosa nei documenti, dillo



FORMATO RISPOSTA:

- Inizia con una risposta DIRETTA alla domanda

- Poi aggiungi dettagli e riferimenti normativi

- Concludi con eventuali note pratiche

- Evita elenchi puntati infiniti - preferisci prosa fluida"""



SYSTEM_PROMPT_LEGAL = """Sei un giurista specializzato in diritto del lavoro e sicurezza occupazionale.



REGOLE:

1. Rispondi in italiano tecnico-giuridico ma comprensibile

2. Cita SEMPRE gli articoli specifici (es. "art. 37 del D.Lgs. 81/2008")

3. Segnala CHIARAMENTE se una norma è stata ABROGATA o MODIFICATA

4. Se ci sono più interpretazioni, menzionale brevemente



⚠️ Le norme cambiano! Se il sistema indica che una legge è ABROGATA, comunicalo chiaramente."""



SYSTEM_PROMPT_PRACTICAL = """Sei un tecnico della prevenzione con esperienza operativa.



REGOLE:

1. Rispondi in modo PRATICO e OPERATIVO

2. Usa esempi concreti quando possibile  

3. Indica le figure professionali coinvolte (RSPP, RLS, MC, ecc.)

4. Specifica i documenti necessari (DVR, POS, DUVRI, ecc.)

5. Dai priorità all'aspetto applicativo"""





# ============================================

# OLLAMA CLIENT

# ============================================



class OllamaClient:

    """Ollama client via REST API"""

    

    def __init__(self, host: str = OLLAMA_HOST, model: str = MODEL_NAME):

        self.host = host

        self.model = model

        self._test_connection()

    

    def _test_connection(self):

        try:

            response = requests.get(f"{self.host}/api/tags", timeout=10)

            if response.status_code == 200:

                models = response.json().get('models', [])

                model_names = [m['name'] for m in models]

                

                if any(self.model in name for name in model_names):

                    print(f"   ✅ Ollama OK - {self.model} disponibile")

                else:

                    print(f"   ⚠️ Model {self.model} non trovato, provo comunque...")

            else:

                raise Exception(f"Ollama returned {response.status_code}")

        except requests.exceptions.ConnectionError:

            print(f"   ❌ Ollama non raggiungibile a {self.host}")

            print(f"      Avvia con: ./ollama serve &")

            raise

    

    def generate(

        self, 

        prompt: str, 

        system: str = SYSTEM_PROMPT_IT,

        temperature: float = 0.2,

        max_tokens: int = 1500,

        stream: bool = True

    ) -> str:

        

        url = f"{self.host}/api/generate"

        

        payload = {

            "model": self.model,

            "prompt": prompt,

            "system": system,

            "stream": stream,

            "options": {

                "temperature": temperature,

                "num_predict": max_tokens,

                "top_p": 0.9,

                "repeat_penalty": 1.1,

            }

        }

        

        try:

            if stream:

                response = requests.post(url, json=payload, stream=True, timeout=180)

                

                full_response = ""

                for line in response.iter_lines():

                    if line:

                        try:

                            data = json.loads(line)

                            chunk = data.get('response', '')

                            full_response += chunk

                            print(chunk, end='', flush=True)

                            

                            if data.get('done', False):

                                break

                        except json.JSONDecodeError:

                            continue

                

                print()

                return full_response

            else:

                response = requests.post(url, json=payload, timeout=180)

                if response.status_code == 200:

                    return response.json().get('response', '')

                else:

                    raise Exception(f"Ollama error: {response.status_code}")

        except requests.exceptions.Timeout:

            return "⚠️ Timeout nella generazione. Riprova con domanda più breve."





# ============================================

# CONTEXT FORMATTER

# ============================================



def format_context_v2(results: List[Dict], max_docs: int = 5) -> str:

    """Formatta contesto per LLM - versione migliorata"""

    

    if not results:

        return "Nessun documento trovato."

    

    context_parts = []

    all_warnings = []

    

    for i, result in enumerate(results[:max_docs], 1):

        payload = result['payload']

        

        # Info base

        title = payload.get('titolo') or payload.get('title') or payload.get('doc_id', 'Documento')

        source = payload.get('source', 'N/A')

        category = payload.get('categoria_principale', '')

        text = payload.get('text', payload.get('article_text', ''))[:2000]

        sintesi = payload.get('sintesi', '')[:400]

        

        # Header documento

        doc_section = f"""

══════════════════════════════════════════════════════════════════

📄 DOCUMENTO {i}: {title[:70]}

══════════════════════════════════════════════════════════════════

Fonte: {source} | Categoria: {category}

Score: {result.get('final_score', result.get('score', 0)):.3f}

"""

        

        # Sintesi

        if sintesi:

            doc_section += f"\n📋 Sintesi: {sintesi}\n"

        

        # Contenuto

        doc_section += f"\n📝 Contenuto:\n{text}\n"

        

        # WARNINGS dal retriever v5 (questi funzionano!)

        warnings = result.get('warnings', [])

        if warnings:

            doc_section += f"\n⚠️ ATTENZIONE NORMATIVA:\n"

            for w in warnings[:5]:

                doc_section += f"   {w}\n"

                all_warnings.append(w)

        

        # Riferimenti normativi con stato

        law_refs = []

        for law_id, ctx in result.get('graph_context', {}).items():

            if ctx.get('is_abrogated'):

                status = "⛔ ABROGATA"

                abr_by = ctx.get('abrogated_by', [])

                if abr_by:

                    status += f" (da {', '.join(abr_by[:2])})"

            elif ctx.get('is_modified'):

                status = "📝 MODIFICATA"

                mod_by = ctx.get('modified_by', [])

                if mod_by:

                    status += f" (da {', '.join(mod_by[:2])})"

            elif ctx.get('is_vigente'):

                status = "✅ VIGENTE"

            else:

                status = ""

            

            if status:

                law_refs.append(f"  • {law_id} {status}")

        

        if law_refs:

            doc_section += f"\n⚖️ Stato normativo:\n" + "\n".join(law_refs[:6])

        

        context_parts.append(doc_section)

    

    # Sezione warnings globale se ci sono

    warning_section = ""

    if all_warnings:

        warning_section = f"""



╔══════════════════════════════════════════════════════════════════╗

║  ⚠️  AVVISI IMPORTANTI SULLO STATO DELLE NORME                    ║

╠══════════════════════════════════════════════════════════════════╣

"""

        for w in set(all_warnings):  # unique

            warning_section += f"║  {w[:65]}\n"

        warning_section += "╚══════════════════════════════════════════════════════════════════╝"

    

    return warning_section + "\n".join(context_parts)





def detect_query_type(question: str) -> str:

    """Rileva tipo domanda per prompt appropriato"""

    

    q = question.lower()

    

    legal_kw = ['articolo', 'comma', 'legge', 'decreto', 'dlgs', 'normativa', 

                'obbligo', 'sanzione', 'abrogat', 'vigente']

    practical_kw = ['come', 'procedura', 'cosa fare', 'devo', 'chi deve',

                    'documento', 'dvr', 'pos', 'cantiere', 'dpi']

    

    legal_score = sum(1 for kw in legal_kw if kw in q)

    practical_score = sum(1 for kw in practical_kw if kw in q)

    

    if legal_score > practical_score + 1:

        return 'legal'

    elif practical_score > legal_score:

        return 'practical'

    return 'general'





# ============================================

# RAG SYSTEM v2

# ============================================



class GraphRAGSystemV2:

    """GraphRAG System v2 - usa hybrid_retriever_v5"""

    

    def __init__(

        self,

        qdrant_port: int = 16333,

        neo4j_port: int = 17687,

        ollama_host: str = OLLAMA_HOST,

        model_name: str = MODEL_NAME

    ):

        print("\n" + "="*70)

        print("🚀 GRAPHRAG SAFETY SYSTEM v2")

        print("="*70 + "\n")

        

        # Set env vars per retriever

        os.environ['QDRANT_PORT'] = str(qdrant_port)

        os.environ['NEO4J_PORT'] = str(neo4j_port)

        

        # Retriever v5 ORIGINALE

        print("📚 Caricamento Retriever v5...")

        self.retriever = HybridRetrieverV5()

        

        # LLM

        print("\n🤖 Connessione Ollama...")

        print(f"   Model: {model_name}")

        self.llm = OllamaClient(host=ollama_host, model=model_name)

        

        print("\n" + "="*70)

        print("✅ SISTEMA PRONTO!")

        print("="*70 + "\n")

    

    def query(

        self, 

        question: str, 

        top_k: int = 5,

        temperature: float = 0.2,

        mode: str = 'auto',

        show_sources: bool = True,

        show_context: bool = False

    ) -> Dict:

        """Query principale"""

        

        print("\n" + "─"*70)

        print(f"❓ {question}")

        print("─"*70)

        

        # Step 1: Retrieve con v5

        print("\n🔍 Ricerca...")

        results = self.retriever.retrieve(question, top_k=top_k, include_graph_context=True)

        

        if not results:

            return {

                'answer': "Non ho trovato documenti rilevanti. Prova a riformulare.",

                'sources': [],

                'status': 'no_results'

            }

        

        # Conta warnings

        total_warnings = sum(len(r.get('warnings', [])) for r in results)

        print(f"   ✅ {len(results)} documenti, {total_warnings} warnings normativi")

        

        # Step 2: Format context

        context = format_context_v2(results, max_docs=top_k)

        

        if show_context:

            print("\n📋 CONTESTO:")

            print(context[:3000])

        

        # Step 3: Select prompt

        if mode == 'auto':

            query_type = detect_query_type(question)

        else:

            query_type = mode

        

        prompts = {

            'legal': SYSTEM_PROMPT_LEGAL,

            'practical': SYSTEM_PROMPT_PRACTICAL,

            'general': SYSTEM_PROMPT_IT

        }

        system_prompt = prompts.get(query_type, SYSTEM_PROMPT_IT)

        print(f"   📝 Modalità: {query_type}")

        

        # Step 4: Build prompt

        prompt = f"""Basandoti ESCLUSIVAMENTE sui seguenti documenti sulla sicurezza sul lavoro, rispondi alla domanda.



IMPORTANTE: Se i documenti indicano che una legge è ABROGATA o MODIFICATA, DEVI segnalarlo chiaramente nella risposta!



{context}



══════════════════════════════════════════════════════════════════

DOMANDA: {question}

══════════════════════════════════════════════════════════════════



Rispondi in italiano. Se citi norme, indica sempre se sono vigenti o abrogate.

Se non trovi informazioni sufficienti, dillo.



RISPOSTA:"""

        

        # Step 5: Generate

        print("\n🤖 Generazione risposta...\n")

        print("─"*70)

        

        answer = self.llm.generate(

            prompt=prompt,

            system=system_prompt,

            temperature=temperature,

            max_tokens=1500

        )

        

        print("─"*70)

        

        # Step 6: Sources

        if show_sources:

            print("\n📚 FONTI:")

            for i, result in enumerate(results[:top_k], 1):

                title = result['payload'].get('titolo') or result['payload'].get('title', 'N/A')

                source = result['payload'].get('source', '')

                score = result.get('final_score', 0)

                n_warn = len(result.get('warnings', []))

                warn_str = f" ⚠️{n_warn}" if n_warn > 0 else ""

                

                print(f"   {i}. [{source}] {str(title)[:50]}... ({score:.3f}){warn_str}")

        

        return {

            'answer': answer,

            'sources': results[:top_k],

            'question': question,

            'query_type': query_type,

            'warnings_count': total_warnings,

            'status': 'success'

        }

    

    def interactive_mode(self):

        """Modalità interattiva"""

        

        print("\n" + "="*70)

        print("💬 MODALITÀ INTERATTIVA")

        print("="*70)

        print("\nComandi: 'quit', 'help', 'mode X', 'debug'")

        

        current_mode = 'auto'

        show_debug = False

        

        while True:

            try:

                question = input("\n🔹 Tu: ").strip()

                

                if not question:

                    continue

                

                if question.lower() in ['quit', 'q', 'exit']:

                    print("\n👋 Arrivederci!")

                    break

                

                if question.lower() == 'help':

                    self._show_help()

                    continue

                

                if question.lower() == 'debug':

                    show_debug = not show_debug

                    print(f"   Debug: {'ON' if show_debug else 'OFF'}")

                    continue

                

                if question.lower().startswith('mode '):

                    m = question[5:].strip().lower()

                    if m in ['legal', 'practical', 'general', 'auto']:

                        current_mode = m

                        print(f"   ✅ Modalità: {current_mode}")

                    continue

                

                self.query(question, top_k=5, mode=current_mode, show_context=show_debug)

                

            except KeyboardInterrupt:

                print("\n\n👋 Bye!")

                break

            except Exception as e:

                print(f"\n❌ Errore: {e}")

    

    def _show_help(self):

        print("""

╔══════════════════════════════════════════════════════════════════╗

║                    ESEMPI DI DOMANDE                              ║

╠══════════════════════════════════════════════════════════════════╣

║  AMIANTO:                                                         ║

║    • Quali sono le procedure per la bonifica dell'amianto?        ║

║    • Chi redige il piano di lavoro per rimozione amianto?         ║

║                                                                   ║

║  ANTINCENDIO:                                                     ║

║    • Cosa prevede il minicodice prevenzione incendi?              ║

║    • Come si valuta il rischio incendio nei luoghi di lavoro?     ║

║                                                                   ║

║  D.LGS. 81/2008:                                                  ║

║    • Quali sono gli obblighi del datore di lavoro?                ║

║    • L'art. 37 cosa prevede sulla formazione?                     ║

║    • Il D.Lgs. 626/1994 è ancora in vigore?                       ║

║                                                                   ║

║  CANTIERI:                                                        ║

║    • Quando serve il coordinatore per la sicurezza?               ║

║    • Differenze tra POS e PSC?                                    ║

╚══════════════════════════════════════════════════════════════════╝

""")

    

    def batch_test(self, queries: List[str], output_file: str = None):

        """Test batch"""

        

        print(f"\n🧪 BATCH TEST - {len(queries)} query\n")

        

        results = []

        for i, q in enumerate(queries, 1):

            print(f"\n{'='*70}")

            print(f"TEST {i}/{len(queries)}")

            

            result = self.query(q, top_k=5)

            results.append({

                'query': q,

                'answer': result.get('answer', '')[:500],

                'warnings': result.get('warnings_count', 0),

                'status': result.get('status', 'unknown')

            })

        

        # Summary

        print(f"\n{'='*70}")

        print("📊 RIEPILOGO")

        print(f"{'='*70}")

        

        ok = sum(1 for r in results if r['status'] == 'success')

        total_warn = sum(r['warnings'] for r in results)

        print(f"   ✅ Successi: {ok}/{len(queries)}")

        print(f"   ⚠️ Warnings totali: {total_warn}")

        

        if output_file:

            with open(output_file, 'w', encoding='utf-8') as f:

                json.dump(results, f, ensure_ascii=False, indent=2)

            print(f"   📁 Salvato: {output_file}")

        

        return results

    

    def close(self):

        self.retriever.close()

        print("\n✅ Sistema chiuso")





# ============================================

# TEST QUERIES

# ============================================



TEST_QUERIES = [

    "Quali sono le procedure per la bonifica dell'amianto?",

    "Il D.Lgs. 626/1994 è ancora in vigore?",

    "Cosa prevede il minicodice prevenzione incendi?",

    "Obblighi del datore di lavoro secondo D.Lgs. 81/2008",

    "Chi nomina il RSPP?",

    "Valori limite esposizione agenti cancerogeni",

    "Quando serve il coordinatore sicurezza in cantiere?",

]





# ============================================

# MAIN

# ============================================



def main():

    import argparse

    

    parser = argparse.ArgumentParser(description='GraphRAG Safety v2')

    parser.add_argument('--query', '-q', type=str, help='Query singola')

    parser.add_argument('--interactive', '-i', action='store_true', help='Interattivo')

    parser.add_argument('--test', '-t', action='store_true', help='Batch test')

    parser.add_argument('--top-k', '-k', type=int, default=5)

    parser.add_argument('--temperature', type=float, default=0.2)

    parser.add_argument('--mode', choices=['auto', 'legal', 'practical', 'general'], default='auto')

    parser.add_argument('--qdrant-port', type=int, default=16333)

    parser.add_argument('--neo4j-port', type=int, default=17687)

    parser.add_argument('--ollama-host', type=str, default=OLLAMA_HOST)

    parser.add_argument('--model', type=str, default=MODEL_NAME)

    parser.add_argument('--debug', '-d', action='store_true')

    

    args = parser.parse_args()

    

    try:

        rag = GraphRAGSystemV2(

            qdrant_port=args.qdrant_port,

            neo4j_port=args.neo4j_port,

            ollama_host=args.ollama_host,

            model_name=args.model

        )

    except Exception as e:

        print(f"\n❌ Errore: {e}")

        print("\n💡 Verifica:")

        print("   1. Tunnel attivo: ssh -N -L 16333:localhost:6333 -L 17687:localhost:7687 dgx003 &")

        print("   2. Ollama: ./ollama serve &")

        sys.exit(1)

    

    try:

        if args.query:

            rag.query(args.query, top_k=args.top_k, temperature=args.temperature,

                     mode=args.mode, show_context=args.debug)

        elif args.interactive:

            rag.interactive_mode()

        elif args.test:

            rag.batch_test(TEST_QUERIES, output_file="test_results.json")

        else:

            # Demo

            print("\n🎬 DEMO\n")

            rag.query("Quali sono le procedure per la bonifica dell'amianto?", top_k=5)

            

            if input("\n⏸️ ENTER per continuare, q per uscire: ").lower() != 'q':

                rag.query("Il D.Lgs. 626/1994 è ancora in vigore?", top_k=5)

    

    finally:

        rag.close()





if __name__ == "__main__":

    main()
