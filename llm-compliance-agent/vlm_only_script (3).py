#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""

Script per applicare VLM a immagini già scaricate

Processa solo le immagini che NON hanno già un file .txt di descrizione

"""



import argparse

import os

import sys

from pathlib import Path

from datetime import datetime

from typing import List, Dict, Any

import json



# Configurazione cache HuggingFace

if 'HF_HOME' not in os.environ:

    hf_cache = Path.home() / '.cache' / 'huggingface'

    os.environ['HF_HOME'] = str(hf_cache)

    os.environ['TRANSFORMERS_CACHE'] = str(hf_cache)



# ============================================================================

# CONFIGURAZIONE VLM

# ============================================================================



class VLMConfig:

    """Configurazione VLM ottimizzata per HPC"""

    

    VLM_MODEL = "HuggingFaceTB/SmolVLM-256M-Instruct"

    VLM_USE_BFLOAT16 = True

    VLM_USE_FLASH_ATTENTION = True

    VLM_USE_8BIT_QUANTIZATION = False  # Cambia a True per risparmiare RAM

    VLM_IMAGE_SIZE = 512

    VLM_MAX_NEW_TOKENS = 300

    VLM_BATCH_SIZE = 8  # Numero immagini processate insieme



    VLM_PROMPT = """Describe this technical document image in detail for academic research. Focus on:

- Main content type (diagram, chart, photo, schematic)

- Key visual elements and their relationships

- Any visible text, labels, or numerical data

- Technical details relevant to workplace safety or industrial processes

Be precise and comprehensive."""



# ============================================================================

# TROVA IMMAGINI

# ============================================================================



def find_all_images(images_dir: Path) -> List[Dict[str, Any]]:

    """

    Trova tutte le immagini PNG nelle sottodirectory

    Returns: Lista di dict con {path, doc_id, filename, has_description}

    """

    

    if not images_dir.exists():

        print(f"[ERROR] Directory immagini non esiste: {images_dir}")

        return []

    

    images = []

    

    # Cerca in tutte le sottodirectory (una per documento)

    for doc_dir in images_dir.iterdir():

        if not doc_dir.is_dir():

            continue

            

        doc_id = doc_dir.name

        

        # Trova tutte le immagini PNG

        for img_path in doc_dir.glob("*.png"):

            txt_path = img_path.with_suffix('.txt')

            

            images.append({

                'path': img_path,

                'txt_path': txt_path,

                'doc_id': doc_id,

                'filename': img_path.name,

                'has_description': txt_path.exists()

            })

    

    return images



def filter_images_without_description(images: List[Dict]) -> List[Dict]:

    """Filtra solo immagini senza descrizione VLM"""

    return [img for img in images if not img['has_description']]



# ============================================================================

# APPLICAZIONE VLM

# ============================================================================



def apply_vlm_to_images(images: List[Dict], config: VLMConfig) -> Dict[str, Any]:

    """

    Applica VLM a lista di immagini con batch processing

    Returns: Statistiche del processing

    """

    

    if not images:

        print("[INFO] Nessuna immagine da processare")

        return {'total': 0, 'success': 0, 'failed': 0}

    

    try:

        from transformers import AutoProcessor, AutoModelForVision2Seq, BitsAndBytesConfig

        from PIL import Image

        import torch

        

        print(f"\n{'='*80}")

        print(f"INIZIALIZZAZIONE VLM")

        print(f"{'='*80}")

        print(f"Modello: {config.VLM_MODEL}")

        print(f"Immagini da processare: {len(images)}")

        print(f"{'='*80}\n")



        device = "cuda" if torch.cuda.is_available() else "cpu"

        print(f"[VLM] Device: {device}")

        

        if device == "cpu":

            print("[WARNING] GPU non rilevata! VLM sarà MOLTO lento su CPU")

            response = input("Continuare comunque? (y/n): ").strip().lower()

            if response != 'y':

                print("Operazione annullata")

                return {'total': 0, 'success': 0, 'failed': 0}



        # Inizializza processor

        print(f"[VLM] Caricamento processor...")

        processor = AutoProcessor.from_pretrained(

            config.VLM_MODEL,

            trust_remote_code=True,

            size={"longest_edge": config.VLM_IMAGE_SIZE}

        )



        # Configura modello

        print(f"[VLM] Caricamento modello...")

        model_kwargs = {

            "trust_remote_code": True,

            "device_map": "auto"

        }



        if config.VLM_USE_8BIT_QUANTIZATION:

            print("[VLM] → Quantizzazione 8-bit (risparmio RAM)")

            quantization_config = BitsAndBytesConfig(load_in_8bit=True)

            model_kwargs["quantization_config"] = quantization_config

        else:

            if config.VLM_USE_BFLOAT16 and torch.cuda.is_available():

                print("[VLM] → bfloat16 precision")

                model_kwargs["torch_dtype"] = torch.bfloat16

            else:

                model_kwargs["torch_dtype"] = torch.float16 if torch.cuda.is_available() else torch.float32



        if config.VLM_USE_FLASH_ATTENTION and device == "cuda":

            print("[VLM] → Flash Attention 2")

            model_kwargs["_attn_implementation"] = "flash_attention_2"

        else:

            model_kwargs["_attn_implementation"] = "eager"



        model = AutoModelForVision2Seq.from_pretrained(

            config.VLM_MODEL,

            **model_kwargs

        )



        print(f"[VLM] ✓ Modello caricato")

        print(f"[VLM] Risoluzione: {config.VLM_IMAGE_SIZE}x{config.VLM_IMAGE_SIZE}")

        print(f"[VLM] Batch size: {config.VLM_BATCH_SIZE}")

        print()



        # Statistiche

        stats = {

            'total': len(images),

            'success': 0,

            'failed': 0,

            'skipped': 0

        }



        # BATCH PROCESSING

        total_images = len(images)

        batch_size = config.VLM_BATCH_SIZE



        for batch_start in range(0, total_images, batch_size):

            batch_end = min(batch_start + batch_size, total_images)

            batch = images[batch_start:batch_end]



            batch_num = batch_start // batch_size + 1

            total_batches = (total_images + batch_size - 1) // batch_size

            

            print(f"{'='*80}")

            print(f"BATCH {batch_num}/{total_batches} ({batch_start+1}-{batch_end} di {total_images})")

            print(f"{'='*80}")



            # Carica immagini del batch

            pil_images = []

            valid_indices = []

            valid_image_data = []



            for idx, img_data in enumerate(batch):

                global_idx = batch_start + idx

                

                try:

                    # Carica immagine

                    image = Image.open(img_data['path']).convert('RGB')



                    # Ridimensiona se troppo grande

                    max_size = config.VLM_IMAGE_SIZE * 4

                    if max(image.size) > max_size:

                        ratio = max_size / max(image.size)

                        new_size = tuple(int(dim * ratio) for dim in image.size)

                        image = image.resize(new_size, Image.Resampling.LANCZOS)

                        print(f"  [{global_idx+1:4d}] Ridimensionata: {image.size[0]}x{image.size[1]} | {img_data['filename']}")

                    else:

                        print(f"  [{global_idx+1:4d}] OK: {image.size[0]}x{image.size[1]} | {img_data['filename']}")



                    pil_images.append(image)

                    valid_indices.append(global_idx)

                    valid_image_data.append(img_data)



                except Exception as e:

                    print(f"  [{global_idx+1:4d}] ✗ Errore caricamento: {e}")

                    stats['failed'] += 1

                    continue



            if not pil_images:

                print("  Nessuna immagine valida in questo batch, skip")

                continue



            # Process batch con VLM

            try:

                print(f"\n  Processing {len(pil_images)} immagini con VLM...")

                

                # Prepara messaggi

                messages_batch = []

                for _ in pil_images:

                    messages_batch.append([{

                        "role": "user",

                        "content": [

                            {"type": "image"},

                            {"type": "text", "text": config.VLM_PROMPT}

                        ]

                    }])



                # Applica chat template

                prompts = [

                    processor.apply_chat_template(msgs, add_generation_prompt=True)

                    for msgs in messages_batch

                ]



                # Process con VLM

                inputs = processor(

                    text=prompts,

                    images=pil_images,

                    return_tensors="pt",

                    padding=True

                )

                inputs = {k: v.to(device) for k, v in inputs.items()}



                # Genera descrizioni

                with torch.no_grad():

                    outputs = model.generate(

                        **inputs,

                        max_new_tokens=config.VLM_MAX_NEW_TOKENS,

                        do_sample=False

                    )



                generated_texts = processor.batch_decode(

                    outputs,

                    skip_special_tokens=True

                )



                # Salva risultati

                print(f"\n  Salvataggio descrizioni...")

                for i, (global_idx, text, img_data) in enumerate(zip(valid_indices, generated_texts, valid_image_data)):

                    

                    # Estrai solo la risposta dell'assistente

                    if "Assistant:" in text:

                        description = text.split("Assistant:")[-1].strip()

                    else:

                        description = text.strip()



                    # Salva file .txt

                    try:

                        with open(img_data['txt_path'], 'w', encoding='utf-8') as f:

                            f.write("="*80 + "\n")

                            f.write(f"VLM DESCRIPTION - {img_data['filename']}\n")

                            f.write("="*80 + "\n\n")

                            f.write("IMAGE METADATA:\n")

                            f.write(f"  - Document ID: {img_data['doc_id']}\n")

                            f.write(f"  - Filename: {img_data['filename']}\n")

                            f.write(f"  - Image path: {img_data['path']}\n")

                            f.write(f"  - Model: {config.VLM_MODEL}\n")

                            f.write(f"  - Resolution: {config.VLM_IMAGE_SIZE}px\n")

                            f.write(f"  - Timestamp: {datetime.now().isoformat()}\n\n")

                            f.write("VLM ANALYSIS:\n")

                            f.write("-"*80 + "\n")

                            f.write(description)

                            f.write("\n" + "-"*80 + "\n")



                        stats['success'] += 1

                        print(f"  [{global_idx+1:4d}] ✓ {description[:60]}...")

                        

                    except Exception as e:

                        print(f"  [{global_idx+1:4d}] ✗ Errore salvataggio: {e}")

                        stats['failed'] += 1



            except Exception as e:

                print(f"\n  ✗ Errore batch processing: {e}")

                import traceback

                traceback.print_exc()

                stats['failed'] += len(valid_indices)



            print()  # Linea vuota tra batch



        # Cleanup

        print("[VLM] Cleanup...")

        del model, processor

        if torch.cuda.is_available():

            torch.cuda.empty_cache()



        print(f"\n{'='*80}")

        print(f"VLM COMPLETATO")

        print(f"{'='*80}")

        print(f"Totale immagini: {stats['total']}")

        print(f"✓ Successo: {stats['success']}")

        print(f"✗ Fallite: {stats['failed']}")

        print(f"{'='*80}\n")



        return stats



    except ImportError as e:

        print(f"[ERROR] Librerie VLM non installate: {e}")

        print("\nInstalla con:")

        print("  pip install transformers torch torchvision accelerate")

        return {'total': 0, 'success': 0, 'failed': 0}

    

    except Exception as e:

        print(f"[ERROR] VLM fallito: {e}")

        import traceback

        traceback.print_exc()

        return {'total': 0, 'success': 0, 'failed': 0}



# ============================================================================

# AGGIORNAMENTO JSON (opzionale)

# ============================================================================



def update_json_with_vlm_descriptions(json_dir: Path, images_dir: Path):

    """

    Aggiorna i file JSON esistenti aggiungendo le descrizioni VLM

    """

    

    print(f"\n{'='*80}")

    print("AGGIORNAMENTO FILE JSON")

    print(f"{'='*80}\n")

    

    if not json_dir.exists():

        print("[WARNING] Directory JSON non trovata, skip update")

        return

    

    json_files = list(json_dir.glob("*.json"))

    json_files = [f for f in json_files if not f.name.startswith('vector_db')]

    

    updated = 0

    skipped = 0

    

    for json_file in json_files:

        try:

            with open(json_file, 'r', encoding='utf-8') as f:

                data = json.load(f)

            

            # Cerca se questo documento ha immagini

            if 'document_content' not in data or not data['document_content']:

                skipped += 1

                continue

            

            figures = data['document_content'].get('figures', [])

            exported_images = data['document_content'].get('exported_images', [])

            

            if not exported_images:

                skipped += 1

                continue

            

            # Aggiorna descrizioni VLM

            changed = False

            for img_data in exported_images:

                img_path = Path(img_data['path'])

                txt_path = img_path.with_suffix('.txt')

                

                if txt_path.exists() and not img_data.get('vlm_description'):

                    # Leggi descrizione

                    try:

                        with open(txt_path, 'r', encoding='utf-8') as f:

                            content = f.read()

                            # Estrai solo la descrizione

                            if "VLM ANALYSIS:" in content:

                                description = content.split("VLM ANALYSIS:")[1].split("-"*80)[1].strip()

                                img_data['vlm_description'] = description

                                img_data['vlm_description_file'] = str(txt_path)

                                changed = True

                    except:

                        pass

            

            # Aggiorna anche figures

            for idx, fig_info in enumerate(figures):

                if idx < len(exported_images):

                    vlm_desc = exported_images[idx].get('vlm_description')

                    if vlm_desc and not fig_info.get('vlm_description'):

                        fig_info['vlm_description'] = vlm_desc

                        changed = True

            

            if changed:

                # Salva JSON aggiornato

                with open(json_file, 'w', encoding='utf-8') as f:

                    json.dump(data, f, ensure_ascii=False, indent=2)

                updated += 1

                print(f"  ✓ {json_file.name}")

            else:

                skipped += 1

        

        except Exception as e:

            print(f"  ✗ {json_file.name}: {e}")

            skipped += 1

    

    print(f"\n[INFO] JSON aggiornati: {updated}, skipped: {skipped}\n")



# ============================================================================

# MAIN

# ============================================================================



def parse_arguments():

    """Parser argomenti command-line"""

    parser = argparse.ArgumentParser(

        description='Applica VLM a immagini già scaricate',

        formatter_class=argparse.RawDescriptionHelpFormatter,

        epilog="""

Esempi:



  # Processa tutte le immagini senza descrizione

  python apply_vlm_to_existing_images.py --images-dir /path/to/data/images



  # Con batch size più grande (se hai molta RAM)

  python apply_vlm_to_existing_images.py --images-dir /path/to/data/images --batch-size 16



  # Con quantizzazione 8-bit (risparmio RAM)

  python apply_vlm_to_existing_images.py --images-dir /path/to/data/images --use-8bit



  # Forza reprocessing anche se descrizioni esistono

  python apply_vlm_to_existing_images.py --images-dir /path/to/data/images --force

        """

    )



    parser.add_argument(

        '--images-dir',

        type=str,

        required=True,

        help='Directory contenente le immagini (es: /path/to/data/images)'

    )



    parser.add_argument(

        '--json-dir',

        type=str,

        default=None,

        help='Directory JSON da aggiornare (opzionale, default: ../json rispetto a images-dir)'

    )



    parser.add_argument(

        '--batch-size',

        type=int,

        default=8,

        help='Numero immagini processate insieme (default: 8)'

    )



    parser.add_argument(

        '--image-size',

        type=int,

        default=512,

        choices=[256, 512, 1024, 1536, 2048],

        help='Risoluzione patches immagini (default: 512)'

    )



    parser.add_argument(

        '--use-8bit',

        action='store_true',

        help='Usa quantizzazione 8-bit per risparmiare RAM'

    )



    parser.add_argument(

        '--force',

        action='store_true',

        help='Riprocessa anche immagini che hanno già descrizioni'

    )



    parser.add_argument(

        '--update-json',

        action='store_true',

        default=True,

        help='Aggiorna file JSON con nuove descrizioni VLM (default: True)'

    )



    parser.add_argument(

        '--no-update-json',

        action='store_true',

        help='NON aggiornare i file JSON'

    )



    return parser.parse_args()



def main():

    """Main function"""

    args = parse_arguments()



    print(f"\n{'='*80}")

    print("VLM APPLICATION TO EXISTING IMAGES")

    print(f"{'='*80}\n")



    # Verifica directory

    images_dir = Path(args.images_dir)

    if not images_dir.exists():

        print(f"[ERROR] Directory immagini non esiste: {images_dir}")

        sys.exit(1)



    # Determina json_dir

    if args.json_dir:

        json_dir = Path(args.json_dir)

    else:

        json_dir = images_dir.parent / 'json'



    # Configura VLM

    config = VLMConfig()

    config.VLM_BATCH_SIZE = args.batch_size

    config.VLM_IMAGE_SIZE = args.image_size

    config.VLM_USE_8BIT_QUANTIZATION = args.use_8bit



    print(f"Images directory: {images_dir}")

    print(f"JSON directory: {json_dir}")

    print(f"Batch size: {config.VLM_BATCH_SIZE}")

    print(f"Image resolution: {config.VLM_IMAGE_SIZE}px")

    print(f"8-bit quantization: {config.VLM_USE_8BIT_QUANTIZATION}")

    print(f"Force reprocess: {args.force}")

    print()



    # Trova immagini

    print("Scansione directory immagini...")

    all_images = find_all_images(images_dir)

    

    if not all_images:

        print("[ERROR] Nessuna immagine trovata!")

        sys.exit(1)



    print(f"[INFO] Trovate {len(all_images)} immagini totali")



    # Filtra quelle senza descrizione

    if args.force:

        images_to_process = all_images

        print(f"[INFO] FORCE mode: processerò tutte le {len(images_to_process)} immagini")

    else:

        images_to_process = filter_images_without_description(all_images)

        already_done = len(all_images) - len(images_to_process)

        print(f"[INFO] Immagini già con descrizione: {already_done}")

        print(f"[INFO] Immagini da processare: {len(images_to_process)}")



    if not images_to_process:

        print("\n[INFO] Tutte le immagini hanno già una descrizione VLM!")

        print("[INFO] Usa --force per riprocessare comunque")

        sys.exit(0)



    # Conferma

    print(f"\n{'='*80}")

    print("RIEPILOGO")

    print(f"{'='*80}")

    print(f"Immagini da processare: {len(images_to_process)}")

    print(f"Batch size: {config.VLM_BATCH_SIZE}")

    print(f"Stima tempo (GPU): ~{len(images_to_process) * 3 / 60:.1f} minuti")

    print(f"Stima tempo (CPU): ~{len(images_to_process) * 30 / 60:.1f} minuti")

    print(f"{'='*80}\n")



    response = input("Procedere? (y/n): ").strip().lower()

    if response != 'y':

        print("Operazione annullata")

        sys.exit(0)



    # Applica VLM

    stats = apply_vlm_to_images(images_to_process, config)



    # Aggiorna JSON se richiesto

    if args.no_update_json:

        print("[INFO] Skip aggiornamento JSON (--no-update-json)")

    elif args.update_json and stats['success'] > 0:

        update_json_with_vlm_descriptions(json_dir, images_dir)



    # Report finale

    print(f"\n{'='*80}")

    print("OPERAZIONE COMPLETATA")

    print(f"{'='*80}")

    print(f"Immagini processate: {stats['success']}/{stats['total']}")

    print(f"Fallite: {stats['failed']}")

    if stats['success'] > 0:

        print(f"\nDescrizioni salvate in: {images_dir}/*/figure_*.txt")

    print(f"{'='*80}\n")



    sys.exit(0 if stats['failed'] == 0 else 1)



if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print("\n\n[INTERRUPTED] Operazione interrotta dall'utente")

        sys.exit(130)

    except Exception as e:

        print(f"\n[FATAL ERROR] {e}")

        import traceback

        traceback.print_exc()

        sys.exit(1)
