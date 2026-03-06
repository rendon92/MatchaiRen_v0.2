#!/usr/bin/env python3
# convertir_masivo_v2.py

import sys

def convertir_masivo(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        lines = f_in.readlines()
        i = 0
        motivos = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('>'):
                # Quitar el '>'
                resto = line[1:]
                
                # Encontrar el último '::' para separar ID compuesto del resto
                if '::' in resto:
                    # Dividir por el último '::'
                    ultimo_doble = resto.rfind('::')
                    id_compuesto = resto[:ultimo_doble]  # Todo hasta el último '::'
                    resto2 = resto[ultimo_doble+2:]      # Lo que sigue después del último '::'
                    
                    # resto2 ahora tiene: "NOMBRE CONSENSO SCORE"
                    partes = resto2.split()
                    if len(partes) >= 3:
                        tf_name = partes[0]
                        consensus = partes[1]
                        score = partes[2]
                    else:
                        # Si no hay 3 partes, asignar valores por defecto
                        tf_name = partes[0] if len(partes) > 0 else "UNKNOWN"
                        consensus = partes[1] if len(partes) > 1 else "NNNNNN"
                        score = partes[2] if len(partes) > 2 else "5.64"
                    
                    # Escribir header con TABS
                    f_out.write(f">{consensus}\t{id_compuesto}::{tf_name}\t{score}\n")
                    
                    # Leer matriz (las siguientes líneas hasta próximo '>')
                    i += 1
                    while i < len(lines) and not lines[i].startswith('>'):
                        mat_line = lines[i].strip()
                        if mat_line and not mat_line.startswith('#'):
                            # Convertir espacios a TABS
                            numeros = mat_line.split()
                            if len(numeros) == 4:
                                f_out.write(f"{numeros[0]}\t{numeros[1]}\t{numeros[2]}\t{numeros[3]}\n")
                        i += 1
                    
                    # Línea en blanco entre motivos
                    f_out.write('\n')
                    motivos += 1
                else:
                    i += 1
            else:
                i += 1
        
        print(f"Motivos convertidos: {motivos}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python convertir_masivo_v2.py <input.motifs> <output.motifs>")
        sys.exit(1)
    
    convertir_masivo(sys.argv[1], sys.argv[2])
