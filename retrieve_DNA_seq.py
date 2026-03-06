## Retrieve DNA sequence from genomic coordinates

import sys
import re
import pyfaidx
from pyfaidx import Fasta

peakfile = sys.argv[1]
out_file = sys.argv[2]
genome = sys.argv[3]

#genes = Fasta('/data/genomes/Homo_sapiens/fasta/Homo_sapiens.GRCh38.dna.primary_assembly_renamed.fa')
genes = Fasta(genome)

peaks = []
with open (peakfile) as file:
	for line in file:
		line = line.strip()
		peaks.append(line)
		

# Verificación: archivo de picos vacío
if not peaks:
    sys.exit("Error: el archivo de picos está vacío o mal formateado.")


with open(out_file, "w") as o:
	for x in peaks:
		seqFile = genes[(re.split(':|-',x)[0])][int(re.split(':|-',x)[1]):int(re.split(':|-',x)[2])]
		o.write('>' + seqFile.fancy_name + '\n')
		o.write(seqFile.seq + '\n')
