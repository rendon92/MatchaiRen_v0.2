##### READ QUERY DNA FASTA SEQUENCE FROM FILE AND CREATE N RANDOM SEQUENCES WITH SIMILAR LENGTH AND SIMILAR GC CONTENT

import re
import random
import contextlib
import sys



output_filepath=sys.argv[1]


# FUNCTION TO GET NUCLEOTIDE FREQUENCY FOR A GIVEN SEQUENCE

def base_frequency(dna):
    d = {}
    for base in 'ATCG':
        d[base] = dna.count(base)/float(len(dna))
    return d


with open (sys.argv[2]) as f:
	lines = f.read().splitlines()
	
	
freqs = base_frequency(lines[1])

### freqs['A']
### freqs['C']
### freqs['G']
### freqs['T']


nucleotide = ['A', 'C', 'T', 'G']



#### Here we create sequence with same nucleotide frequencies as the query and with the same nucleotide length



with open(output_filepath, "w") as o:

	for i in range(10):

		a = random.choices(nucleotide, weights = (freqs['A'], freqs['C'], freqs['G'], freqs['T']), k=len(lines[1]))
		b = ''.join(a)
		
		with contextlib.redirect_stdout(o):
	
			print('>', i, '\n', b, sep='')
	

#run homer
#homer2 known -i fasta_test.tsv -b background_test.fa -m /data/GRN_software/homer/data/knownTFs/vertebrates/all.motifs >  test


