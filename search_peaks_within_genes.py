import re
import sys
import contextlib

docr = sys.argv[1]
deg = sys.argv[2]

with open (docr) as file:
     seq = file.read().splitlines()
     
     
with open (deg) as f:
	degs = f.read().splitlines()     
     
    
#>>> seq[0]
#'>25:34595225-34596108\tCAATGCGCCATGAGCTCTGTAATGTCATTATATTTATGC
#re.split(r'>|:|-', seq[0].split()[0])


#chr_peak=re.split(r'>|:|-', seq[0].split()[0])[1]
#start_peak=re.split(r'>|:|-', seq[0].split()[0])[2]
#end_peak=re.split(r'>|:|-', seq[0].split()[0])[3]
#sequence_peak=re.split(r'>|:|-', seq[0].split()[1])



#>>> degs[0]
#'21 18445942 20484289 slc45a2'


#chr_aug_gene=degs[0].split()[0]
#start_aug_gene=degs[0].split()[1]
#end_aug_gene=degs[0].split()[2]
#name_aug_gene=degs[0].split()[3]


matches = []

for gene in degs:


	chr_aug_gene=gene.split()[0]
	start_aug_gene=gene.split()[1]
	end_aug_gene=gene.split()[2]
	name_aug_gene=gene.split()[3]


	for peak in seq:
		chr_peak=re.split(r'>|:|-', peak.split()[0])[1]
		start_peak=re.split(r'>|:|-', peak.split()[0])[2]
		end_peak=re.split(r'>|:|-', peak.split()[0])[3]
		sequence_peak=re.split(r'>|:|-', peak.split()[1])

		if ((chr_peak == chr_aug_gene) & (int(start_peak) >= int(start_aug_gene)) & (int(end_peak) <= int(end_aug_gene))):
			mat =  chr_peak + ":" + start_peak + "-" + end_peak + " peak within reach of gene " + name_aug_gene
		
			matches.append(mat) 


with open(sys.argv[3], 'w') as f:

	f.write("\n".join(map(str, matches)))




