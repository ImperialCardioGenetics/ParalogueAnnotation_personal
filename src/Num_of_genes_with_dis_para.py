import os, sys, subprocess, re

# script to find Number of disease protein coding genes that have at least 1 paralogue and at least 1 P/LP variant aligning to a non-gap position of another disease gene

input_file = "clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_paralogs"

with open(input_file, "r", encoding="utf-8") as f:
	for line in f:
		if line != "\n" and not line.startswith("#"):
			chrom = line.split()[0]
			position = line.split()[1]
			ID = line.split()[2]
			Ref = line.split()[3]
			Alt = line.split()[4]
			
			line = line.rstrip()
			paralog_codons = line.split("PARALOGS->")[1]
			gene = paralog_codons.split("|")[3]
			
			# out_file.write(str(chrom)+" "+str(position)+"\t"+str(ID)+"\t"+str(gene)+"\t"+str(Ref)+"\t"+str(Alt)+"\t")

			paralog_codons = paralog_codons.split("|&")[1]
			paralog_codons = paralog_codons.split("&")
			paralog_codons.remove("")
			

			for codon in paralog_codons:

				refID = codon.split(":")[-1]
				paralog_gene = codon.split(":")[0]
				codon = codon.split(":")[1]
				paralog_chrom = codon.split("_")[0].lstrip("chr")
				paralog_position = int(codon.split("-")[1])

			# 	out_file.write(
			# 		str(paralog_chrom)+" "+str(paralog_position-2)+"\t"+
			# 		str(paralog_chrom)+" "+str(paralog_position-1)+"\t"+
			# 		str(paralog_chrom)+" "+str(paralog_position)+"\t"
			# 		)

			# out_file.write("\n")

			


