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
			
			out_file.write(str(chrom)+" "+str(position)+"\t"+str(ID)+"\t"+str(gene)+"\t"+str(Ref)+"\t"+str(Alt)+"\t")

			paralog_codons = paralog_codons.split("|&")[1]
			paralog_codons = paralog_codons.split("&")
			paralog_codons.remove("")
			refid_flav2_check = 0
			if refid_flavour == "all_con":
				if all(refids.endswith("REFID=1") for refids in paralog_codons) == True:
					refid_flav2_check = 1

			for codon in paralog_codons:

				if refid_flavour == "noQC" or refid_flav2_check == 1:
					refID = codon.split(":")[-1]
					codon = codon.split(":")[1]
					paralog_chrom = codon.split("_")[0].lstrip("chr")
					paralog_position = int(codon.split("-")[1])

					out_file.write(
						str(paralog_chrom)+" "+str(paralog_position-2)+"\t"+
						str(paralog_chrom)+" "+str(paralog_position-1)+"\t"+
						str(paralog_chrom)+" "+str(paralog_position)+"\t"
						)

				if refid_flavour == "para_con":
					refID = codon.split(":")[-1]
					if refID == "REFID=1":
						codon = codon.split(":")[1]
						paralog_chrom = codon.split("_")[0].lstrip("chr")
						paralog_position = int(codon.split("-")[1])

						out_file.write(
							str(paralog_chrom)+" "+str(paralog_position-2)+"\t"+
							str(paralog_chrom)+" "+str(paralog_position-1)+"\t"+
							str(paralog_chrom)+" "+str(paralog_position)+"\t"
							)

			out_file.write("\n")

			if len(paralog_codons)*3 + 3 > max_split:
				max_split = len(paralog_codons)*3 + 3
				max_paralog_codons = paralog_codons
	# print("THE MAX NUMBER OF COLUMNS: ",max_split, max_paralog_codons)	#number of columns for dataframe in R; + 2 for (chrom poistion) and (ID)
out_file.close()

#Remove last "\n" from file so R doesn't think its a separate column when reading in
# infile = input_file+"2.noQC" #paralogs2 file

out_file_name = infile+"woLastCol"
out_file = open(infile+"woLastCol", "w")
with open(infile) as f:
	for line in f:
		line = line.rstrip()
		out_file.write(str(line)+"\n")

out_file.close()
os.system("mv " + out_file_name + " " + infile)

print("File_prep_for_R Done!")


