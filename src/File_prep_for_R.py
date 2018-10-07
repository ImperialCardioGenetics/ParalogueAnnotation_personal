import os, sys, subprocess, re

#SCRIPT FOR SEPARATING EACH LINE IN VCF OUTPUT FROM VEP_PARALOGANNO.PY INTO TAB SEP FILE READY FOR JOIN IN R
#don't use if VEP_ParalogAnno.py has been edited to already do this, check to see if same code exists in that script
###BELOW ONLY WORKS FOR CLINVAR VCF FILES (EXAC HAS NO INFO IN ID VARIABLE)

def R_file_prep(input_file, refid_flavour="noQC"):
	#input_file = sys.argv[1]	#path of VEP_ParalogAnno.py paralogs output file
	dir1 = input_file.rsplit("/", 1)[0]
	#refid_flavour = sys.argv[3] #how to consider the refID: "noQC", "para_con" or "all_con". "noQC" - don't do any quality control and just take all alignments as annotations. "para_con" - only consider paralogous pairwise alignments where refID = 1. "all_con" - only consider total gene family paralogous alignments where all refID in that family = 1

	if refid_flavour == "noQC":
		out_file = open(input_file+"2.noQC", "w")
		infile = input_file+"2.noQC"
	elif refid_flavour == "para_con":
		out_file = open(input_file+"2.para_con", "w")
		infile = input_file+"2.para_con"
	elif refid_flavour == "all_con":
		out_file = open(input_file+"2.all_con", "w")
		infile = input_file+"2.all_con"

	max_split = 0	#check to see max number of columns for R import as wide format
	max_paralog_codons = ""
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


def R_file_prep2(input_file2):	#REDUNDANT DUE TO TABLEIZE! Do not use!
	#input_file2 = sys.argv[2]	#path of vcf file containing only orginal query variants list used as input for VEP_ParalogAnno.py. I.e. the very first input file! Can be "null" if using tableize instead
	if input_file2 != "null":
		if input_file2.endswith(".vcf"):
			out_file2 = open(input_file2+"_onlyVariantslist", "w")
			with open(input_file2, "r") as f:
				for line in f:
					if not line.startswith("#"):
						chrom = line.split()[0]
						position = line.split()[1]
						ID = line.split()[2]
						out_file2.write(str(chrom)+" "+str(position)+"\t"+str(ID)+"\n")

		out_file2.close()