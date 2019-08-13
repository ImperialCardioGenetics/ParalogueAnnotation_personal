import  os, sys, subprocess, re,  shlex, subprocess, codecs, pysam, gzip, argparse

#SCRIPT FOR FILTERING OUT RARE VARIANTS DEFINED AS NOT BEING IN discovEHR

#to tabix index files - tabix -p vcf input.vcf.gz

gnomad_file = sys.argv[1]	#e.g. /work/nyl112/data/gnomad_split_by_line/gnomad.r2.1.1/gnomad.exomes.r2.1.1.sites.vcf_split00.gz
out_file_name = gnomad_file.rsplit(".",1)[0]+"_notin_discovEHR.vcf"

out_file = open(out_file_name, "w")
out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

discovEHR_file = pysam.TabixFile("/work/nyl112/data/DiscovEHR/discovEHR_GRCh37.vcf.gz") #path of tabix indexed file (.tab.gz file) #NOT tbi file

# def get_info_value(db_file,chrom,pos,ref_nt,alt_nt,ref_aa,alt_aa):
#     # position_found = False
#     REVEL= ""
#     for row in db_file.fetch(chrom, pos - 1, pos):
#         row_column = row.split('\t')
#         if str(pos) != row_column[1]:
#             continue
#         # position_found = True
#         row_ref_nt = row_column[2]
#         row_alt_nt = row_column[3]
#         row_ref_aa = row_column[4]
#         row_alt_aa = row_column[5]
#         #print row_ref_nt,row_alt_nt,row_ref_aa,row_alt_aa
#         #print ref_nt,alt_nt,ref_aa,alt_aa
#         if row_ref_nt == ref_nt and row_alt_nt == alt_nt and row_ref_aa == ref_aa and row_alt_aa ==alt_aa:
#             REVEL = row_column[6]
#         elif row_ref_nt == ref_nt and row_alt_nt == alt_nt:
#             REVEL = row_column[6]
#     # print(REVEL)
#     while not isinstance(REVEL, str):
#         REVEL = REVEL[0]
#     REVEL = [REVEL]
#     return REVEL

# chrom = 1
# pos = 69088

with gzip.open(gnomad_file, "rt") as f:
	for line in f:
		if not line.startswith("#"):
			split_line = line.split()
			FILTER = split_line[6]
			if FILTER == "PASS":
				CHROM = split_line[0]
				POS = split_line[1]
				REF = split_line[3]
				ALT = split_line[4]
				
				for row in discovEHR_file.fetch(CHROM, int(POS) - 1, int(POS)):
					row = row.split()
					disc_CHROM = row[0]
					disc_POS = row[1]
					# disc_ID = row[2]
					disc_REF = row[3]
					disc_ALT = row[4]
					if CHROM == disc_CHROM and POS == disc_POS and REF == disc_REF and ALT == disc_ALT:
						

						ID = split_line[2]
						QUAL = split_line[5]
						INFO = split_line[7]
						AC = line.split(";AN=")[0].split("AC=")[1]
						AN = line.split(";AN=")[1].split(";AF=")[0]
						AF = line.split(";AF=")[1].split(";rf_tp_probability=")[0]
						nhomalt = line.split(";nhomalt=")[1].split(";")[0]

						out_file.write()






out_file.close()






