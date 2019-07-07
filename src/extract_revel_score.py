##Goal:
#given a list of genomic coordinates of variants, return the requested field infomation
#e.g given the genomic coordinates of variants (chr_pos_ntref_ntalt_aaref_aaalt), return the requested field information

#Annotate query variant file with REVEL scores

import os, sys, subprocess, re, shlex, subprocess, codecs, pysam, gzip, argparse

input_file1 = sys.argv[1]   #path of query variants tableized input file
input_file2 = sys.argv[2]   #path of tabix indexed file (.tab.gz file)
tabix_file = pysam.TabixFile(input_file2)
out_file_name = input_file1+"_w_REVEL"

out_file = open(out_file_name, "w")

def get_info_value(db_file,chrom,pos,ref_nt,alt_nt,ref_aa,alt_aa):
    # position_found = False
    REVEL= ""
    for row in db_file.fetch(chrom, pos - 1, pos):
        row_column = row.split('\t')
        if str(pos) != row_column[1]:
            continue
        # position_found = True
        row_ref_nt = row_column[2]
        row_alt_nt = row_column[3]
        row_ref_aa = row_column[4]
        row_alt_aa = row_column[5]
        #print row_ref_nt,row_alt_nt,row_ref_aa,row_alt_aa
        #print ref_nt,alt_nt,ref_aa,alt_aa
        if row_ref_nt == ref_nt and row_alt_nt == alt_nt and row_ref_aa == ref_aa and row_alt_aa ==alt_aa:
            REVEL = row_column[6]
        elif row_ref_nt == ref_nt and row_alt_nt == alt_nt:
            REVEL = row_column[6]
    # print(REVEL)
    while not isinstance(REVEL, str):
        REVEL = REVEL[0]
    REVEL = [REVEL]
    return REVEL

with open(input_file1, "r") as f:
    for line in f:
        if line.startswith("CHROM"):
        	out_file.write(line.rstrip()+"\tREVEL_score\n")
        if not line.startswith("CHROM"):
            line = line.split()
            if not "/" in line[8]:
            	out_file.write(str("\t".join(line+["NA"]))+"\n")
            	continue
            chrom = line[0]
            pos = int(line[1])
            ntref = line[2]
            ntalt = line[3]
            aaref = line[8].split("/")[0]
            aaalt = line[8].split("/")[1]
            print(chrom,pos,ntref,ntalt,aaref,aaalt)
            REVEL_score = get_info_value(tabix_file,chrom,pos,ntref,ntalt,aaref,aaalt)
            print(REVEL_score)
            new_line = "\t".join(line+REVEL_score)
            out_file.write(str(new_line)+"\n")

tabix_file.close()
out_file.close()