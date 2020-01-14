import os, sys, subprocess, re, shlex, subprocess, codecs
from VEP_ParalogAnno import *
from File_prep_for_R import *
from Tableize_wrapper_with_ParalogueVars import *

input_file = sys.argv[1]	#path of input file

build = sys.argv[2]	#genome build: either 37 or 38
if int(build) == 37:
	genome_build = "GRCh37"
elif int(build) == 38:
	genome_build = "GRCh38"
else:
	sys.exit("ERROR: genome build not recognized")

flavour = sys.argv[3] #flavour of plugin
if flavour == "base":
	flavour = 0
elif flavour == "variant":
	flavour = 1
elif flavour == "paraloc": #NORMALLY paraloc IS WHAT YOU WANT
	flavour = 2
else:
	sys.exit("ERROR: plugin mode selected not recognized")

refid_flavour = sys.argv[4] #flavour of refid filtering - how conserved do you want the reference alleles to be?
if refid_flavour == "noQC":
	refid_flavour = "noQC"
elif refid_flavour == "para_con":
	refid_flavour = "para_con"
elif refid_flavour == "all_con":
	refid_flavour = "all_con"
else:
	sys.exit("ERROR: refid flavour selected not recognized")

tableize_dir = sys.argv[5] #path of where tableize_vcf.py is located, e.g. on HPC - "/work/nyl112/loftee/src/"; can have "/" at end or not; on RBH - "/data/Share/nick/Paralog_Anno/loftee/src/"

para_zscore_dir = sys.argv[6] #path of where para_zscore folder is located, e.g. on HPC - "/work/nyl112/data"; can have "/" at end or not; on RBH - "/data/Share/nick/Paralog_Anno/data_files/"

print(len(sys.argv), sys.argv)


VEP_Plugin_run(input_file, genome_build=genome_build, flavour=flavour, VEPversion=93, offline=1)#, output_filename=sys.argv[2])
Tableize_wrap(input_file.rsplit(".",1)[0]+".out_paraloc", tableize_dir, para_zscore_dir)
R_file_prep(input_file.rsplit(".",1)[0]+".out_paraloc_paralogs", refid_flavour=refid_flavour)
