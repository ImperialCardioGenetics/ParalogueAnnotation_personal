import os, sys, subprocess, re, shlex, subprocess, codecs
from VEP_ParalogAnno import *
from File_prep_for_R import *
from Tableize_wrapper import *

input_file = sys.argv[1]	#path of input file

build = sys.argv[2]	#genome build: either 37 or 38
if int(build) == 37:
	genome_build = "GRCh37"
elif int(build) == 38:
	genome_build = "GRCh38"
else:
	sys.exit("ERROR: genome build not recognized")

flavour = sys.argv[3]
if flavour == "base":
	flavour = 0
elif flavour == "variant":
	flavour = 1
elif flavour == "paraloc":
	flavour = 2
else:
	sys.exit("ERROR: mode selected not recognized")

print(len(sys.argv), sys.argv)


VEP_Plugin_run(input_file, genome_build=genome_build, flavour=flavour, VEPversion=93, offline=1)#, output_filename=sys.argv[2])
Tableize_wrap(input_file.rsplit(".",1)[0]+".out_paraloc")
R_file_prep(input_file.rsplit(".",1)[0]+".out_paraloc_paralogs", "noQC")
