import  os, sys, subprocess, re, shlex, subprocess, codecs
from VEP_ParalogAnno import *
from File_prep_for_R import *
from Tableize_wrapper import *

input_file = sys.argv[1]	#name of input file

print(len(sys.argv), sys.argv)

if len(sys.argv) == 4:
	flavour = int(sys.argv[3])
else:
	flavour = 2

VEP_Plugin_run(input_file, genome_build="GRCh38", flavour=flavour, VEPversion=93, offline=1)#, output_filename=sys.argv[2])
Tableize_wrap(input_file.rsplit(".",1)[0]+".out_paraloc")
R_file_prep(input_file.rsplit(".",1)[0]+".out_paraloc_paralogs", "noQC")
