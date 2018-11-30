#One-time use of Tableize_wrapper and File_prep_for_R to compare Paralogue annotation with SIFT and REVEL scores. Not intended to be part of the pipeline

import os, sys, subprocess, re, shlex, subprocess, codecs

input_file1 = sys.argv[1] #Tableized file
input_file2 = sys.argv[2] #paralogs file