#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	This script will run one blastn per file in the query_directory. In our
        study it is used to run blastn for pG4r vs rG4 and rG4 vs pG4r.

Usage:
    python ~/PATH/RunBlast.py query_directory fasta_target output_Directory
"""

import os
import sys

def do_blast (query_dir, db, outdir):
    filelist = os.listdir(query_dir)
    for queryfile in filelist:
		#~ print str(queryfile)
        out = outdir + queryfile.split('.')[0] + '.blastresult.xml'
        cmd = "blastn -query " + query_dir + queryfile + " -db " + db +\
         " -out "+ out + " -outfmt 5" + " -word_size 11"
        #~ print cmd
        os.system(cmd)

if __name__ == '__main__':
    do_blast(sys.argv[1], sys.argv[2], sys.argv[3])
