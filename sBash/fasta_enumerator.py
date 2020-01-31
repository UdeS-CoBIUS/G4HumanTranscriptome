#!/usr/bin/env python
from Bio import SeqIO
import sys
import re

def splitter(fasta_file, output, limit, large_handling=False):
    """
    Splits a large fasta_file in sub-files created in a given directory.
    """
    file_ = open(fasta_file, 'r')
    file_count = 1
    outfile = open(output.rstrip("/")+"/%s_%05d.fas"%(
        fasta_file.split('/')[-1].split('.')[0],file_count),'w')
    nt_count = 0
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        if large_handling == True and len(str(seq.seq)) >= int(limit):
            file_count += 1
            largefile = open(output.rstrip("/")+"/%s_%05d_XL.fas"%(
                fasta_file.split('/')[-1].split('.')[0],file_count),'w')
            largefile.write(">"+str(seq.description)+"\n"+"\n".join(
                str(seq.seq)[i:i+50]for i in xrange(0,len(seq.seq),50))+"\n")
            largefile.close()
        else:
            nt_count += len(str(seq.seq))
            outfile.write(">"+str(seq.description)+"\n"+"\n".join(
                str(seq.seq)[i:i+50]for i in xrange(0,len(seq.seq),50))+"\n")            
            if nt_count >= int(limit):
                outfile.close()
                file_count += 1
                nt_count = 0
                outfile = open(output.rstrip("/")+"/%s_%05d.fas"%(
                    fasta_file.split('/')[-1].split('.')[0],file_count),'w')
    outfile.close()


def usage(error_message=False):
    """
    Provide the user with instructions to use fasta_enumerator.py.    
    """
    print "Usage: python PATH/TO/fasta_enumerator.py [OPTIONS...]"
    print "Use -h, -? or --help to show this message\n"
    print "Options:"
    print "  -f, --fasta-file\tSupply a fasta file"
    print "  -n, --nt-limit  \tNumber of nts to set the split limit"
    print "  -o, --output    \tOutput directory"
    print "  -l, --large     \tTake large sequence to a single file"
    print "                  \tLarge sequences are larger than provided limit"
    print "  -e, --error     \tRaise errors and exceptions\n"
    if error_message:
        print "UsageError:", error_message
        sys.exit(500)
    else:
        sys.exit(0)

def main():
    """
    Handles arguments.
    """
    option_dict = {}
    for no, arg in enumerate(sys.argv):
        if arg[0] == "-":
            if arg in ["-h","-?","--help"]:
                usage()
            elif arg in ["-l","--large",
                    "-v","--verbose",
                    "-e","--error"]:
                option_dict[arg] = True
            elif arg in ["-f","--fasta-file",
                    "-n","--nt-limit",
                    "-o","--output"]:
                option_dict[arg] = sys.argv[no+1]
            else:
                usage('Argument "%s" not recognized'%arg)
    if len(sys.argv) == 1:
        usage()
    try:
        splitter(option_dict.get("-f") or option_dict.get("--fasta-file"),
                option_dict.get("-o") or option_dict.get("--output"),
                option_dict.get("-n") or option_dict.get("--nt-limit"),
                option_dict.get("-l") or option_dict.get("--large"))
    except:
        if "-e" in option_dict.keys() or "--error" in option_dict.keys():
            raise
        else:
            usage('An option is missing, incorrect or not authorized')

if __name__ == '__main__':
    main()
