#!/usr/bin/python
# Creator: Kathryn Zamiela
# Description: Takes command-line input of a file containing DNA sequences in multi-FASTA format,
#               extracts the data, and performs operations on the data
# Questions for instructor: Do you ever run into sequences where you encounter a start codon and then another start codon, or a stop codon and then another stop codon?  Or are all DNA sequences start, stuff, stop?

import sys
import getopt
import seqstats

# Input: none
# Output: string with the file path
def check_args():
    o, a = getopt.getopt(sys.argv[1:], 'h')

    opts = {}

    for k, v in o:
        opts[k] = v

    # if user asked for help, print usage and exit
    if '-h' in opts.keys():
        usage(); sys.exit() # multiple instructions on same line if sep by ;

    # check for requirements
    if len(a) < 1:
        usage(); sys.exit("***Error***: Input FASTA file is missing")
    elif len(a) == 1:
        # if the correct number of argument are present, return file name
        return a[0] # only returning first element in the list since we know the list length is 1
    else:
        usage(); sys.exit("***Error***: Incorrect number of arguments provided")
    

# Input: a string path to a file containing DNA sequences in multi-FASTA format
# Output: file handler
def open_file(filename):

    try:
        f = open(filename)
    except IOError:
        print(filename); sys.exit("File does not exist")

    return f


# Input: none
# Output: text of how to use the file
def usage():
    print("""
    pygenfinal.py : reads a FASTA file and performs operations on extracted data

    pygenfinal.py [-h] <filename>

    -h              print this message
    
    <filename>      the file has to be in FASTA format
    """)


## MAIN ##
# check if arguments provided are valid
fname = check_args()

# open file
ufile = open_file(fname)

# create dictionary of sequences
seqs = seqstats.create_sequence_dict(ufile)

# create SeqStats class and get metrics
all_stats = seqstats.SeqStats(seqs)

# print stats
all_stats.print_seq_stats()

# identify all ORFs
# prompt for reading frame
reading_frame = 0

while reading_frame not in ['1', '2', '3']:
    reading_frame = input("At which index would you like to start the reading frame: 1, 2, or 3? ")

# what is the longest ORF in the file? what is the identifier of the containing sequence?
all_stats.find_orfs(int(reading_frame))
all_stats.print_orf_stats()

# For a given sequence identifier, what is the longest ORF contained in the sequence?
orf_seq = ''
while orf_seq != 'y' and orf_seq != 'n':
    orf_seq = input("Do you want to calculate metrics for the longest ORF of a particular sequence? (y/n) ")

if orf_seq == 'y':
    test_seq_id = input("Please input the sequence ID you would like to test: ")

    all_stats.get_longest_orf(test_seq_id)
    all_stats.print_orf_stats(test_seq_id)

n = ""

while not n.isnumeric():
    n = input("Please enter repeat length as an integer: ")

# Given a length n, identify all repeats of length n in all sequences in the FASTA file
all_stats.get_repeats(int(n))

all_stats.print_repeat_stats()