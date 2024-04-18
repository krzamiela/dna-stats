# dna-stats
Final project for Coursera's Python for Genomics course, offered through Johns Hopkins.
Reads a multi-FASTA file and performs operations on extracted data.

## Requirements
- Python 3

## Usage
    pygenfinal.py [-h] <filename>

    -h              print this message
    
    <filename>      the file has to be in FASTA format

## Summary
This program was developed as a final project, and parts of it could be improved for better usage.  For example, I chose to wrap all statistics about the given multi-FASTA file into a class.  I could have split the class into a separate module and provided an example Python file showing how to interface with the class.

## SeqStats Class Basics
To create a new instance:
    foo = SeqStats(seq_dict)

**Note:** The variable seq_dict must be formatted such that each key represents a sequence ID, and each value is a list of sequences matching that sequence ID

Upon creation, the SeqStats class automatically calculates some statistics for the given seq_dict, which can be summarized by calling print_seq_stats():
    foo.print_seq_stats()

This prints:

- Number of records
- A dictionary of sequence lengths
- The longest sequence length
- The total number of sequences of the longest sequence length
- A list of all sequences with the longest sequence length
- The shortest sequence length
- The total number of sequence of the shortest sequence length
- A list of all sequences with the shortest sequence length

To calculate ORF statistics, invoke the find_orfs() method:
    foo.find_orfs(1)

**Note:** find_orfs() expects an integer of either 1, 2, or 3, to specify the reading frame.

To print ORF statistics based on the current reading frame, use print_orf_stats():
    foo.print_orf_stats()

This prints:

- Length of the longest ORF in the file
- Longest ORF sequence in the file
- Sequence ID of longest ORF

To calculate ORF statistics for a particular sequence ID only, enter that sequence ID as an optional second argument:
    foo.get_longest_orf('bar')

To print ORF statistics based on the current reading frame for a particular sequence ID only, enter that sequence ID as an optional second argument:
    foo.print_orf_stats('bar')

This prints:

- Length of the longest ORF in the given sequence
- Longest ORF in the given sequence
- Starting position in sequence of the longest ORF

To find all repeats of length n in the file, use get_repeats():
    foo.get_repeats(n)

**Note:** get_repeats() expects a positive integer value

To print repeat statistics of length n, use print_repeat_stats():
    foo.print_repeat_stats()

This prints:

- A dictionary of all repeats and their frequencies
- A dictionary of the highest repeat frequency and all sequences with this frequency