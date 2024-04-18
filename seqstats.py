# Class to store the sequence data and calculate statistics about it
# Takes a dictionary of DNA sequences.  Each key is a sequence ID, each value is a list of sequences.
class SeqStats:
    # on initialization, calculate metrics
    def __init__(self, seq_dict):
        self.data = seq_dict
        self.seq_len_dict = self.len_sequences()
        self.longest_seq_len = 0
        self.get_longest_seq()
        self.longest_seq_list = []
        self.create_longest_seq_list()
        self.shortest_seq_len = 0
        self.get_shortest_seq()
        self.shortest_seq_list = []
        self.create_shortest_seq_list()
        self.start_codon = 'ATG'
        self.stop_codons = ['TAA', 'TAG', 'TGA']
        self.orf_dict = {}
        self.longest_orf_length = 0
        self.longest_orf = 0
        self.longest_orf_seq_id = 0
        self.longest_orf_seq_length = 0
        self.longest_orf_seq = 0
        self.longest_orf_seq_start = -1
        self.repeat_dict = {}
        self.highest_repeat_dict = {}


    # get number of records
    def num_records(self):
        return len(self.data.keys())


    # create dictionary to store sequence lengths
    def len_sequences(self):
        len_dict = {}

        for key in self.data.keys():
            len_dict[key] = len(self.data[key])

        # return dictionary
        return len_dict


    # get longest sequence in file
    def get_longest_seq(self):
        # find
        for key in self.seq_len_dict.keys():
            # if sequence length greater than stored sequence length
            if self.seq_len_dict[key] > self.longest_seq_len:
                # save max sequence length
                self.longest_seq_len = self.seq_len_dict[key]


    # save all keys with this longest sequence length
    def create_longest_seq_list(self):
        for key in self.seq_len_dict.keys():
            # if length is max seq length
            if self.seq_len_dict[key] == self.longest_seq_len:
                self.longest_seq_list.append(key)


    # get shortest sequence in file
    def get_shortest_seq(self):
        # find shortest seq
        for key in self.seq_len_dict.keys():
            # if a length hasn't been defined yet, use the current length as a base
            if self.shortest_seq_len == 0:
                self.shortest_seq_len = self.seq_len_dict[key]
            elif self.seq_len_dict[key] < self.shortest_seq_len:
                # if sequence length smaller than stored sequence length
                # save shortest sequence length
                self.shortest_seq_len = self.seq_len_dict[key]


    # save all keys with this shortest sequence length
    def create_shortest_seq_list(self):
        for key in self.seq_len_dict.keys():
            if self.seq_len_dict[key] == self.shortest_seq_len:
                self.shortest_seq_list.append(key)


    # find all ORFs
    def find_orfs(self, reading_frame):
        # if reading frame is a weird value, print usage and end
        if reading_frame < 1 or reading_frame > 3:
            print("Please enter a reading frame of 1, 2, or 3")
            return

        for seq_id, seq in self.data.items():
            start_idx = reading_frame - 1

            # start blank list of strings
            self.orf_dict[seq_id] = []

            for i in range(start_idx, (len(seq) - 2), 3):
                codon = seq[i:i + 3]

                # if start codon, store
                if codon == self.start_codon:
                    exon = codon

                    # advance i
                    i = i + 3

                    # get the rest of the sequence
                    for j in range(i, (len(seq) - 2), 3):
                        # Question for instructor: do we need to handle a start codon in a sequence without a stop codon?  Does this ever occur?
                        codon = seq[j:j + 3]
                        
                        # add codon
                        exon = exon + codon

                        # if stop codon, break
                        if codon in self.stop_codons:                            
                            # save string to list
                            self.orf_dict[seq_id].append(exon)

                            # increment read index
                            i = j + 3

                            # exit loop
                            break

        # get longest ORF stats for entire file
        self.get_longest_orf()


    # get the starting position in the sequence string of a sub-sequence
    def get_subseq_start_pos(self, seq_id, seq):
        # get full sequence from seq_id
        full_seq = self.data[seq_id]

        # find subseq in seq
        start_pos = full_seq.find(seq)

        # if sequence not found, print error but still return
        if start_pos == -1:
            print("Sequence ", seq, " not found for sequence ID: ", seq_id)

        # return starting position (we want character position, not index)
        return start_pos + 1


    # get longest ORF in file
    def get_longest_orf(self, seq_id=0):
        # if given seq id, calculate ORF metrics for that sequence only
        if seq_id != 0:
            seq_list = self.orf_dict[seq_id]

            for seq in seq_list:
                # if current seq longer than max seq, update longest_orf_seq
                if len(seq) > self.longest_orf_seq_length:
                    self.longest_orf_seq_length = len(seq)
                    self.longest_orf_seq = seq
                    # find starting character in original sequence
                    self.longest_orf_seq_start = self.get_subseq_start_pos(seq_id, seq)
        else: # calculate ORF metrics for entire file
            # iterate through orf dict list
            for key in self.orf_dict.keys():
                seq_list = self.orf_dict[key]

                for seq in seq_list:
                    # if current seq longer than max seq, update longest_orf and note position
                    if len(seq) > self.longest_orf_length:
                        self.longest_orf_length = len(seq)
                        self.longest_orf = seq
                        self.longest_orf_seq_id = key

    
    # get all repeats (overlapping) of length n
    # calculate how many times each repeat occurs in the file
    def get_repeats(self, n):
        # get DNA sequence
        for key in self.data.keys():
            seq = self.data[key]

            # iterate through sequence
            for i in range(0, (len(seq) - n)):
                subseq = seq[i:i+n]

                # if sub-sequence is already a known repeat, increment counter
                if subseq in self.repeat_dict.keys():
                    self.repeat_dict[subseq] = self.repeat_dict[subseq] + 1
                else:
                    # add subsequence to dict
                    self.repeat_dict[subseq] = 1

        # delete all dictionary entries where repeat count is 1
        new_rep_dict = { key:value for key, value in self.repeat_dict.items() if value > 1 }
        self.repeat_dict = new_rep_dict

        # calculate highest repeat frequency
        self.get_highest_repeats()


    # most frequent repeat of a given length
    def get_highest_repeats(self):
        max_repeat_freq = 0

        for key in self.repeat_dict.keys():
            # if repeat frequency greater than stored repeat frequency
            if self.repeat_dict[key] > max_repeat_freq:
                # update max repeat frequency
                max_repeat_freq = self.repeat_dict[key]

        # save all repeats of max frequency in a dictionary
        for key in self.repeat_dict.keys():
            if self.repeat_dict[key] == max_repeat_freq:
                self.highest_repeat_dict[key] = self.repeat_dict[key]


    # print stats
    def print_seq_stats(self):
        # number of records
        print("Number of records: ", self.num_records())
        # all sequence lengths
        print("Sequence lengths: ", self.seq_len_dict)
        # longest sequence length
        print("Longest sequence length: ", self.longest_seq_len)
        # seqs with max length
        print("Total sequences of longest length: ", len(self.longest_seq_list))
        # keys of seqs with max length
        print("Longest sequences: ", self.longest_seq_list)
        # shortest sequence length
        print("Shortest sequence length: ", self.shortest_seq_len)
        # seqs with min length
        print("Total sequences of shortest length: ", len(self.shortest_seq_list))
        # keys of seqs with min length
        print("Shortest sequences: ", self.shortest_seq_list)


    # print ORF-specific stats
    def print_orf_stats(self, seq_id=0):
        # if only self is entered, assume stats for all ORFs in file
        if seq_id == 0:
            # longest ORF length
            print("Longest ORF in file: ", self.longest_orf_length, " characters")
            # longest ORF
            print("Longest ORF sequence in file: ", self.longest_orf)
            # longest ORF seq ID
            print("Sequence ID of longest ORF: ", self.longest_orf_seq_id)
        else: # print stats for specific sequence ID only (must have already been created)
            # longest ORF in seq length
            print("Longest ORF in sequence ", seq_id, ": ", self.longest_orf_seq_length, " characters")
            # longest ORF in seq
            print("Longest ORF in sequence ", seq_id, ": ", self.longest_orf_seq)
            # longest ORF starting position
            print("Longest ORF starting position in sequence  ", seq_id, ": ", self.longest_orf_seq_start)


    # print repeat-specific stats
    def print_repeat_stats(self):
        # all repeat frequencies of given length
        print("All repeat frequencies: ", self.repeat_dict)
        # highest repeat frequencies
        print("Highest repeat frequencies: ", self.highest_repeat_dict)


# Input: opened file
# Output: dictionary of sequences
def create_sequence_dict(f):
    seq_dict = {}

    for line in f:
        # discard newlines
        line = line.rstrip()

        if line[0] == '>': # if header
            # split on whitespace
            words = line.split()

            # name is after the >
            name = words[0][1:]

            # create new entry in dict
            seq_dict[name] = ''
        else: # not header
            # update seq in dict
            seq_dict[name] = seq_dict[name] + line

    # check if there are any sequences in the dict
    if len(seq_dict.keys()) <= 0:
        sys.exit("***Error***: No DNA sequences found in file.  Please ensure that your file is in FASTA format.")
    else:
        return seq_dict