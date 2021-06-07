# Module for reading restriction enzyme file and FASTA DNA sequence files

# used to remove text between brackets in imported sequences
import re


def import_restriction_enzymes(filename):
    """Reads the file with restriction enzyme recognition sequences and names, parses them and adds them to a
    dictionary. Returns the dictionary and the length of the longest sequence"""
    sequence_dict = {}
    count = 0
    tree_depth = 0
    try:
        with open(filename, "r") as RE_file:
            for line in RE_file:
                if count == 0:
                    count += 1
                    continue
                #line = line.replace('/', '')  # remove the cutting site in the sequence
                line = line.replace('\n', '')  # remove the newline character
                #line = re.sub("[\(\[].*?[\)\]]", "",
                #              line)  # removes brackets and the contents as some sequences have extra
                # offset info in brackets for example ACCTGC(4/8),BspMI/BfuAI
                # TODO: Split the line then remove brackets so the / between RE names are preserved
                # could also try this s = ''.join(filter(str.isalnum, s)) to strip out everything that's not a character
                re_seq = line.split(',')[0]
                re_seq = ''.join(filter(str.isalpha, re_seq))
                re_name = line.split(',')[1]
                # sequence_dict[re_seq] = re_name
                sequence_dict[re_name] = re_seq
                sequence_length = len(line.split(',')[0])
                if tree_depth < sequence_length:
                    tree_depth = sequence_length
                count += 1
        RE_file.close()
    except Exception as err:
        print(f"\nERROR - FileHandler.import_restriction_enzymes() had a problem with the file: {filename}.\nError was: ", err)
    return sequence_dict, tree_depth


def import_seq_file(seq_filename):
    # TODO: implement the FASTA file read
    pass
    return
