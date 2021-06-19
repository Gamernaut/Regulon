# Module for reading restriction enzyme file and FASTA DNA sequence files

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
                line = line.replace('\n', '')  # remove the newline character
                re_seq = line.split(',')[0]
                # Filter out the non alphabet characters such as brackets, slash and numbers etc.
                re_seq = ''.join(filter(str.isalpha, re_seq))
                re_name = line.split(',')[1]
                # sequence_dict[re_seq] = re_name
                sequence_dict[re_name] = re_seq
                sequence_length = len(re_seq)
                if tree_depth < sequence_length:
                    tree_depth = sequence_length
                count += 1
        RE_file.close()
    except Exception as err:
        print(f"\nERROR - FileHandler.import_restriction_enzymes() had a problem with the file: {filename}.\nError was: ", err)
    return sequence_dict, tree_depth


def import_seq_file(filename):
    """Reads the file with the DNA sequence. Assumes a standard FASTA file format which is info on the first line
    marked with a greater than symbol and then the DNA sequence is all the following lines concatenated together.
    This function returns the DNA sequence as a string."""
    dna_sequence = ""

    try:
        with open(filename, "r") as seq_file:
            for line in seq_file:
                if line.startswith(">"):
                    continue
                else:
                    line = line.replace('\n', '')  # remove the newline character
                    dna_sequence += line
        seq_file.close()
    except Exception as err:
        print(f"\nERROR - FileHandler.import_seq_file() had a problem with the file: {filename}.\nError was: ", err)

    return dna_sequence
