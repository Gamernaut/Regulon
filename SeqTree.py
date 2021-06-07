# Module to create and search the tree based data structure of restriction enzyme sequences

# Import standard modules
import sys
import time
import pprint
# import utility module to handle reading the source files
import FileHandler


class Node(object):
    """DNA Tree node with support for the 4 basic nucleotides, a list of enzymes that match this sequence
            and knowledge of whether it's the last node in the branch or not"""

    def __init__(self, base=None):
        self.A = None
        self.C = None
        self.G = None
        self.T = None
        self.nucleotide = base
        self.RE_list = []  # supports several RE names for the same sequence
        self.is_branch_end = None


class RESeqTree:
    """This is the class that represents the restriction enzyme (RE) sequence tree which is made up of the Node class"""

    def __init__(self):
        self.root = Node()
        self.tree_width = 0
        self.tree_depth = 0

        # Next 2 attributes are used by the Print Tree function
        self.branch_sequence = ""
        self.tree_sequences = []

        # Following is to support new dict based approach to print tree
        self.re_seq_dict = {}


        # Used to calculate progress in loading Restriction Enzyme definition file and creating tree
        self.sequence_count = 0
        self.unique_sequence_count = 0

        # If a new request to create tree with same file then we can just pass the root back (singleton type idea)
        self.re_filename = ""
        # More robust test in case the content of the file has changed but not the name, MD5 will pick up the difference
        self.re_file_MD5_checksum = 0
        # Currently pulled out as Result_Manager as a separate implementation,
        # something that manages the matches for this combination and run
        # self.matches = dict()

    def get_root(self):
        return self.root

    def get_tree_width(self):
        return self.tree_width

    def get_tree_depth(self):
        return self.tree_depth

    def get_sequence_count(self):
        return self.sequence_count

    def build_tree(self, filename):
        # Check if the file being passed to us is the same as the existing source for the current tree.
        # If so just pass back the root which is a simple way to support multiple runs in a future version.
        # TODO: Using a hashcode to check the content of the file is a much better way to ensure that the file
        # contents have not been changed even if the name stays the same

        if filename == self.re_filename:
            return self.get_root()

        # It's a new file so we need a new tree, set the name so the check above will catch the reuse of the file
        self.re_filename = filename

        # Call the function to load the sequences from the restriction enzyme definition file
        self.re_seq_dict, self.tree_depth = FileHandler.import_restriction_enzymes(filename)

        # TODO: can we use the re_seq_dict as part of the print tree method?? Maybe print out the RE name then the sequences

        # Calculate the progress of the tree build and display this to the user
        current_seq_count = 1
        sequences_to_build = len(self.re_seq_dict)
        for re_name, seq in self.re_seq_dict.items():
            pct_complete = (current_seq_count/sequences_to_build) * 100
            sys.stdout.write(f"\rProcessing Restriction Enzymes sequences from {filename}: %d%% completed " % pct_complete)  # print on the same line
            self.insert_sequence(self.root, seq, re_name)
            current_seq_count += 1
            time.sleep(0.02)
        print()
        # print(re_seq_dict)  # DEBUG

    def insert_sequence(self, node, sequence, name):
        if len(sequence) <= 0:
            node.RE_list.append(name)
            node.is_branch_end = True
            self.unique_sequence_count += 1
            self.tree_width += 1
            return
        current_base = sequence[0]
        new_sequence = sequence[1:]
        # These are the 4 basic nucleotides
        if current_base == "A":
            self.insert_A(current_base, name, new_sequence, node)
        elif current_base == "C":
            self.insert_C(current_base, name, new_sequence, node)
        elif current_base == "G":
            self.insert_G(current_base, name, new_sequence, node)
        elif current_base == "T":
            self.insert_T(current_base, name, new_sequence, node)

        # Support for DNA ambiguity codes from https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
        # The following are combinations of 2 nucleotides
        elif current_base == "Y":  # Pyrimidine (C or T)
            self.insert_C("C", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)
        elif current_base == "R":  # Purine (A or G)
            self.insert_A("A", name, new_sequence, node)
            self.insert_G("G", name, new_sequence, node)
        elif current_base == "W":  # Weak (A or T)
            self.insert_A("A", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)
        elif current_base == "S":  # Strong (G or C)
            self.insert_C("C", name, new_sequence, node)
            self.insert_G("G", name, new_sequence, node)
        elif current_base == "K":  # Keto (T or G)
            self.insert_G("G", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)
        elif current_base == "M":  # Amino (A or C)
            self.insert_A("A", name, new_sequence, node)
            self.insert_C("C", name, new_sequence, node)

        # The following are combinations of 3 nucleotides
        elif current_base == "D":  # A, G and T (not C)
            self.insert_A("A", name, new_sequence, node)
            self.insert_G("G", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)
        elif current_base == "V":  # A, C and G (not T)
            self.insert_A("A", name, new_sequence, node)
            self.insert_C("C", name, new_sequence, node)
            self.insert_G("G", name, new_sequence, node)
        elif current_base == "H":  # A, C and T (not G)
            self.insert_A("A", name, new_sequence, node)
            self.insert_C("C", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)
        elif current_base == "B":  # C, G and T (not A)
            self.insert_C("C", name, new_sequence, node)
            self.insert_G("G", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)

        # The following is all 4 nucleotides
        elif current_base == "N" or current_base == "X":  # All nucleotides (A, C, G and T)
            self.insert_A("A", name, new_sequence, node)
            self.insert_C("C", name, new_sequence, node)
            self.insert_G("G", name, new_sequence, node)
            self.insert_T("T", name, new_sequence, node)
        else:
            print("Nucleotide", current_base, "not supported in SeqTree.insert_sequence()")

    def insert_A(self, current_base, name, new_sequence, node):
        if node.A is None:
            new_node = Node(current_base)
            node.A = new_node
            self.insert_sequence(new_node, new_sequence, name)
        else:
            self.insert_sequence(node.A, new_sequence, name)

    def insert_C(self, current_base, name, new_sequence, node):
        if node.C is None:
            new_node = Node(current_base)
            node.C = new_node
            self.insert_sequence(new_node, new_sequence, name)
        else:
            self.insert_sequence(node.C, new_sequence, name)

    def insert_G(self, current_base, name, new_sequence, node):
        if node.G is None:
            new_node = Node(current_base)
            node.G = new_node
            self.insert_sequence(new_node, new_sequence, name)
        else:
            self.insert_sequence(node.G, new_sequence, name)

    def insert_T(self, current_base, name, new_sequence, node):
        if node.T is None:
            new_node = Node(current_base)
            node.T = new_node
            self.insert_sequence(new_node, new_sequence, name)
        else:
            self.insert_sequence(node.T, new_sequence, name)

    # TODO: Implement search function
    def find_matches(self, filename, result_manager):
        pass

    def print_branch(self, node):
        if node is None:
            print("No valid tree node provided to SeqTree.print_branch()")
            return
        # Add the nucleotide in the current node to the sequence constructed so far (from previous calls)
        # If it's the root there is no nucleotide so skip this
        if node.nucleotide:
            self.branch_sequence = self.branch_sequence + node.nucleotide
        if node.is_branch_end:
            # TODO: Change code to have RE name first then list sequences
            # eg. BsaJI
            #       CCCGGG
            #       CCCTGG
            temp = self.branch_sequence+" -> "+', '.join(node.RE_list)  # replaces [] from std List print with ","
            self.tree_sequences.append(temp)

        if node.A:
            self.print_branch(node.A)
            # as you come up the call stack remove this nucleotide
            self.branch_sequence = self.branch_sequence[:-1]
        if node.C:
            self.print_branch(node.C)
            self.branch_sequence = self.branch_sequence[:-1]
        if node.G:
            self.print_branch(node.G)
            self.branch_sequence = self.branch_sequence[:-1]
        if node.T:
            self.print_branch(node.T)
            self.branch_sequence = self.branch_sequence[:-1]

    def print_tree(self):
        root = self.get_root()
        if root:
            self.print_branch(root)
        else:
            print("ERROR - Can't get root for tree")
            return
        print(f"\nThis tree has {self.unique_sequence_count} unique restriction enzyme sequences")
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(self.tree_sequences)
