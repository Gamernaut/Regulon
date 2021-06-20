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
        # Using a hashcode to check the content of the file is a much better way to ensure that the file
        # contents have not been changed even if the name stays the same

        if filename == self.re_filename:
            return self.get_root()

        # It's a new file so we need a new tree, set the name so the check above will catch the reuse of the file
        self.re_filename = filename

        # Call the function to load the sequences from the restriction enzyme definition file
        self.re_seq_dict, self.tree_depth = FileHandler.import_restriction_enzymes(filename)

        # Calculate the progress of the tree build and display this to the user
        current_seq_count = 1
        sequences_to_build = len(self.re_seq_dict)
        for re_name, seq in self.re_seq_dict.items():
            pct_complete = (current_seq_count / sequences_to_build) * 100
            sys.stdout.write(f"\rProcessing Restriction Enzymes sequences from {filename}: %d%% completed"
                             % pct_complete)  # Prints on the same line
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

    def find_matches(self, filename, result_manager):
        dna_sequence = FileHandler.import_seq_file(filename)
        dna_seq_length = len(dna_sequence)
        print(f"DNA sequence is {dna_seq_length} bases long")
        ref_seq_length = self.get_tree_depth()
        print(f"Longest reference sequence is {ref_seq_length} bases long")
        node = self.get_root()

        # To optimise search limit the number of nucleotides beyond the current sequence read to the length of the
        # longest reference sequence in the tree i.e tree_depth.
        for pos in range(dna_seq_length):
            if pos + ref_seq_length > dna_seq_length:
                search_window_end = dna_seq_length
            else:
                search_window_end = pos + ref_seq_length
            dna_sub_sequence = dna_sequence[pos:search_window_end]
            if dna_sub_sequence[0] == "A" and node.A:
                self.search_branch(node.A, dna_sub_sequence, pos)
            if dna_sub_sequence[0] == "C" and node.C:
                self.search_branch(node.C, dna_sub_sequence, pos)
            if dna_sub_sequence[0] == "G" and node.G:
                self.search_branch(node.G, dna_sub_sequence, pos)
            if dna_sub_sequence[0] == "T" and node.T:
                self.search_branch(node.T, dna_sub_sequence, pos)

    def search_branch(self, node, dna_sub_sequence, position):
        if node is None:
            print("No valid tree node provided to SeqTree.search_branch()")
            return
        if node.is_branch_end:
            # TODO: make it work with ResultManager
            # Add the current
            #if position in self.matches:
            #    print("Duplicate position match at ", position)
            #    print("Dictionary contains", self.matches[position])
            #    print("trying to add", node.RE_list)
            #    # get the existing list for this key
            #    re_list = self.matches[position]
            #    # Use extend to combine the 2 lists but keep them flat. Append would nest one list inside the other.
            #    re_list.extend(node.RE_list)
            #    # convert to set to remove duplicates the back to list
            #    re_list = list(set(re_list))
            #    # overwrite existing entry with new list
            #    self.matches[position] = re_list
            #    print("Dict now reads ->", self.matches[position])
            #else:
            # Add the names of the reference sequences and the original position to the results manager dictionary
            # self.matches[position] = node.RE_list

            # need to add 1 to position because loop starts from 0
            print(f"found match at {position + 1} for {node.RE_list}")
            return
        dna_sub_sequence = dna_sub_sequence[1:]
        if dna_sub_sequence[0] == "A" and node.A:
            self.search_branch(node.A, dna_sub_sequence, position)
            return
        if dna_sub_sequence[0] == "C" and node.C:
            self.search_branch(node.C, dna_sub_sequence, position)
            return
        if dna_sub_sequence[0] == "G" and node.G:
            self.search_branch(node.G, dna_sub_sequence, position)
            return
        if dna_sub_sequence[0] == "T" and node.T:
            self.search_branch(node.T, dna_sub_sequence, position)
            return

    def print_branch(self, node):
        if node is None:
            print("No valid tree node provided to SeqTree.print_branch()")
            return
        # Add the nucleotide in the current node to the sequence constructed so far (from previous calls)
        # If it's the root there is no nucleotide so skip this
        if node.nucleotide:
            self.branch_sequence = self.branch_sequence + node.nucleotide
        if node.is_branch_end:
            temp = self.branch_sequence + " -> " + ', '.join(node.RE_list)  # replaces [] from std List print with ","
            self.tree_sequences.append(temp)
        if node.A:
            self.print_branch(node.A)
            # as you come back up the call stack remove this nucleotide
            self.branch_sequence = self.branch_sequence[:-1]
        if node.C:
            self.print_branch(node.C)
            # as you come back up the call stack remove this nucleotide
            self.branch_sequence = self.branch_sequence[:-1]
        if node.G:
            self.print_branch(node.G)
            # as you come back up the call stack remove this nucleotide
            self.branch_sequence = self.branch_sequence[:-1]
        if node.T:
            self.print_branch(node.T)
            # as you come back up the call stack remove this nucleotide
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
