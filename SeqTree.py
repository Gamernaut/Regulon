# Module to create and search the tree based data structure of restriction enzyme sequences

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
        self.is_leaf = None


class RE_seq_tree():
    """This is the class that represents the restriction enzyme (RE) sequence tree which is made up of the Node class
    above """

    def __init__(self):
        self.root = Node()
        self.tree_width = 0
        self.tree_depth = 0
        self.DNA_sequence = ""
        self.sequence_count = 0
        self.re_filename = ""   # If a new request to create with same file then we can just pass the root back (singelton type idea)
        self.re_file_hashcode = 0  # More robust test in case the content of the file has changed but not the name, the hashcode will pick up the difference.
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



if __name__ == '__main__':
    root = Node()
    insert_sequence(root, "TACGTT", "AclI")
    insert_sequence(root, "AAGCTT", "HindIII")
    insert_sequence(root, "AATATT", "SspI")
    insert_sequence(root, "AATY", "MluCI_Y")
    insert_sequence(root, "AATR", "MluCI_R")
    insert_sequence(root, "AAGCTT", "HindIII V2")
    print_tree(root)
