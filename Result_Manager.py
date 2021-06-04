# module to hold the matches that were found

class Result_Manager():
    """ This class holds the results of a particular search and is specific to a combination of RE sequences
     and DNA sequence"""

    def __init__(self, re_file, seq_file):
        self.re_file = re_file
        self.seq_file = seq_file
        self.matches = dict()

    def print_files_used(self):
        print(f"RE file: {self.re_file}")
        print(f"Seq file: {self.seq_file}")

    def print_matches(self):
        """Function to either print on screen or output to file the results of the matches found"""
