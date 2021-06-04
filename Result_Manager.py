# module to hold the matches that were found

class Result_Manager():
    """ This class holds the results of a particular search and is specific to a combination of RE sequences
     and DNA sequence"""

    def __init__(self, re_file, seq_file, out_file=None):
        self.re_file = re_file
        self.seq_file = seq_file
        self.out_file = out_file
        self.matches = dict()

    def print_files_used(self):
        print(f"RE file: {self.re_file}")
        print(f"Seq file: {self.seq_file}")

    def print_matches(self):
        """Function to output the matches to the file specified"""
        if self.out_file is not None:
            print(f"Printing results to file: {self.out_file}")
        else:
            print(f"Printing results to screen")
