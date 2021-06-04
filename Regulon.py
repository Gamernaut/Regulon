# Main entry point for program

import sys
import SeqTree
import Result_Manager

# used to hold the names of the 2 files
files = {}
display_option = "screen"

def display_useage_info():
    print("\ncommand line usage:  Regulon re_file=filename seq_file=filename out_file=filename -options\n")
    print("where re_file is the name of the text file containing the list of names and sequences of the restriction enzymes")
    print("where seq_file is the name of the text file containing the DNA sequence to be searched in FASTA format")
    print("and the options can be any of -h/-help (this help page) or -o which requires a filename to be provided\n")


def parse_command_line():
    opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

    # deal with the options first
    if len(opts) > 0:
        for opt in opts:
            if opt.__contains__("-h"):
                display_useage_info()
            elif opts.__contains__("-help"):
                display_useage_info()
            elif opts.__contains__("-o"):
                display_option="file"
            else:
                print(f"Unrecognised option: {opt}")
                display_useage_info()

    # now lets extract the filename parameters
    if len(args) != 2:
        print("Sorry didn't understand the format of the filenames you used.")
        display_useage_info()
        return False
    else:
        for arg in args:
            filename = arg.split("=")
            if filename[0] == "re_file":
                files["re_file"] = filename[1]
            if filename[0] == "seq_file":
                files["seq_file"] = filename[1]
            if filename[0] == "out_file" and display_option == "file":
                files["out_file"] = filename[1]
        if len(files) < 2:
            print("Sorry didn't understand the format of the filenames you used.")
            display_useage_info()
            return False
        return True


if __name__ == "__main__":
    parse_success = parse_command_line()
    if parse_success:
        print(files)  # DEBUG
        # For each run:
        # 1) Initialise a results manager which holds the results of a particular search with a specific combination of RE seqeunces and a DNA sequence
        my_result_manager = Result_Manager.Result_Manager(files["re_file"],files["seq_file"])
        my_result_manager.print_files_used()   # DEBUG
        # 2) If a tree doesn't already exist for the current RE file (based on file name) then create one (maybe update to use Hashcode to allow changes in file)
        # RE_search_tree = create_re_seq_tree(files["re_file"])
        # 3) Search the tree using a specific sequence file in FASTA format, storing any matches in the results manager
        # RE_search_tree.find_seq_matches(files["seq_file"], result_manager)
        # 4) Results manager dispays the results or saves to file depending on command line option
        print(f"Display option set to: {display_option}")
        # result_manager.print_matches()

