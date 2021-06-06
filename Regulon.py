# Main entry point for program

import sys
import SeqTree
import ResultManager

# used to hold the names of the 2 files and set the default output mode to display on the screen
files = {"re_file" : "" , "seq_file" : "" , "out_file" : "screen"}

def display_useage_info():
    print("\ncommand line usage:  Regulon re_file=filename seq_file=filename [out_file=filename] -options\n")
    print("Mandatory: re_file is the name of the text file containing the list of names and sequences of the restriction enzymes")
    print("Mandatory: seq_file is the name of the text file containing the DNA sequence to be searched in FASTA format")
    print("Optional:  out_file is the name of results file, if not used the rersults will be diplayed on the screen")
    print("and the options can be any of -h/-help (this help page)\n")


def parse_command_line():
    global display_mode
    opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
    args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

    # deal with the options first
    if len(opts) > 0:
        for opt in opts:
            if opt.__contains__("-h"):
                display_useage_info()
            elif opts.__contains__("-help"):
                display_useage_info()
            else:
                print(f"Unrecognised option: {opt}")
                display_useage_info()

    # now lets extract the filename parameters
    if len(args) < 2:
        print("Sorry didn't understand the format of the command line arguments you used.")
        display_useage_info()
        return False
    else:
        for arg in args:
            filename = arg.split("=")
            if filename[0] == "re_file":
                files["re_file"] = filename[1]
            if filename[0] == "seq_file":
                files["seq_file"] = filename[1]
            if filename[0] == "out_file":
                files["out_file"] = filename[1] # Output file was specificied on command line so update dictionary

        # Check to make sure the required filenames have been added to the dictionary
        if not files["re_file"] or not files["seq_file"]:
            print("\nSorry didn't understand the format of the filenames you used.")
            display_useage_info()
            return False
        return True


if __name__ == "__main__":
    parse_success = parse_command_line()
    if parse_success:
        print(files)  # DEBUG
        # For each run:
        # 1) Initialise a results manager which holds the results of a particular search with a specific combination of RE seqeunces and a DNA sequence
        result_manager = ResultManager.Result_Manager(files["re_file"], files["seq_file"], files["out_file"])
        # my_result_manager.print_files_used()   # DEBUG

        # 2) Step that requests the root of a newly created tree if one doesn't already exist for the current RE
        # file (based on file name or maybe Hashcode in future version).
        # If one does exist it just gets the root of that tree.
        current_search_tree = SeqTree.RESeqTree()
        if not current_search_tree:
            print("Couldn't create search tree")
            exit(-1)

        # We have a valid root so build the tree with the sequences in the Restriction Enzyme Definition file
        current_search_tree.build_tree(files["re_file"])

        # 3) Search the tree using a specific sequence file in FASTA format, storing any matches in the results manager
        current_search_tree.find_matches(files["seq_file"], result_manager)

        # 4) Results manager dispays the results or saves to file depending on command line arguments used
        # display_mode = files["out_file"] # DEBUG
        # print(f"Display option set to: {display_mode}") #DEBUG
        result_manager.print_matches()
