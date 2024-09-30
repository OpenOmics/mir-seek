#!/usr/bin/env python
from __future__ import print_function
import sys


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)


def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def index_file(file_name, join_index, header=False):
    """Converts a file to a dictionary to create an index
    for quick lookup later.
    @param file_name <str>: 
        A file or a path to a file to index.
    @param join_index <int>: 
        Index of a column to join on (zero-based).
    @param header <bool>:
        Boolean to indicate whether the file has a header.
    @returns file_dict[<str>] = list[<str>,<str>,...]: 
        A dictionary where the join index column are keys
        and the values are the entire line as a list.
    """
    filedict = {}
    fh = open(file_name, 'r')
    
    if header:
        firstline = next(fh).strip().replace('"', '').split("\t")
        filedict["_______internal-header-str"] = firstline

    for line in fh:
        linelist = line.strip().replace('"', '').split("\t")
        joinon = linelist[join_index]
        filedict[joinon] = linelist

    fh.close()
    
    return filedict


def intersect(file_dict, file2, join_index, header=False):
    """Intersects two files based on a join index.
    @param file_dict <dict[str]=list[<str>,<str>, ...]>:
        A dictionary of the first file to intersect.
        This file was indexed by the join index of the
        first file.
    @param file2 <str>:
        A path to the second file to intersect.
    @param join_index <int>:
        The index of the column to join on.
    @param header <bool>:
        Boolean to indicate whether the file has a header.
    returns None:
        Prints the intersection of the two files
        to standard output.
    """
    fh2 = open(file2, 'r')

    if header:
        firstline = next(fh2).strip().replace('"', '').split("\t")
        headerline = "\t".join(
            file_dict["_______internal-header-str"]
        ).rstrip("\n") + "\t" + "\t".join(firstline)
        print(headerline)

    for line in fh2:
        linelist = line.strip().replace('"', '').split("\t")
        joinon = linelist[join_index]
        try: file_dict[joinon]
        except KeyError: continue  # intersect not found, go to next line
        intersection = "\t".join(
            file_dict[joinon]
        ).rstrip("\n") + "\t" + "\t".join(linelist)
        try:
            print(intersection)
        except IOError:
            fh2.close()
            sys.exit()		

    fh2.close()


def main():
    
    header = False  # By default, join on column count starts at 0
    if "--header" in sys.argv:
        header = True
    
    try:
        f1_path = sys.argv[1]
        f2_path = sys.argv[2]
        f1_idx  = int(sys.argv[3])
        f2_idx  = int(sys.argv[4])
    except IndexError:
        # Bad usage, missing some of the required arguments
        err("USAGE:\nintersect <file1> <file2> <file1_col_idx> <file2_col_idx>")
        fatal("\tExample. intersect WT.bed KO.bed 0 0")

    # Create an index of the first file
    indexed_f1 = index_file(
        f1_path, 
        f1_idx, 
        header
    )

    # Intersect the second file with the
    # index of the first file, using this
    # approach we can handle very large files
    # and we only need to loop through each
    # file once.
    intersect(
        indexed_f1, 
        f2_path, 
        f2_idx, 
        header
    )

if __name__ == "__main__":
    main()
