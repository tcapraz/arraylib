#!/usr/bin/env python3

import os

def writing_arg_decider(path):
    """
    whether to create or append to file.

    """
    arg = "w"
    if os.path.exists(path):
        arg = "a+"
    return arg

def txt_writer(path, outfile):
    
    """ writes the indicated outfile into a .txt file in the directory"""
    
    arg = writing_arg_decider(path)

    with open(path, arg) as output: #writes the output
        for line in outfile:
            output.write(line)
        
def file_compiler(outputed,inputed):
    
    """ compiles several files into one, deleting the temp files"""
    
    with open(outputed,"w") as out:
        for path in inputed:
            with open(path) as current:
                for line in current:
                    out.write(line)
            os.remove(path)
        out.close()