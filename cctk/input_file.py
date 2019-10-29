import sys
import os 
import re
import numpy as np

class InputFile():
    def __init__(self):
        pass    

    def write_file(self, filename,text, overwrite_existing=True):
        """
        Writes output text to a file.

        Args:
            filename (str): path to file, including name (e.g. "path/to/input.gjf")
            text (str): desired contents of file
            overwrite_existing (Bool): whether any existing files should be overwritten or not
        
        Returns:
            `True` if write succeeded, `False` otherwise
        """
        if not isinstance(text, str):
            raise TypeError("cannot write non-string to file!")
            return False

        if os.path.exists(filename) and overwrite_existing==False:
            raise ValueError("file already exists; overwrite_existing set to False!")
            return False
        else:
            try:   
                with open(filename, 'w+') as output_file:
                    output_file.write(text) 
                return True
            except OSError:
                return False
