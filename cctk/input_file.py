import sys
import re
import numpy as np

class InputFile():
    def __init__(self):
        pass    

    def write_file(self, filename,text):
        """
        Writes output text to a file.

        Args:
            filename (str): path to file, including name (e.g. "path/to/input.gjf")
            text (str): desired contents of file
        """
        with open(filename, 'w+') as output_file:
            output_file.write(text) 
        print(text)

 
