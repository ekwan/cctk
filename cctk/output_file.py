import sys
import re
import numpy as np

class OutputFile():
    def __init__(self):
        pass    

    def read_file(self, filename):
        '''
        Read a file.
        '''
        with open(filename, 'r') as filehandle:
            lines = filehandle.read().splitlines()
            return lines
             
