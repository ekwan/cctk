import sys
import re
import numpy as np

class OutputFile():
    '''
    Class for generic computational chemistry output files. 

    Contains generic methods used for reading and parsing all output files.
    '''
    
    def __init__(self):
        pass    

    def read_file(self, filename):
        '''
        Reads a file and outputs a list of the lines. .
        
        Args:
            filename (str): The path to the file. 

        Returns:
            A list containing all the lines in the file.
        '''
        with open(filename, 'r') as filehandle:
            lines = filehandle.read().splitlines()
            return lines
             
