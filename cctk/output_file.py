import sys
import re
import numpy as np
from abc import ABC, abstractmethod

class OutputFile(ABC):
    '''
    Represents an output file from another program like Gaussian.
    This class is abstract.
    Class for generic computational chemistry output files.

    Contains generic methods used for reading and parsing all output files.
    '''

    @abstractmethod
    def __init__(self):
        pass

    def read_file(self, filename):
        '''
        Reads a file and parses into lines.

        Args:
            filename (str): The path to the file.

        Returns:
            A list containing all the lines in the file.
        '''
        with open(filename, 'r') as filehandle:
            lines = filehandle.read().splitlines()
            return lines
