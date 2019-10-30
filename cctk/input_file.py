import sys
import os
import re
import numpy as np
from abc import ABC, abstractmethod

class InputFile(ABC):
    """
    Represents a text file for use as input to another program like Gaussian.
    This class is abstract.
    """
    @abstractmethod
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

        if not overwrite_existing and os.path.exists(filename):
            raise ValueError(f"{filename} already exists but not allowed to overwrite")
        else:
            try:
                with open(filename, 'w+') as output_file:
                    output_file.write(text)
                return True
            except OSError as e:
                print(e)
                return False
