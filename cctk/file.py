import sys
import os
import re
import numpy as np

from abc import ABC, abstractmethod


class File(ABC):
    """
    Represents a text file for use as input or output to another program like Gaussian.
    This class is abstract.
    """

    @abstractmethod
    def __init__(self):
        pass

    @staticmethod
    def write_file(filename, text, overwrite_existing=True):
        """
        Writes output text to a file.

        Args:
            filename (str): path to file, including name (e.g. ``path/to/input.gjf``)
            text (str): desired contents of file
            overwrite_existing (Bool): whether any existing files should be overwritten or not

        Returns:
            ``True`` if write succeeded, ``False`` otherwise
        """
        if not isinstance(text, str):
            raise TypeError("cannot write non-string to file!")
            return False

        if not overwrite_existing and os.path.exists(filename):
            raise ValueError(f"{filename} already exists but not allowed to overwrite")
        else:
            try:
                with open(filename, "w+") as output_file:
                    output_file.write(text)
                return True
            except OSError as e:
                print(e)
                return False

    @staticmethod
    def append_to_file(filename, text):
        """
        Appends output text to a file.

        Args:
            filename (str): path to file, including name (e.g. ``path/to/input.gjf``)
            text (str): desired contents of file

        Returns:
            ``True`` if write succeeded, ``False`` otherwise
        """
        if not isinstance(text, str):
            raise TypeError("cannot write non-string to file!")
            return False

        if os.path.exists(filename):
            try:
                with open(filename, "a+") as output_file:
                    output_file.write(text)
                return True
            except OSError as e:
                print(e)
                return False
        else:
            raise ValueError(f"{filename} does not exist")

    @staticmethod
    def read_file(filename):
        """
        Reads a file and parses into lines.

        Args:
            filename (str): The path to the file.

        Returns:
            A list containing all the lines in the file.
        """
        with open(filename, "r") as filehandle:
            lines = filehandle.read().splitlines()
            return lines
