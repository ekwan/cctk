import sys
import re
import numpy as np

from cctk import InputFile

class GaussianInputFile(InputFile):
    def __init__(self):
        pass    

    def write_file(filename,text):
        super.write_file()
 
