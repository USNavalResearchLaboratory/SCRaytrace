#Imports
import pandas as pd
import numpy as np
import os
from array import array
import sys


#Function to write the data to a dat file in the form: number of items... items... number of items... items...
def write_to_dat(data, nd): #nd is the number of dimensions of data
    path = os.path.dirname(os.path.dirname(os.path.dirname(__file__))) + "/data/VSFvariable.dat"
    file = open(path, "wb")
    # print(np.int32(len(data)).tobytes())
    file.write(np.int32(len(data[0])).tobytes())
    for i in range(nd):
        file.write(bytes(data[i]))
    file.close()


#VSF Functions
def testfunc():
    data = []
    data.append(np.array([1, 2, 3, 4, 5], dtype='f4'))
    write_to_dat(data, 1)



#Main for testing
if __name__ == "__main__":
    # testfunc()
    pass