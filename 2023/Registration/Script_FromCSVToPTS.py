from importlib.resources import path
import os
import csv
import numpy as np
from fnmatch import fnmatch

rootdir = "Registration/E11.0_1"
pattern = "*.csv"
listing = []
filesep = '/'

for path, subdirs, files in os.walk(rootdir):
    for name in files:
        if fnmatch(name, pattern):
            listing.append(os.path.join(path,name))

listing.sort()
nFiles = len(listing)

for i in range(nFiles):
    os.chdir(rootdir)
    pathOrig = os.path.abspath(listing[i])
    print(pathOrig)

    with open(pathOrig) as file:
        ncols = len(file.readline().split(','))
        matrixCSV = np.loadtxt(listing[i], delimiter="," , skiprows=1, usecols=range(1,ncols))
    file.close()

    nPoints, dim = np.shape(matrixCSV)

    #save PTS
    pathDest1 = listing[i][:-3]
    pathCSV = pathDest1 + "pts"
    print(pathCSV)
    
    #os.chdir(filedir)
    with open(pathCSV, 'w') as fileID:
        fileID.write("point\n" + str(nPoints) + "\n")
        for row in matrixCSV:
            fileID.write(" ".join(str(elem) for elem in row) + '\n')
    fileID.close()

