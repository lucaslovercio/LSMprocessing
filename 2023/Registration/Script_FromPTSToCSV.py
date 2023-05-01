import os
import numpy as np
from fnmatch import fnmatch
import csv

rootdir = "Registration/E11.0_1"
pattern = "*.pts"
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

    with open(pathOrig) as file:
        matrixPTS = np.loadtxt(listing[i], delimiter=' ', skiprows=2)
    file.close()

    #save CSV
    pathDest1 = listing[i][:-3]
    pathCSV = pathDest1 + "csv"
    print(pathCSV)

    with open(pathCSV, 'w') as csvfile:
        csvfile.write(" ,X,Y,Z\n")
        count = 1
        for row in matrixPTS:
            csvfile.write(str(count) + ',' + ','.join(str(elem) for elem in row) + '\n')
            count += 1
    csvfile.close()
    

