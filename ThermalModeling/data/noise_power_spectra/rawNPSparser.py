import numpy as np

allData = ""
dataset = str(2148)

f = open('raw/ds'+dataset,'r')

for line in f:
    allData += line
allData = allData.replace('\n', '')
allDataList = np.transpose(map(float, allData[1:-2].split(',')))

import csv

with open('ds'+dataset+'.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(allDataList)
