#f = open('blue_CompoundInterestNumberE.txt','r')
#for line in f:
#    if '#' not in line:
#        print line

import glob

bar = []
for filename in sorted(glob.glob('Canvas_1.txt')):
    f = open(filename,'r')
    for line in f:
        if "graph->SetPoint" in line:
            dataArray = []
            crystalEnergy = float(line.split(',',1)[1].split(',',1)[0])
            ntdEnergy = float(line.split(',',1)[1].split(',',1)[1].split(')',1)[0])
            if crystalEnergy+ntdEnergy>4.0:
                dataArray.append(crystalEnergy)
                dataArray.append(ntdEnergy)
                bar.append(dataArray)
    #iteration end

#bar = zip(*bar)

import csv
with open("output.csv", "wb") as g:
    writer = csv.writer(g)
    writer.writerows(bar)
