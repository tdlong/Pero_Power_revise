import sys
import re

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

myfolder=sys.argv[1]
part=sys.argv[2]

tfin= myfolder + "/temppow_LOD_" + part + ".txt"
tfout= myfolder + "/temppow_pow_" + part + ".txt"
fin = open(tfin, 'r')
fout = open(tfout, 'w')
fnamesin = open("colnames.yr.348.small.txt", 'r')
colnames = fnamesin.readline()
colnames = colnames.rstrip().replace('"','').split(" ")

ncol = len(colnames)
Cmax =  [0] * ncol
Cchr = ['NA'] * ncol
Cpos = [0] * ncol
Chit = [0] * ncol

for line in fin:
	row = line.rstrip().split(" ")
	chr=row[0]
	pos=row[1]
	for i in range(2, ncol+2):
		lod = row[i]
		if isfloat(lod) and float(lod) > 7.5:
			Chit[i-2] = Chit[i-2] + 1
		if isfloat(lod) and float(lod) > Cmax[i-2]:
			Cmax[i-2] = float(lod)
			Cchr[i-2] = chr
			Cpos[i-2] = pos

# now print out for each phenotype
for i in range(0, ncol):
	fout.write(colnames[i]+"\t"+Cchr[i]+"\t"+str(Cpos[i])+"\t"+str(Cmax[i])+"\t"+str(Chit[i])+"\n")


