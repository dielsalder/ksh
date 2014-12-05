import re
import sys

inp_pdb = '3MVA.pdb'
#inp_pdb = sys.argv[1]

crds = []
with open(inp_pdb) as pdb:
    for line in pdb:
        if re.search('(^ATOM)\s*\S*\s*\S*\s*(DA5|DA3|DA|DT5|DT3|DT|DG5|DG3|DG|DC5|DC3|DC)', line):
            row = line.split()
            crds.append(row[2:4] + row[6:9])

with open('nucleic.dat', 'w') as out:
    for row in crds:
        out.write(('{} {} {} {} {}\n').format(row[0], row[1], row[2], row[3], row[4]))
