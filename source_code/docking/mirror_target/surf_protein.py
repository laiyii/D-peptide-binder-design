import argparse

def blank(s):
    s = s.split(' ')
    ss = []
    for i in s:
        if len(i) != 0:
            ss.append(i)
    return ss

parser = argparse.ArgumentParser(description='Process some files.')
parser.add_argument('-i_asa', type=str, required=True, help='Input ASA file')
parser.add_argument('-i', type=str, required=True, help='Input PDB file')
parser.add_argument('-o', type=str, required=True, help='Output PDB file')

args = parser.parse_args()

asa = open(args.i_asa, 'r')
asas = asa.readlines()
mono = open(args.i, 'r')
monoline = mono.readlines()
new_pdb = open(args.o, 'w')

l = ''
for i in range(len(asas)):
    area = float(blank(asas[i])[-2])
    monol = list(monoline[i])
    if area > 1:
        monol[69] = '1'
    else:
        monol[69] = '0'
    for j in range(70, len(monol)):
        monol[j] = ' '
    for j in range(56, 66):
        monol[j] = ' '
    monol = ''.join(monol)
    monol = monol + '\n'
    l += monol

new_pdb.write(l)