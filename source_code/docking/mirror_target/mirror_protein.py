import argparse
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

# Set up argument parser
parser = argparse.ArgumentParser(description='Mirror the origin protein and resn unchanged')
parser.add_argument('-i', '--input', required=True, help='Input PDB file')
parser.add_argument('-o', '--output', required=True, help='Output PDB file')
args = parser.parse_args()
print(args)

target = args.input
output = args.output
p = PDBParser()
s = p.get_structure(target[:-4], target)

for chains in s:
    for chain in chains:
        for residue in chain:
            for atom in residue:
                coordinate = atom.get_vector()
                x = coordinate[0]
                y = coordinate[1]
                z = -coordinate[2]
                newcoord = x, y, z
                atom.set_coord(newcoord)

io = PDBIO()
io.set_structure(s)
io.save(output)