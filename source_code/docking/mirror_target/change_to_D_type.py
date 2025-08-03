import argparse
import os

flip_dict = {"GLY":"GLY","DAL":"ALA","DVA":"VAL","DIL":"ILE","DLE":"LEU",
                "DSE":"SER","DTH":"THR","DHI":"HIS","DLY":"LYS","DAR":"ARG",
                "DTR":"TRP","DPH":"PHE","DTY":"TYR","DME":"MET","DCY":"CYS",
                "DGU":"GLU","DGN":"GLN","DPR":"PRO","DAS":"ASP","DAN":"ASN",
                "ALA":"DAL","VAL":"DVA","ILE":"DIL","LEU":"DLE","SER":"DSE",
                "THR":"DTH","HIS":"DHI","LYS":"DLY","ARG":"DAR","TRP":"DTR",
                "PHE":"DPH","TYR":"DTY","MET":"DME","CYS":"DCY","GLU":"DGU",
                "GLN":"DGN","PRO":"DPR","ASP":"DAS","ASN":"DAN"}

def changetoD(input_file, output_file):
    with open(input_file, 'r') as pdbfile:
        lines = pdbfile.readlines()
    
    data = ''
    for line in lines:
        line = list(line)
        if line[:4] == list('ATOM'):
            line[17:20] = list(flip_dict[''.join(line[17:20])])
        if line[:3] == list('TER'):
            line[17:20] = list(flip_dict[''.join(line[17:20])])
        data += ''.join(line)
    
    with open(output_file, 'w') as pdbfile:
        pdbfile.write(data)

# Set up argument parser
parser = argparse.ArgumentParser(description='Change protein residues to D-type')
parser.add_argument('-i', '--input', required=True, help='Input PDB file')
parser.add_argument('-o', '--output', required=True, help='Output PDB file')
args = parser.parse_args()

changetoD(args.input, args.output)