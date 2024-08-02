import os
import argparse

flip_dict = {"GLY":"GLY","DAL":"ALA","DVA":"VAL","DIL":"ILE","DLE":"LEU",
                "DSE":"SER","DTH":"THR","DHI":"HIS","DLY":"LYS","DAR":"ARG",
                "DTR":"TRP","DPH":"PHE","DTY":"TYR","DME":"MET","DCY":"CYS",
                "DGU":"GLU","DGN":"GLN","DPR":"PRO","DAS":"ASP","DAN":"ASN",
                "ALA":"DAL","VAL":"DVA","ILE":"DIL","LEU":"DLE","SER":"DSE",
                "THR":"DTH","HIS":"DHI","LYS":"DLY","ARG":"DAR","TRP":"DTR",
                "PHE":"DPH","TYR":"DTY","MET":"DME","CYS":"DCY","GLU":"DGU",
                "GLN":"DGN","PRO":"DPR","ASP":"DAS","ASN":"DAN","DCS":"CYS"}

def changetoD(pdb):
    with open(pdb, 'r') as pdbfile:
        lines = pdbfile.readlines()
    data = ''
    for line in lines:
        line = list(line)
        if len(line) > 21 and line[21] == 'B':
            line[17:20] = list(flip_dict[''.join(line[17:20])])
        data += ''.join(line)
    with open(f'D_{pdb}', 'w') as pdbfile:
        pdbfile.write(data)

def main():
    parser = argparse.ArgumentParser(description='Process PDB files.')
    parser.add_argument('-i', '--input', nargs='+', help='Input PDB file(s)', required=True)
    args = parser.parse_args()
    pdb = args.input
    changetoD(pdb)

if __name__ == "__main__":
    main()