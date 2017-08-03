from Bio import SeqIO
from collections import defaultdict

import os
import sys
sys.path.append('../SPAdes-Contig-Graph')
import spades_contig_graph


__contigsFile = "../../runs/contigs_corrector.fasta"
__snpsFile = "../../runs/contigs_corrector.used_snps"
__linksFile = "../../runs/contigs.fastg"
__pathsFile = "../../runs/contigs.paths"


#__contigsFile = "../../runs/pilon.fasta"
#__snpsFile = "../../runs/pilon.used_snps"

#__contigsFile = "../../runs/ecoli.contigs.fasta"
#__snpsFile = "../../runs/ecoli-contigs.used_snps"

class Pos:
    def __init__(self, name, pos, char):
        self.name = name
        self.pos = int(pos) - 1
        char = str(char)
        if len(char) != 1:
            raise ValueError('3nd arg is not char')
        self.char = char

    def __str__(self):
        return "(" + str(self.pos) + " " + str(self.char) + ")"


class Snp:
    def __init__(self, referencePos, contigPos):
        self.referencePos = referencePos
        self.contigPos = contigPos

    def key(self):
        return self.contigPos.name

    def __str__(self):
        return "<ref: " + str(self.referencePos) + ", contig: " + str(self.contigPos)+ ">"


class Parser:
    __refFile = "../../runs/MG1655-K12.fasta"
    reference = list(SeqIO.parse(__refFile, "fasta"))[0]
        
    @staticmethod
    def parseSnps(snpsFile):
        def toSnp(line):
            tokens = line.split()
            tokens = [x for x in tokens if x]

            (refName, contigName, refPos, refChar, contigChar, contigPos) = tokens

            return Snp(Pos(refName, refPos, refChar), Pos(contigName, contigPos, contigChar))

        def splitSnps(snps):
            maxDist = 300
                
            res = []
            lastPos = -maxDist
            for snp in snps:
                pos = snp.contigPos.pos
                if snp.contigPos.pos - maxDist >= lastPos :
                    res.append([])
                res[-1].append(snp)
                lastPos = pos
                
            return res

        
        snpsDict = defaultdict(list)
        with open(snpsFile) as f:
            snps = [toSnp(line) for line in f]
            for snp in snps:
                snpsDict[snp.key()].append(snp)

        for key, value in snpsDict.items():
            snpsDict[key] = splitSnps(sorted(value, key = lambda snp: snp.contigPos.pos))
        
        return snpsDict

    @staticmethod
    def parseContigs(contigsFile):
        return SeqIO.to_dict(SeqIO.parse(contigsFile, "fasta"))


def printPart(snps, contigsDict):
    def substr(s, getPos):
        s = str(s)
        res = [str(getPos(snps[0]).char)]
        
        for i in range(len(snps) - 1):
            cur = getPos(snps[i])    
            ne = getPos(snps[i + 1])
            res.append(s[cur.pos+1:ne.pos+1])
            if ne.char == ".":
                res.append(".")
        return ''.join(res)

    if snps and snps[0].referencePos.pos <= snps[-1].referencePos.pos:
        first = snps[0]
        contig = contigsDict[first.key()]
        last = snps[-1]

        conitgStr = []
        refStr = []
        
        sContig = substr(contig.seq, lambda snp: snp.contigPos)
        sRef = substr(Parser.reference.seq, lambda snp: snp.referencePos)

        diff = []
        cnt = 0
        for a, b in zip(list(sContig), list(sRef)):
            if a != b:
                cnt += 1
            if a == b:
                diff.append('=')
            elif a == '.' and b == '.':
                assert true
            elif a == '.':
                diff.append('i');
            elif b == '.':
                diff.append('d')
            else:
                diff.append('m')

        assert cnt == len(snps)
        assert len(sContig) == len(sRef)
        
        print(sContig)
        print(sRef)
        
        print(''.join(diff))

def printSnps(snps, contigsDict):
    for snp in snps:
        print("    " + str(snp))
    printPart(snps, contigsDict)
    print()

def getComplementNode(node_name):
    if (node_name[-1] == "'"):
        return node_name[:-1]
    else:
        return node_name + "'";

links = spades_contig_graph.load_graph_links(__linksFile)
contigs = spades_contig_graph.load_contigs(__contigsFile)
paths = spades_contig_graph.load_paths(os.path.abspath(__pathsFile), links)

for key, value in links.items():
    print(str(key) + ": " + str(value))
    print()


for contig in contigs:
    print(type(contig))
    print()

for key, value in paths.items():
    print(str(key) + ": " + str(value))
    print()

def main():
    print(len(Parser.reference.seq))
    
    contigsDict = Parser.parseContigs(__contigsFile)
    
    for key, value in contigsDict.items():
        print(str(key) + " " + str(len(value.seq)))
    print()

    snpsDict = Parser.parseSnps(__snpsFile)
    for key, value in snpsDict.items():
        print(str(key) + " {")
        for snps in value:
            printSnps(snps, contigsDict)
        print("}")
        print()
        
if __name__ == "__main__":
    pass
    #main()
