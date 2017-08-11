from collections import defaultdict
import os

__contigsFile = "../../runs/contigs_corrector.fasta"
__snpsFile = "../../runs/contigs_corrector.used_snps"
__linksFile = "../../runs/contigs.fastg"
__pathsFile = "../../runs/contigs.paths"


#__contigsFile = "../../runs/pilon.fasta"
__snpsFile_before = "/home/m_vinnichenko/intership/output/R.sphaeroides/assembled_contigs.used_snps"
__snpsFile_after = "/home/m_vinnichenko/intership/output/R.sphaeroides/corrected_contigs.used_snps"

__logFile = "/home/m_vinnichenko/intership/output/R.sphaeroides/run27.log"
__snpsDir = "/home/m_vinnichenko/intership/output/R.sphaeroides"


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

class CorrectorLog:
    class Node:
        def __init__(self, interesting, changed, coverage):
            self.interesting = interesting
            self.changed = changed
            self.coverage = coverage

        def __hash__(self):
            return hash((self.interesting, self.changed))

        def __eq__(self, other):
            return (self.interesting, self.changed) == (other.interesting, other.changed)

        def __ne__(self, other):
            return not (self == other)

        def __str__(self):
            return "(interesting {}, changed {})".format(self.interesting, self.changed)

    def __init__(self):
        self.changed = defaultdict(bool)
        self.interesting = defaultdict(bool)
        self.coverage = defaultdict(int)

    def get(self, key):
        if key in self.coverage:
            return CorrectorLog.Node(self.interesting[key], self.changed[key], self.coverage[key])
        return None

class Parser:
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
            snpsDict[key] = sorted(value, key = lambda snp: snp.contigPos.pos)
        return snpsDict

    @staticmethod
    def simpleName(name):
        if not "NODE_" in name:
            return name
        (pref, id) = name.split("_")[:2]
        return "_".join([pref, id])

    @staticmethod
    def toKey(contig, pos):
        return (Parser.simpleName(contig), int(pos))

    @staticmethod
    def parseLog(logFile):
        interesting_tag = "inited_interesting_pos"
        changed_tag = "change_pos"
        tags = [interesting_tag, changed_tag]
        low_cov = 20

        def checkLine(line):
            if (not "INFO" in line) and (not "DEBUG" in line):
                return False
            for tag in tags:
                if tag in line:
                    return len(line.split(")   ")) == 2
            return False

        res = CorrectorLog()
        with open(logFile) as f:
            for line in f:
                if not checkLine(line):
                    continue
                _, line = line.split(")   ")
                line = line.strip()
                tag = [tag for tag in tags if tag in line][0]

                if (tag == interesting_tag):
                    (_, pos, _, coverage, _, contig) = line.split()
                    key = Parser.toKey(contig, pos)
                    res.interesting[key] = True

                if (tag == changed_tag):
                    (tag, pos, _, _, _, _, coverage, _, contig) = line.split()
                    key = Parser.toKey(contig, pos)
                    res.changed[key] = True


                res.coverage[key] = int(coverage)

                #print(tag, key, res.coverage[key], res.interesting[key], res.changed[key], end='\n')
        return res

def printSnps(snps):
    for snp in snps:
        print("    " + str(snp))

def getLogStats(snpsFile, logFile):
    snpsDict = Parser.parseSnps(snpsFile)
    correctorLog = Parser.parseLog(logFile)
    logStats = defaultdict(set)

    for (contig, snps) in snpsDict.items():
        for snp in snps:
            pos = snp.contigPos.pos
            key = Parser.toKey(contig, pos)
            logStats[correctorLog.get(key)].add(key)
            #if (key in correctorLog.coverage):
            #    print("{} interesting:{} changed:{} coverage:{}".format(key, correctorLog.interesting[key], correctorLog.changed[key], correctorLog.coverage[key]))
            #else:
            #    print("skipped {}".format(key))
    return logStats

def printLogStats(logStats):
    size = sum([len(vals) for (_, vals) in logStats.items()])
    for (key, vals) in logStats.items():
        print("{1} {0}".format(key, len(vals) / size))
#        for val in vals:
#            print("   " + str(val))
#        print()

def defaultdictSetMinus(dict1, dict2):
    res = defaultdict(set)
    for key in dict1:
        res[key] = dict1[key] - dict2[key]
    for key in dict2:
        res[key] = dict1[key] - dict2[key]
    return res

def main():
    log_stats_before = getLogStats(__snpsFile_before, __logFile)
    for __snpsFile_after in os.listdir(__snpsDir):
        if not __snpsFile_after.endswith(".used_snps"):
            continue

        __snpsFile_after = __snpsDir + "/" + __snpsFile_after
        log_stats_after = getLogStats(__snpsFile_after, __logFile)

        print("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+= " + __snpsFile_after)

        print("========== BEFORE ==========")
        printLogStats(log_stats_before)

        print("========== AFTER ==========")
        printLogStats(log_stats_after)

        print("========== BEFORE - AFTER ==========")
        printLogStats(defaultdictSetMinus(log_stats_before, log_stats_after))

        print("========== AFTER - BEFORE ==========")
        printLogStats(defaultdictSetMinus(log_stats_before, log_stats_after))
        print()
        print()

if __name__ == "__main__":
    pass
    main()
