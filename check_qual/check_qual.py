from Bio import SeqIO
from Bio import Seq
import shutil
import copy
import random
import os
import numpy as np

alf = ['A', 'C', 'G', 'T']
__contigs_file0 = "/Johnny/data/input/Bacteria/E.coli/K12/MG1655-K12.first400K.fasta"
__reads_1="/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped/s_6.first400000_1.fastq.gz"
__reads_2="/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped/s_6.first400000_2.fastq.gz"
__reads = (__reads_1, __reads_2)

spades_dir = "../../algorithmic-biology/SPAdes-3.11.0-dev"
pwd = os.path.dirname(__file__)
mismatch_rate = 0.3
indel_rate = 0.3
contigs_file = pwd + "/contigs.fasta"

broken_dir = pwd + "/broken"
if not os.path.exists(broken_dir):
    os.makedirs(broken_dir)
contigs_broken_file = broken_dir + "/contigs.fasta"


class ErrorGenerator:
    def __init__(self, mismatch_rate, insert_rate, delete_rate):
        self.mismatch_rate = mismatch_rate
        self.insert_rate = insert_rate
        self.delele_rate = delete_rate
        self.ok_rate = 1 - insert_rate - mismatch_rate - delete_rate

    def make_mismatch(self, char):
        return random.choice([new_char for new_char in alf if new_char != char])

    def gen_error(self, char):
        indel_str = char + random.choice(alf)
        mismatch_str = self.make_mismatch(char)
        return np.random.choice([indel_str, '', mismatch_str, char],
                                p=[self.insert_rate, self.delele_rate, self.mismatch_rate, self.ok_rate])

    def gen_errors(self, chars):
        return list(''.join(map(self.gen_error, chars)))

    def gen_errors_in_record(self, record, f):
        record = copy.deepcopy(record)
        chars = list(str(record.seq))
        chars = f(chars)
        chars = ''.join(chars)
        record.seq = Seq.Seq(chars)
        return record

    def gen_errors_in_contig(self, contig):
        return self.gen_errors_in_record(contig, self.gen_errors)

def gen_mismatches():
    gen = ErrorGenerator(mismatch_rate, indel_rate / 2, indel_rate / 2)
    contigs = list(SeqIO.parse(contigs_file, "fasta"))
    contigs_broken = list(map(gen.gen_errors_in_contig, contigs))
    SeqIO.write(contigs_broken, contigs_broken_file, "fasta")

class Corrector:
    def __init__(self, contigs_path, reads, outdir, corrected_path):
        self.contigs_path = contigs_path
        self.outdir = outdir
        self.corrected_path = corrected_path
        self.reads_1, self.reads_2 = reads

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def get_corrected_path(self):
        return self.corrected_path

    def show_snps(self, reference_path):
        __nucmer_path = "/home/maxim/Documents/intership/MUMmer3.23"

        out_file = self.outdir + "/contigs_mismatch_{}_indel_{}".format(mismatch_rate, indel_rate)
        nucmer_cmd = "{}/nucmer --prefix={} {} {}".format(__nucmer_path, out_file, reference_path, self.corrected_path)
        snps_cmd = "{}/show-snps -Clr {}.delta > {}.snps".format(__nucmer_path, out_file, out_file)
        os.system(nucmer_cmd)
        os.system(snps_cmd)

class Pilon(Corrector):
    def __init__(self, contigs_path, reads, outdir):
        Corrector.__init__(self, contigs_path, reads, outdir, outdir + "/pilon.fasta")

    def run(self):
        run_pilon_cmd = "/bin/bash run_pilon.sh {} {} {} {}".format(self.contigs_path, self.reads_1, self.reads_2, self.outdir)
        os.system(run_pilon_cmd)

class SpadesCorrector(Corrector):
    def __init__(self, contigs_path, reads, outdir):
        Corrector.__init__(self, contigs_path, reads, outdir, outdir + "/corrected_contigs.fasta")
        self._corrector_tmp = self.outdir + "/tmp"

        yml_file = pwd + "/input_dataset.yaml"
        input_dataset = "- left reads: [{}]\n".format(self.reads_1) + \
                        "  orientation: fr\n" + \
                        "  right reads: [{}]\n".format(self.reads_2) + \
                        "  type: paired-end\n"
        with open(yml_file, "w") as f:
            f.write(input_dataset)

        if not os.path.exists(self._corrector_tmp):
            os.makedirs(self._corrector_tmp)

        self._corrector_config = self.outdir + "/config.info"

        with open(self._corrector_config, 'w') as config:
            config.write('"bwa": "{}/bin/bwa-spades"\n'.format(spades_dir))
            config.write('"dataset": "{}"\n'.format(yml_file))
            config.write('"max_nthreads": !!int "16"\n')
            config.write('"output_dir": "{}"\n'.format(self.outdir))
            config.write('"strategy": "mapped_squared"\n')
            config.write('"work_dir": "{}"\n'.format(self._corrector_tmp))

    def run(self):
        cmd_run = "{}/bin/corrector {} {}".format(spades_dir, self._corrector_config, self.contigs_path)
        print(cmd_run)
        os.system(cmd_run)

        shutil.rmtree(self._corrector_tmp)

def main():
    gen_mismatches()
    corrector1 = SpadesCorrector(contigs_broken_file, __reads, pwd + "/corrector1")
    corrector1.run()
    corrector1.show_snps(contigs_file)

    corrector2 = SpadesCorrector(corrector1.get_corrected_path(), __reads, pwd + "/corrector2")
    corrector2.run()
    corrector2.show_snps(contigs_file)

    pilon1 = Pilon(contigs_broken_file, __reads, pwd + "/pilon1")
    pilon1.run()
    pilon1.show_snps(contigs_file)

    pilon2 = Pilon(pilon1.get_corrected_path(), __reads, pwd + "/pilon2")
    pilon2.run()
    pilon2.show_snps(contigs_file)

    #correctorBroken = Corrector(contigs_broken_file, __reads, broken_dir, contigs_broken_file)
    #correctorBroken.show_snps(contigs_file)
#main()

def test():
    chars = alf[:] + alf[:]
    def subtest(chars, f):
        chars2 = chars[:]
        chars2 = f(chars2)
        print(chars2)

    for i in range(2):
        gen = ErrorGenerator(mismatch_rate, indel_rate / 2, indel_rate / 2)
        subtest(chars, gen.gen_errors)

    print()

test()
