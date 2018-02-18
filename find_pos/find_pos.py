import re
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import subprocess

len = 600
file_from = "/media/maxim/DATA/intership/output/SRR5713995/ref/contigs.fasta"
contig_from_name = "CP011636.1"
pos_from = 3 * (10 ** 6)

file_where = "/media/maxim/DATA/intership/output/SRR5713995/canu/SRR5713995.contigs_fixed_names.fasta"
tmp_file = "tmp.fasta"

def make_tmp():
    contigs_dict = SeqIO.to_dict(SeqIO.parse(file_from, "fasta"))
    contig_from = contigs_dict[contig_from_name]
    fragment = contig_from[pos_from : pos_from + len]
    with open(tmp_file, "w") as output_handle:
        SeqIO.write([fragment], output_handle, "fasta")

def run_bwa():
    subprocess.call(["bwa", "index", file_where])
    output = subprocess.check_output(["bwa", "mem", file_where, tmp_file]).strip()
    last_line = output.splitlines()[-1].decode("utf-8")
    read, flag, contig, pos = last_line.split()[:4]
    print(pos)

if __name__ == "__main__":
    make_tmp()
    run_bwa()