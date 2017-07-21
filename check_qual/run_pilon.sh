contigs=$1
file1=$2
file2=$3
outDir=$4
pilonBin=/home/maxim/Documents/intership/pilon/pilon-1.22.jar

mkdir -p $outDir
cd $outDir
runDir=$(pwd)

bwa index $contigs
bwa mem $contigs $file1 $file2 -t 16 -a > aln.sam

samtools view -Sb aln.sam -@ 16 > aln.bam
samtools sort aln.bam aln.sorted -@ 16
samtools index aln.sorted.bam

echo $runDir
java -Xmx16G -jar $pilonBin --genome $contigs --bam aln.sorted.bam --threads 16 --outdir $runDir