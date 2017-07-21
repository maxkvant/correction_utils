path=/Johnny/data/input/Bacteria/E.coli/K12/is220/cropped
len=400000
fileName1=s_6.first${len}_1.fastq.gz
fileName2=s_6.first${len}_2.fastq.gz

file1=$path/$fileName1
file2=$path/$fileName2

bin=/home/maxim/Documents/intership/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar

echo $file1
echo $file2

java -jar $bin PE -phred33 $file1 $file2 -baseout $path/trimmed_s_6.first${len}.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
