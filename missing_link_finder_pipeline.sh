###preparing the k-mers sets

#trimming only for the reference sets
java -jar $TRIMMOMATIC_JAR PE -threads 16 -phred33 raw/accession1_1.fastq raw/accession1_2.fastq clean/accession1_1.paired.fastq.gz clean/accession1_1.unpaired.fastq.gz clean/accession1_2.paired.fastq.gz clean/accession1_2.unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:75

#counting the kmers with jellyfish from paired reads
zcat clean/accession1*.paired* | jellyfish count -C -m 51 -s 50G -t 32 /dev/fd/0 -o kmers/accession1_kmers.jf

#sorting 
jellyfish dump -L 3 -c kmers/accession1_kmers.jf | sed 's/ .*//' | sort > kmers/accession1_kmers.txt

#counting the number of kmers in the reference set
uno= wc -l reference/reference_kmers.txt


###pairwise comparison
for i in $(cat lists/list_of_accessions_names.txt)
do
#counting the number of k-mers in the accession set
due=$(/usr/bin/nice -n 20 /usr/bin/ionice -c 3 wc -l kmers_sets/${i} | sed 's/ //g' | cut -f 1 -d "/")
#counting the number of common k-mers
tre=$(/usr/bin/nice -n 20 /usr/bin/ionice -c 3 comm -12 reference/reference_kmers.txt kmers_sets/${i} | /usr/bin/nice -n 20 /usr/bin/ionice -c 3 wc -l | sed 's/ //g' | cut -f 1 -d "/" )
#computing Jaccard's distance
sette=$(( $due + $uno ));otto=$(( $sette - $tre ))
distan=$(echo "scale=10;$tre/$otto" | bc -l )
#saving the result
echo ${i}","${due}","${tre}","$distan >> results/reference_vs_accessions.csv
done












