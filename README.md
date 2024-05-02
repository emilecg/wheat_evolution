Scripts used in the publication: [Origin and evolution of the bread wheat D genome](https://doi.org/10.1101/2023.11.29.568958)

# Missing Link Finder

Missing Link Finder is a pipeline developed to screen a large number of samples sequenced at a shallow depth for detecting content of specific reference k-mers.

The pipeline includes in the following steps:
- Trimming of reference raw reads
- Create reference k-mer set
- Create k-mer sets for samples to investigate (sample k-mer sets)
- Pairwise comparison between reference k-mer set and sample k-mer sets for samples
- Analysis of the output file

### Trimming of reference raw reads
To create a specific reference k-mer set, the first step consists in trimming the fastq files from which the k-mers needs to be counted. Trimming is a step required to reduce the amount of k-mers coming from sequencing errors. The pipeline might work without trimming, but the higher amount of k-mers would raise the running time.
The trimming can be performed with any trimming tool, in our study we used [Trimmomatic](https://github.com/usadellab/Trimmomatic/tree/main) v. 0.38.
```bash
# Example:
java -jar $TRIMMOMATIC_JAR PE -threads 16 -phred33 raw/accession1_1.fastq raw/accession1_2.fastq clean/accession1_1.paired.fastq.gz clean/accession1_1.unpaired.fastq.gz clean/accession1_2.paired.fastq.gz clean/accession1_2.unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:75
```
### Create reference k-mer set
The next step is to generate the reference k-mer set. We used [Jellyfish](https://github.com/gmarcais/Jellyfish) to generate canonical 51-mers, but any other software to count k-mers could be used. 
The final form of the reference k-mer file needs to be a single column uncompressed text file containing all the sorted 51-mers.
```bash
# Example:
zcat clean/reference*.paired* | jellyfish count -C -m 51 -s 50G -t 32 /dev/fd/0 -o reference_kmers.jf
jellyfish dump -L 3 -c reference_kmers.jf | sed 's/ .*//' | sort > reference_kmers.txt
```
### Create k-mers sets for samples to investigate (sample k-mer set)
Since the reference file and the many sample files need to be compared using the comm bash command, the sample files also need to be in the same text sorted format. DArTseq data used in the study were provided in FASTQCOL format. In order to prepare an input fasta file for Jellyfish, a combination of cut and awk was used.
No filtering for the occurrences of k-mers is applied in sample files due to the low coverage.
```bash
# Example:
for i in $(seq 1 80000)
do
zcat sample${i}.FASTQCOL.gz | cut -f 3 --delimiter="," | awk 'BEGIN{cont=0}{printf ">seq_%d\n",cont; print $0;cont++}' | jellyfish count -C -m 51 -s 50M -t 4 /dev/fd/0 -o sample${i}_kmers.jf
jellyfish dump -c sample${i}_kmers.jf | sed 's/ .*//' | sort > sample${i}_kmers.txt
done
```

A k-mers sorted file should look like this:
```bash
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTCTC
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCCG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCTTGC
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCTAGAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAAGAC
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGTGTTT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTAGATA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATATATCT
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAAACAC
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGTTGGG
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGTGTTTA
```

### Pairwise comparison between reference k-mer set and sample k-mer sets
The comparison is performed with comm bash command counting the amount of common lines (= common k-mers) between the reference and each sample file. The value of similarity between the sets is estimated by counting the Jaccard index. The bash command wc using -l option is used to count the k-mers in the files since each k-mer occupies a single line.
For the comparisons we achieved a good performance improvement after loading the full reference file to the memory.
The for loop is structured to load the reference file once and after screen over multiple sample k-mer files.
To parallelize the work, it is possible to make multiple lists of names of sample k-mer files and run them over multiple processes. If parallelized, consider writing the output in separated .csv files to avoid multiple processes writing in the same file at the same time.
```bash
# Example:
#counting the number of k-mers in the reference set
uno= wc -l reference/reference_kmers.txt

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
```

The list_of_accessions_names.txt file should look like this:
```bash
sample1_kmers.txt
sample2_kmers.txt
sample3_kmers.txt
sample4_kmers.txt
sample5_kmers.txt
sample6_kmers.txt
sample7_kmers.txt
sample8_kmers.txt
sample9_kmers.txt
```
### Analysis of the output file
The output file "reference_vs_accessions.csv" will have the name of the sample k-mer set in the first column, the number of k-mers in the sample k-mer set in the second column, the number of common k-mers between sample and reference k-mer sets in the third column and the computed Jaccard distance in the fourth column.
```bash
sample21443_kmers.txt,5896384,89460,2.24205E-05
sample104185_kmers.txt,5799752,97678,2.44807E-05
sample115166_kmers.txt,2128868,202,5.06E-08
sample11301_kmers.txt,8467296,72506,1.81597E-05
sample495185_kmers.txt,6903174,112303,2.81385E-05
sample519213_kmers.txt,6836817,115916,2.90443E-05
sample14108_kmers.txt,3680032,63682,1.59688E-05
sample11577_kmers.txt,8596722,126819,3.17622E-05
```


# Haplotype analysis

The subpopulation finder is a tool developed to analyse the haplotype patterning and painting a reference genome using the haplotype from wild donors. This prediction is not based on the single accessions but on the population structures of the donor. 50 kb windows of the genomes are assigned to haplotypes based on the occurrence of the haplotype in the subpopulations defined by snmf analysis. The values of identity are determined with [IBSpy](https://github.com/Uauy-Lab/IBSpy). 
The script is in Python (v 3.8). A single run takes around 20 minutes to be completed on a MacBook for a wheat genome. 
Input files defined inside the script are:
- tvs file with all the variations scores of IBSpy runs against the reference genome
- A file with assigned accessions to each subpopulation (only un-admixed accession included)
- Three files, each one containing the lineage information for each one of the three *Aegilops tauschii* samples.

The files are available in subpopulations_finder_data directory of this repository.
The script is designed to run through 3 rounds of the population structure with K=9, K=15 and K=22. We reported and analysed only the results from K=9.
Details about the functioning of the script can be found in the supplementary notes of the [manuscript](https://doi.org/10.1101/2023.11.29.568958).


