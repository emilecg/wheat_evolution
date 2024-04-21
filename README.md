# Missing link finder pipeline

Missing link finder pipeline is a tool developed to screen a large amount of samples sequenced at a shallow depth for detecting content of specific reference k-mers.

The pipeline consists in two main steps:
- k-mers sets preparation
- Pairwise comparison

### k-mers sets preparation
To create a specific reference k-mers set the first step is trimming the fastq files from which the k-mers needs to be counted. Trimming is a step required to reduce the amount of k-mers coming from sequencing errors, the pipeline might work as well without trimming, but the higher amount of k-mers would also raise the running time.
The trimming can be performed with any trimming tool, in our study we used Trimmomatic v. 0.38.
The next step is to generate the reference k-mers set. We used Jellyfish to generate canonical 51-mers, but any other software to count k-mers could be used. 
The final form of the reference k-mers file needs to be a single column uncompressed text file containing all the sorted 51-mers.

Since the reference file and the many sample files need to be compared using the comm bash command also the sample files need to be in the same text sorted format.

### Pairwise comparison
The comparison is performed with comm bash command counting the amount of common lines (= common k-mers) between the reference and each sample file, the value of similarity between the sets is estimated by counting the Jaccard index.
For the comparisons we had a great improvement of the performance after loading to the memory the full reference file, reading such big files (around 100 Gb) can be quite intense reading directly from the disk.
The for loop is structured to load the reference file once and after screen over multiple sample k-mers files.
To parallelise the work it's possible to make multiple lists of names of sample k-mers files and run them over multiple processes, if parallelised consider writing the output in separated .csv files to avoid multiple processes writing in the same file at the same time.


# Haplotype analysis

The subpopulation finder is a tool developed to analyse the haplotype patterning and painting a reference genome using the haplotype from wild donors. This prediction is not based on the single accessions but on the population structures of the donor, 50 Kb windows of the genomes are assigned to haplotypes based on the occurrence of the haplotype in the subpopulations defined by snmf analysis. The values of identity are determined with IBSpy. 
The script is in python (v 3.8). A single run takes around 20 minutes to be completed on a MacBook. 
Input files defined inside the script are:
- tvs file with all the variations scores of IBSpy runs against the reference genome
- A file with assign each accessions to each subpopulation (only not admixed accession included)
- Three files each one containing the lineage information for each one of the *Aegilops tauschii* samples.

The script is designed to run through 3 rounds of the population structure with K=9, K=15 and K=22. We reported and analyze only the results from K=9.
Population structure files are available in this repository. Variations scores files will be available in DRYAD after the publication of the paper.
Details about the functioning of the script can be found in the supplementary notes of the manuscript.


Preprint available at: https://doi.org/10.1101/2023.11.29.568958
