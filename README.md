# Rbec: a tool for analysis of amplicon sequencing data from synthetic microbial communities

Our manuscript is currently available as a preprint in [*bioRxiv*](https://doi.org/10.1101/2021.01.15.426834).

Rbec is a tool for analysing amplicon sequencing data from synthetic communities (SynComs), where the reference sequences for each strain are already available. Rbec can accurately correct PCR and sequencing errors, identify intra-species polymorphic variation, and detect contaminations in SynCom amplicon data.


<p align="center"><img src="https://github.com/PengfanZhang/Rbec/blob/master/Rbec_workflow.png" height="400" /></p>

Rbec firstly align the dereplicatred reads to the reference sequences and calculate the *k*-mer distance between each read and reference. The reference sequence shows the lowest *k*-mer distance is chosen as a candidate error-producing reference sequence for a given read. To calculate the error matrix, reads are subsampled and are aligned to their corresponding error-producing reference sequences. The error matrix is then calculated based on the occurrence frequency of each transition case under different sequencing qualities. Finally, a *P-value* is calculated for each dereplicated read under the *Poisson* distribution. A *P-value* threshold is used to correct the errorneous reads and the expectation value is used to identify the polymorphic paralogues originating from the same strain.

Content
---

[1. Installation](#Installation)

[2. Usages](#Usages)

   [&nbsp; &nbsp; --Parameters](#Parameters)
    
   [&nbsp; &nbsp; --Output files](#Output-files)
    
   [&nbsp; &nbsp; --Characterizing microbial communities in SynComs](#Characterizing-microbial-communities-in-SynComs)
    
   [&nbsp; &nbsp; --Detecting contaminants](#Detecting-contaminants)
    
   [&nbsp; &nbsp; --Dependence on accurate reference sequences](#Dependence-on-accurate-reference-sequences)

[3. Demo starting from raw data](#Demo-starting-from-raw-data)

[4. Credits](#Credits)


## Installation

Rbec is a free R package. To install Rbec, one needs to have the 'dada2' package (https://benjjneb.github.io/dada2/dada-installation.html) installed beforehand manually with the following command:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```

Then you can directly download the newest version with devtools:
```
devtools::install_github("PengfanZhang/Rbec", dependencies = TRUE)
```

Alternatively, you can download from Bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rbec")
```

## Usages


### Parameters


`fastq`: path of the FASTQ file containing the merged amplicon sequencing reads (Ns are not allowed in the reads).

`reference`: path to the nucleotide FASTA file containing the unique reference sequences. Each sequence must be in one line (Ns are not allowed).

`outdir`: path to the output directory.

`threads`: number of threads used, default 1.

`sampling_size`: sampling size for calculating the error matrix, default 5,000.

`ascii`: ASCII characters used to encode phred scores (33 or 64), default 33.

`min_cont_obs_abd`: minimum observed abundance for a non-corrected unique tag to be considered a potential contaminant (default 200).

`min_cont_abd`: minimum relative abundance for a non-corrected unique tag to be considered a potential contaminant (default 0.03).

`min_E`: minimum value expectation of the Poisson distribution for identification of paralogues (default 0.05).

`min_P`: minimum P-value of the Poisson distribution required to correct a sequencing read (default 1e-40).

`ref_seeker`: method used for finding the candidate error-producing reference sequence for a tag showing identical lowest k-mer distance to multiple references. 1 for the abundance-based method; 2 for the transition probability-based method, default 1.

`cn`: path to a tab-separated table containing the copy number of the marker gene for each strain. The first column of the table should contain the strain IDs, and the second the copy number estimates (real values are allowed). If no table is provided, Rbec will normalize the abundance based on the internally inferred copy number (which tends to be un underestimate). Rbec will output the normalized as well as the non-normalized abundance tables. Default NULL.


### Output files

`strain_table.txt`: tab-separated file containing the strain-level community composition table.

`strain_table_normalized.txt`: tab-separated file strain-level community composition table, normalized using copy-number information (see also documentation on the parameter ‘cn’).

`contamination_seq.fna`: nucleotide FASTA file containing the sequences of potential contaminants.

`rbec.log`: text file containing the percentage of corrected reads per sample, which can be used to predict potentially contaminated (see also parameter ‘min_cont_abd’).

`paralogue_seq.fna`: nucleotide FASTA file including sequences of polymorphic paralogues sequences found for each strain.

`lambda_final.out`: test file containing the lambda and P-values of the Poisson distribution for each unique tag.

`error_matrix_final.out`: text file containing the error matrix used in the final iteration.


### Characterizing microbial communities in SynComs

You can use the following commands to profile the microbial composition and the test data is a subset of amplicon sequencing data from individual strain:

```
fq <- system.file("extdata", "test_raw_merged_reads.fastq.gz", package="Rbec")

ref <- system.file("extdata", "test_ref.fasta", package="Rbec")

Rbec(fq, ref, tempdir(), 1, 500, 33)
```

Users can merge the profiling tables from different samples into a single one by using 'mergeSequenceTables' function in DADA2.

### Detecting contaminants

One of the main sources of technical variation in gnotobiotic experiments is caused by microbial contaminations occurring during the development of the experiment or already present during input SynCom preparation. One of the features of Rbec is the assessment of likely contaminated samples based the recruitment ratio of sequencing reads across samples. When analysing data with Rbec, a separate log file is provided as an output for each sample. To predict contaminated samples, a text file containing a list of all paths to the log files needs to be provided before running the following command:

```
log_path <- list.files(paste(path.package("Rbec"), "inst/extdata/contamination_test", sep="/"), recursive=TRUE, full.names=TRUE)

log_file <- tempfile()

writeLines(log_path, log_file)

Contam_detect(log_file, tempdir())
```

This command will generate a plot showing the distribution of percentages of corrected reads across the whole sample set and a log file with predicted contaminated samples are generated. As a general rule, 90% or more of reads should be corrected in clean SynCom samples.

### Dependence on accurate reference sequences


Since Rbec is a reference-based method for error correction in sequencing reads, the accuracy of the reference sequence would critically influence the result. Inference of marker gene sequences from draft genome assemblies or via Sanger sequencing might lead to errors that may negatively impact the accuracy of the results. If errors are present in the reference sequence of a certain strain and no reads can be perfectly aligned to that reference (with an initial abundance of 0 for that reference), Rbec flags this strain as absent with an abundance of 0. One tricky way to overcome this issue in which references are not 100% accurate is that you can look up the contamination sequences outputted by the Rbec function at the first round and align the contamination sequences to the references of strains that are missing in the community. When this occurs, the correct sequence will be flagged by Rbec as a putative contaminant. The user can use this information to replace the erroneous entry in the reference database and re-run Rbec. This should be done with care to avoid true contaminants to be mistaken by members of the original input SynCom when their phylogenetic distance is low.

## Demo starting from raw data

In the `testdata` folder, you can find the rawdata and the corresponding reference file from a real synthetic community. In this section, we'll guide you through the preprocess of the rawdata and the reference sequences.

Rbec can use as an input FASTQ files from either single- or pair-end sequencing data. In addition, a database of reference sequences in FASTA format needs to be provided.

### Preparation of reference sequences
Sequences in the reference database should be already truncated to the region exactly matching the amplicon reads. For example, If we sequence the V5-V7 region of the synthetic bacterial community, the reference sequence of each strain in the reference database should also be truncated to retain the V5-V7 region only rather than the full-length of *16S* rRNA sequences. To truncate the reference sequences into a specific region, a tool such as cutadapt can be used:

```
cutadapt -g AACMGGATTAGATACCC -a GGAAGGTGGGGATGACGT -n 2 -o testdata/reference_V5V7.fasta testdata/reference_full.fasta
```

Rbec currently only supports the reference database in a non-wrapped format. Towards this end, the following Unix command can be used to transform the wrapped sequences into non-wrapped ones.

```
awk 'BEGIN{seqs=""}{if(/^>/){if(seqs!=""){print seqs;seqs=""}; print $0} else{seqs=seqs""$0}}END{print seqs}' testdata/reference_V5V7.fasta > testdata/reference_V5V7_nonwrapped.fasta
```
or use `Seqkit`
```
seqkit seq -w 0 -o testdata/reference_V5V7_nonwrapped.fasta testdata/reference_V5V7.fasta
```

### Preprocess of the rawdata

Firstly, pair-end raw reads should be merged by using any reads merging tools:
```
flash2 testdata/test_R1.fastq.gz testdata/test_R2.fastq.gz -M 250 -o test -x 0.25 -d testdata
```
usearch -fastq_filter test.extendedFrags.fastq -fastqout test.extendedFrags.filtered.fastq -fastq_maxns 0 -fastq_minlen 200



Amplicon reads should be manually filtered by excluding reads with ambiguous reads with USEARCH or a similar software before running Rbec and short reads are removed to avoid dimers. Following is an example of removing reads with ambiguous reads with USEARCH:
```
usearch -fastq_filter testdata/test.extendedFrags.fastq -fastqout testdata/test.extendedFrags.filtered.fastq -fastq_maxns 0 -fastq_minlen 200
```

### Start with Rbec

Now, with the two files in hands, you can invoke Rbec in your R environment to explore your samples:
```
Rbec("testdata/test.extendedFrags.filtered.fastq", "testdata/reference_V5V7_nonwrapped.fasta", 1, 2000)
```

## Credits

Pengfan Zhang (pzhang@mpipz.mpg.de); Ruben Garrido-Oter (garridoo@mpipz.mpg.de)
