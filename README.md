# Rbec: a tool for analysis of amplicon sequencing data from synthetic microbial communities

Our manuscript is currently available as a preprint in [*bioRxiv*](https://doi.org/10.1101/2021.01.15.426834).

Rbec is a tool for analysing amplicon sequencing data from synthetic communities (SynComs), where the reference sequences for each strain are already available. Rbec can accurately correct PCR and sequencing errors, identify intra-species polymorphic variation, and detect contaminations in SynCom amplicon data.

Installation
---

Rbec is a free R package. To install Rbec, one needs to have the 'dada2' package (https://benjjneb.github.io/dada2/dada-installation.html) installed beforehand manually with the following command:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2")
```

Then you can directly download the source package (Rbec_0.1.0.tar.gz) and install it locally:
```
install.packages("Rbec_0.1.0.tar.gz", type='source')
```

The package is under review in Bioconductor and soon it will be available for downloading there.

How to use Rbec?
---

1. Characterizing microbial communities in SynComs
---

Rbec can use as an input FASTQ files from either single- or pair-end sequencing data. In addition, a database of reference sequences in FASTA format needs to be provided.

Amplicon reads should be manually filtered by excluding reads with ambiguous reads with USEARCH or a similar software before running Rbec. Following is an example of removing reads with ambiguous reads with USEARCH:
```
usearch -fastq_filter <input_fastq> -fastqout <output_fastq> -fastq_maxns 0
```

Sequences in the reference database should be already truncated to the region matching the amplicon reads, which should be flanked by the primer sequences. For example, If we sequence the V5-V7 region of the synthetic bacterial community, the reference sequence of each strain in the reference database should also be truncated to retain the V5-V7 region only rather than the full-length of *16S* rRNA sequences. To truncate the reference sequences into a specific region, a tool such as cutadapt can be used:

```
cutadapt -g <forward_sequencing_primer> -a <reverse_sequencing_primer> -n 2 -o <output> <input>
```

Rbec currently only supports the reference database in a non-wrapped format. Towards this end, the following Unix command can be used to transform the wrapped sequences into non-wrapped ones.

```
awk 'BEGIN{seqs=""}{if(/^>/){if(seqs!=""){print seqs;seqs=""}; print $0} else{seqs=seqs""$0}}END{print seqs}' <input_reference_database> > <output_reference_database>
```

After preparation of the input files, Rbec can be run to estimate SynCom community profiles. An example with a small test dataset from a single bacterial strain illustrates this process:

```
fq <- system.file("extdata", "test_raw_merged_reads.fastq", package="Rbec")

ref <- system.file("extdata", "test_ref.fasta", package="Rbec")

Rbec(fq, ref, "./", 1, 500, 33)
```

2. Detecting contaminants
---

One of the main sources of technical variation in gnotobiotic experiments is caused by microbial contaminations occurring during the development of the experiment or already present during input SynCom preparation. One of the features of Rbec is the assessment of likely contaminated samples based the recruitment ratio of sequencing reads across samples. When analysing data with Rbec, a separate log file is provided as an output for each sample. To predict contaminated samples, a text file containing a list of all paths to the log files needs to be provided before running the following command:

```
Contam_detect(log_file, outdir)
```

This command will generate a plot showing the distribution of percentages of corrected reads across the whole sample set and a log file with predicted contaminated samples are generated. As a general rule, 90% or more of reads should be corrected in clean SynCom samples.

3. Rbec dependence on accurate reference sequences
---

Since Rbec is a reference-based method for error correction in sequencing reads, the accuracy of the reference sequence would critically influence the result. Inference of marker gene sequences from draft genome assemblies or via Sanger sequencing might lead to errors that may negatively impact the accuracy of the results. If errors are present in the reference sequence of a certain strain and no reads can be perfectly aligned to that reference (with an initial abundance of 0 for that reference), Rbec flags this strain as absent with an abundance of 0. One tricky way to overcome this issue in which references are not 100% accurate is that you can look up the contamination sequences outputted by the Rbec function at the first round and align the contamination sequences to the references of strains that are missing in the community. When this occurs, the correct sequence will be flagged by Rbec as a putative contaminant. The user can use this information to replace the erroneous entry in the reference database and re-run Rbec. This should be done with care to avoid true contaminants to be mistaken by members of the original input SynCom when their phylogenetic distance is low.

Credits
---

Pengfan Zhang (pzhang@mpipz.mpg.de); Ruben Garrido-Oter (garridoo@mpipz.mpg.de)
