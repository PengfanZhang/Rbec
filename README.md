# Rbec

# Reference-based error correction of amplicon sequencing data from synthetic communities

Rbec will be available for downloading soon......

Rbec is an adapted version of DADA2 for analyzing amplicon sequencing data from synthetic communities (SynComs), where the reference sequences for each strain exists. Rbec can not only accurately profile the microbial compositions in SynComs, but also predict the contaminants in SynCom samples.

Installation
---

Rbec is a free R package. To install Rbec, one needs to have the 'dada2' package (https://benjjneb.github.io/dada2/dada-installation.html) installed beforehand manually with the following command:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")
```

Then you can directly download the source package (Rbec_0.1.0.tar.gz) and install it locally:
```
install.packages("Rbec_0.1.0.tar.gz", type='source')
```


How to use Rbec?
---

1.Characterizing microbial communities in SynComs
---

Fastq file from single-end sequencing or merged fastq file from pair-end sequencing and reference database comprising the input data can be fed into Rbec.

Amplicon reads should be manually filtered by excluding reads with ambiguous reads with Usearch or other software before being fed into Rbec. Following is an example of removing reads with ambiguous reads with Usearch:
```
usearch -fastq_filter <input_fastq> -fastqout <output_fastq> -fastq_maxns 0
```

Sequences in the reference database should be already truncated to the region perfectly matching the amplicon reads. For example, If we sequence the V5-V7 region of the synthetic bacterial community, the reference sequence of each strain in the reference database should also be truncated to retain the V5-V7 region only rather than the full-length of 16S sequences. To truncate the reference sequences into a specific region, cutadapt can be used:
```
cutadapt -g <forward_sequencing_primer> -a <reverse_sequencing_primer> -n 2 -o <output> <input>
```

What's more, Rbec currently only supports the reference database in a non-wrapped format, which means a sequence is not splited up to take up multiple lines. Towards this end, you can use the following Unix command to transform the wrapped sequences into non-wrapped ones.
```
awk 'BEGIN{seqs=""}{if(/^>/){if(seqs!=""){print seqs;seqs=""}; print $0} else{seqs=seqs""$0}}END{print seqs}' <input_reference_database> > <output_reference_database>
```

After both input files are well prepared, you can start to run Rbec to decipher the microbial compositions in SynComs. The test data provided in the package only includes amplicon data from one strain.
```
fq <- system.file("extdata", "test_raw_merged_reads.fastq", package = "Rbec")

ref <- system.file("extdata", "test_ref.fasta", package = "Rbec")

Rbec(fq, ref, "./", 1, 500, 33)
```

2.Detecting comtaminants
---

For SynCom samples, we have to be cautious about contaminations that might confound the experimental outputs. If you run a bunch of SynCom samples with Rbec, you can evaluate whether the samples are contaminated with unknown taxa by sorting the recruitment ratio of reads by Rbec for each sample. The idea is that if some samples with extremely low recruitment ratios of reads and deviate from other samples, we can rationally assume those samples were contaminated.

Any analyzing the data with the Rbec function, a log file was outputted for each sample, which documents the percentage of corrected reads in the sample. To predict the contaminated samples, you need to intergrate the paths to all the log files into one file and run the following command:
```
Contam_detect(log_file, outdir)
```

After running the above command, a plot showing the distribution of percentages of corrected reads across the whole sample set and a log file recording contaminated samples are generated. Emperically, almost 90% or more of reads can be corrected in clean SynCom samples.

3.Dependence on accurate reference sequences
---

Since Rbec is a reference-based method for error correction in sequencing reads, the accuracy of the reference sequence would profoundly influence the result. We could sequence the marker gene of strains in the SynCom beforehand with sanger sequencing or predict them from their whole genome sequencing, however, neither of the methods is error-free. If errors are present in the reference sequence of a certain strain and no reads can be peferctly aligned to that reference (with an initial abundance of 0 for that reference), Rbec flags this strain as absence with an abundance of 0. 

One tricky way to overcome this issue in which references are not 100% accurate is that you can look up the contamination sequences outputted by the Rbec function at the first round and align the contamination sequences to the references of strains that are missing in the community. If you find one contamination sequence is extremely close to one reference sequence you have, that contamination sequence should be the accurate reference for the corresponding strain. Then you can replace the errorneous reference with the new one and re-run Rbec on all the samples. We also acknowledge that the missing strains are attribute to the low abundance and low possibility to be detected by sequencing and the close contamination sequence might come from a true contaminant, but we advocate that the possibility to randomly include contaminants that are close to any one of the strains in the SynCom is low.

Credits
---

Pengfan Zhang (pzhang@mpipz.mpg.de); Ruben Garrido-Oter (garridoo@mpipz.mpg.de)

