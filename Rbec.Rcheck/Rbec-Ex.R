pkgname <- "Rbec"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Rbec')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rbec")
### * Rbec

flush(stderr()); flush(stdout())

### Name: Rbec
### Title: Reference-based error correction of amplicon sequencing data
### Aliases: Rbec

### ** Examples

fastq <- system.file("extdata", "test_raw_merged_reads.fastq.gz", package = "Rbec")

ref <- system.file("extdata", "test_ref.fasta", package = "Rbec")

Rbec(fastq=fastq, reference=ref, outdir="./", threads=1, sampling_size=500, ascii=33)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
