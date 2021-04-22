#' Reference-based error correction of amplicon sequencing data
#'
#' @description
#' This function calculate the error matrix
#'
#' @details Ruben Garrido-Oter's group, Plant-Microbe interaction, Max Planck Institute for Plant Breeding Research
#' @author  Pengfan Zhang
#'
#' @param fq the path of merged amplicon sequencing reads in fastq format (Ns are not allowed in the reads)
#' @param ref the unique reference sequences of the reference seqeunces, each sequence must be in one line (Ns are not allowed in the sequences)
#' @param sample_size the sampling size of reads to generate the transition matrix
#' @param threads the number of threads used to align the query reads to reference sequences
#' @param ascii ascii characters used to encode phred scores
#'
#'
#' @import dada2
#' @import foreach
#' @import doParallel
#'
#' @usage error_m(fq, ref, sample_size, threads, ascii)
#'
#' @return The output is a 20 by 43 transition probability matrix
#'
#' @noRd
#'

error_m <- function(fq, ref, sample_size=10000, threads, ascii) {

    x <- NA
    t1 <- Sys.time()

    # read FASTQ file
    if (grepl(".gz", fq, fixed = TRUE)) {
        raw_data <- read_lines(gzfile(fq))
    }
    else {
        raw_data <- read_lines(fq)
    }
    if (sample_size > length(raw_data)/4){
      stop("The sampling size ", sample_size ," exceeds the total number of reads ", length(raw_data)/4 ," in the input file")
    }

    # calculate the initial abundance of reference sequences
    derep1 <- derepFastq(fq)
    ref$abd <- derep1[["uniques"]][match(ref$ref_seq, names(derep1[["uniques"]]))]
    ref <- ref[!is.na(ref$abd), ]
    ref$est_abd <- 0

    # find the reference sequences with highest similarities for each unique sequence
    uniseqs <- names(derep1[["uniques"]])

    findbest <- function(que) {
        # use the kmer distance function implemented in DADA2
        kdist <- function(que, ref, kmer=7) {
            dis <- kmer_dist(que, ref, kmer)
            return(dis)
        }
        # calculate k-mer distances for all sequences
        kdist_res <- unlist(lapply(ref$ref_seq, function(x) kdist(que, x)))
        best_match <- which.min(kdist_res)
        # if there are more than 1 hits showing the same maximum identity,
        # the reference with the highest abundance is chosen as the best reference sequence
        if (length(best_match) > 1) {
          ref_candidate <- ref[best_match, ]
          ref_best <- ref_candidate[which.max(ref_candidate[, 2]), 1]
        }
        else {
          ref_best <- ref[best_match, 1]
        }
        return(c(que, ref_best))
    }

    #cl <- makeCluster(threads)
    registerDoParallel(cores=threads)
    que2ref <- foreach::foreach(x=uniseqs, .combine='rbind') %dopar% findbest(x)
    stopImplicitCluster()
    que2ref <- as.character(que2ref[, 2])
    message("Finished finding the best reference sequences for each unique sequences.")
    t2 <- Sys.time()
    message(t2-t1)

    # sample raw sequences for transition matrix generation
    sample_label <- sample(seq_len(length(raw_data)/4), sample_size)
    seqs <- raw_data[sample_label*4-2]
    seqs <- toupper(seqs)
    qual <- raw_data[sample_label*4]
    bestref <- que2ref[derep1[["map"]][sample_label]]
    query <- data.frame(seqs=seqs, qual=qual, refs=bestref)

    # calculate the transition matrix and error matrix
    trans_matrix <- trans_m(query, ascii)
    error_matrix <- loessErr(trans_matrix)

    derep1[["bestref"]] <- que2ref

    error_ref_matrix <- list(ref=ref, err=error_matrix, derep=derep1, total_reads=length(raw_data)/4)
    return(error_ref_matrix)

}

