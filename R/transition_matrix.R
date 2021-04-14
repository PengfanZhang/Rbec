
#' Reference-based error correction of amplicon sequencing data
#'
#' @description
#' This function count the transition matrix
#'
#' @details Ruben Garrido-Oter's group, Plant-Microbe interaction, Max Planck Institute for Plant Breeding Research
#' @author Pengfan Zhang
#'
#' @param query list containing subsampled amplicon sequencing reads, quality scores, and reference sequences showing the highest identity for each read (Ns are not allowed in the reads)
#' @param ascii ascii characters used to encode phred scores
#'
#' @usage trans_m(query, ascii)
#'
#' @return The output is a 20 by 43 matrix containing the counts for different kinds of transitions
#'
#' @noRd
#'

trans_m <- function(query, ascii) {

    # initialize the transition matrix
    trans_matrix <- matrix(0, ncol=43, nrow=20)
    trans_kinds <- c(paste(rep(c("A","T","G","C","-"), each=4),
                           rep(2, 20), rep(c("A","T","G","C"), 5), sep=""))
    rownames(trans_matrix) <- trans_kinds
    colnames(trans_matrix) <- as.character(seq(0, 42))

    # TODO: define outside
    tm_cal <- function(que, qual, ref){

    # perform global alignment between query and reference using NW
    align_out <- nwalign(que, ref)
    Seq1align <- align_out[1]
    Seq2align <- align_out[2]
    Seq1align <- unlist(strsplit(Seq1align, ""))
    Seq2align <- unlist(strsplit(Seq2align, ""))
    align_pairs <- paste(Seq2align, Seq1align, sep="2")
    align_pairs <- align_pairs[!endsWith(align_pairs, "-")]

    # check if the quality system is phred 33 or 64
    qual <- as.character(strtoi(charToRaw(paste((qual), collapse="")), 16L)-ascii)

    align_data <- data.frame(alig=align_pairs, qual=qual, stringsAsFactors=FALSE)
    return(align_data)

    }

    # for all query sequences, calculate global alignments to reference sequence and corresponding transitions
    align_data <- do.call(rbind, apply(query, 1, function(x) tm_cal(x[1], x[2], x[3])))

    # transform transitions for each aligned querity to transition matrix
    align_data <- data.frame(table(align_data), stringsAsFactors=FALSE)
    for(i in 1:nrow(align_data)){
      trans_matrix[as.character(align_data[i, 1]), as.character(align_data[i, 2])] <- align_data[i, 3]
    }

    return(trans_matrix)

}

