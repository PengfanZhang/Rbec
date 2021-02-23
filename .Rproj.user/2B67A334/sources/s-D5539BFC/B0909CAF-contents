#' Reference-based error correction of amplicon sequencing data
#'
#' @description
#' This function calculates the abundance probabilities for each reads using poisson distribution
#'
#' @details Ruben Garrido-Oter's group, Plant-Microbe interaction, Max Planck Institute for Plant Breeding Research
#' @author  Pengfan Zhang
#'
#' @param derep dereplicated reads (Ns are not allowed in the reads)
#' @param ref the unique reference sequences of the reference seqeunces, each sequence must be in one line (Ns are not allowed in the sequences)
#' @param error_matrix The error matrix from the former iteration
#'
#' @usage abd_prob(derep, ref, error_matrix)
#' @examples
#'
#' @return Returns the lambda value and pvalue for each reads
#'
#' @export
#'

abd_prob <- function(derep, ref, error_matrix) {

  lambda_out <- data.frame(numeric(0))
  seqs <- names(derep[["uniques"]])
  abd <- derep[["uniques"]]
  seqs <- toupper(seqs)
  bestref <- derep[["bestref"]]

  input_data <- data.frame(seqs, abd, bestref, stringsAsFactors=F)

    # calculate the estimated abundance of reference sequences by summing up the abundances of
    # reads potentailly produced from that reference sequence
    ref_abd <- data.frame(bestref, abd)
    ref_abd <- data.frame(tapply(ref_abd[, -1], INDEX=list(ref_abd[,1]), sum))
    ref$est_abd[match(rownames(ref_abd), ref[, 1])] <- ref_abd[, 1]

    # function for lambda calculation
    lambda_calculation <- function(seq, abd, bestref) {
      best_match <- match(bestref, ref[, 1])
      align_out <- nwalign(seq, bestref)
      Seq1align <- align_out[1]
      Seq2align <- align_out[2]
      Seq1align <- unlist(strsplit(Seq1align, ""))
      Seq2align <- unlist(strsplit(Seq2align, ""))
      align_pairs <- paste(Seq2align, Seq1align, sep="2")
      align_pairs <- align_pairs[! endsWith(align_pairs, "-")]
      qual <- round(derep[["quals"]][match(seq, names(derep[["uniques"]])), ], digits=0)
      qual <- as.character(qual[!is.na(qual)])
      retain_index <- which(!align_pairs %in% c("A2A", "C2C", "G2G", "T2T"))
      align_pairs <- align_pairs[retain_index]
      qual <- qual[retain_index]
      data4lambda <- data.frame(align_pairs, qual)
      tp <- apply(data4lambda, 1, function(x) error_matrix[x[1], x[2]])
      lambda0 <- prod(tp)
      lambda_out <- rbind(c(lambda0, best_match, abd))
      return(lambda_out)
    }

    # calculate lambdas for each dereplicated sequence
    lambda_out <- apply(input_data, 1, function(x) lambda_calculation(x[1], x[2], x[3]))
    lambda_out <- t(lambda_out)
    lambda_out <- apply(lambda_out, 2, function(x) as.numeric(x))
    lambda_out <- as.data.frame(lambda_out, stringsAsFactors=F)

    # calculate abundance probabilities
    names(lambda_out) <- c("Lambda", "Ref_ID", "Obs_abd")
    lambda_out$ref_abd <- ref[lambda_out[, 2], "est_abd"]
    lambda_out$E <- lambda_out$Lambda*lambda_out$ref_abd
    lambda_out$pval <- apply(lambda_out, 1, function(x) 1-sum(dpois(0:(x[3]-1), x[5])))
    return(lambda_out)
}

