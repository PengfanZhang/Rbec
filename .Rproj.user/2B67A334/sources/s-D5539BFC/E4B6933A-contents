#' Reference-based error correction of amplicon sequencing data
#'
#' @description
#' This function fits the loess regression to the error matrix
#'
#' @details Ruben Garrido-Oter's group, Plant-Microbe interaction, Max Planck Institute for Plant Breeding Research
#' @author  Pengfan Zhang
#'
#' @param trans the transition matrix
#' @param min_err_rate the minimum transition probability for each substitution or insertion case
#'
#' @usage loessErr(trans, min_err_rate=1e-07)
#' @examples
#'
#' @return Returns the loess fitted error matrix
#'


loessErr <- function (trans, min_err_rate=1e-07) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A", "C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        mod.lo <- loess(rlogp ~ q, df, weights = tot)
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }

  # add the insertion probabilities (-2A/(-2A+A2A+T2A+C2A+G2A))
  for (ntj in c("A", "C", "G", "T")) {
    errs <- trans[paste0("-", "2", ntj), ]
    tot <- colSums(trans[paste0(c("A", "C", "G", "T", "-"), "2", ntj), ])
    rlogp <- log10((errs + 1)/tot)
    rlogp[is.infinite(rlogp)] <- NA
    df <- data.frame(q = qq, errs = errs, tot = tot,
                     rlogp = rlogp)
    mod.lo <- loess(rlogp ~ q, df, weights = tot)
    pred <- predict(mod.lo, qq)
    maxrli <- max(which(!is.na(pred)))
    minrli <- min(which(!is.na(pred)))
    pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
    pred[seq_along(pred) < minrli] <- pred[[minrli]]
    est <- rbind(est, 10^pred)
  }

  est[est < min_err_rate] <- min_err_rate

  # calculate the transition probabilities for ('A2A', 'T2T', 'C2C', 'G2G') and fill out the error matrix
  err <- rbind(1-colSums(est[1:3, ]), est[1:3, ], est[4, ], 1-colSums(est[4:6, ]), est[5:6, ], est[7:8, ],
               1-colSums(est[7:9, ]), est[9, ], est[10:12, ], 1-colSums(est[10:12, ]), est[13:16, ])

  rownames(err) <- paste0(rep(c("A", "C", "G", "T", "-"), each = 4),
                          "2", c("A", "C", "G", "T"))

  colnames(err) <- colnames(trans)
  return(err)

}

