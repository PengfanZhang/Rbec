#include "dada.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


//' @name kmer_dist
//' @title DADA2
//' @description Calculate the kmer distance between two sequences. This function is originally developed by DADA2 in 2012 and we adapt it here.
//'
//' @param s1 A \code{character(1)} of DNA sequence 1.
//' @param s2 A \code{character(1)} of DNA sequence 2.
//' @param kmer_size Kmer size.
//' @return The kmer distance between two sequences
//' @noRd
//'
// [[Rcpp::export]]

NumericVector kmer_dist(std::vector< std::string > s1, std::vector< std::string > s2, int kmer_size) {
  size_t len1 = 0, len2 = 0;
  char *seq1, *seq2;
  size_t n_kmers = (1 << (2*kmer_size));  // 4^k kmers

  size_t nseqs = s1.size();
  if(nseqs != s2.size()) { Rcpp::stop("Mismatched numbers of sequences."); }

  Rcpp::NumericVector kdist(nseqs);
  uint16_t *kv1 = (uint16_t *) malloc(n_kmers * sizeof(uint16_t)); //E
  uint16_t *kv2 = (uint16_t *) malloc(n_kmers * sizeof(uint16_t)); //E
  if(kv1 == NULL || kv2 == NULL) Rcpp::stop("Memory allocation failed.");

  for(int i=0;i<nseqs;i++) {
    seq1 = intstr(s1[i].c_str());
    len1 = s1[i].size();
    assign_kmer(kv1, seq1, kmer_size);
    seq2 = intstr(s2[i].c_str());
    len2 = s2[i].size();
    assign_kmer(kv2, seq2, kmer_size);
    kdist[i] = kmer_dist(kv1, len1, kv2, len2, kmer_size);
    free(seq2);
    free(seq1);
  }

  free(kv1);
  free(kv2);
  return(kdist);
}
