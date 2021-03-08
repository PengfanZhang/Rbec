#include <stdlib.h>
#include "dada.h"
#include "emmintrin.h"
using namespace std;
// [[Rcpp::interfaces(cpp)]]



void nt2int(char *oseq, const char *iseq) {
  int i, len = strlen(iseq);

  for (i = 0; i < len; i++, iseq++, oseq++) {
    switch (*iseq) {
    case 'A':
      *oseq = 1;
      break;
    case 'C':
      *oseq = 2;
      break;
    case 'G':
      *oseq = 3;
      break;
    case 'T':
      *oseq = 4;
      break;
    case 'N':
      *oseq = 5;
      break;
    case '-':
      *oseq = '-';
      break;
    default:
      Rprintf("invalid character in input:%c.\n",*iseq);
      *oseq = '\0';
    }
  }
  *oseq = '\0';
  return;
}

char *intstr(const char *iseq) {
  char *oseq = (char *) malloc(strlen(iseq)+1); //E
  if (oseq == NULL)  Rcpp::stop("Memory allocation failed!\n");

  strcpy(oseq, iseq);
  nt2int(oseq, oseq);
  return oseq;
}


double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  int i;
  int n_kmer = 1 << (2*k); // 4^k kmers
  uint16_t dotsum = 0;
  double dot = 0.0;

  for(i=0;i<n_kmer; i++) {
    dotsum += (kv1[i] < kv2[i] ? kv1[i] : kv2[i]);
  }

  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);

  return (1. - dot);
}


void assign_kmer(uint16_t *kvec, const char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <= 0 || len > SEQLEN) { Rcpp::stop("Unexpected sequence length."); }
  if(k >= len || k < 3 || k > 8) { Rcpp::stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  for(kmer=0;kmer<n_kmers;kmer++) { kvec[kmer] = 0; }

  if(len <=0 || len > SEQLEN) {
    Rcpp::stop("Unexpected sequence length.");
  }

  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        Rcpp::stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }

    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; }
    else if(kmer >= n_kmers) {
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
}

