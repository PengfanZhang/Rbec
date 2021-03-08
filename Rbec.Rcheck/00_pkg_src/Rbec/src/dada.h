#ifndef _Rbec_H_
#define _Rbec_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
//#include <gsl/gsl_cdf.h>

#define SEQLEN 9999

void int2nt(char *oseq, const char *iseq);
char *intstr(const char *iseq);
void assign_kmer(uint16_t *kvec, const char *seq, int k);
double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k);

#endif
