/****************************************************************************
 * Copyright (C) 2017 by Paula Perez Rubio                                  *
 * Copyright (C) 2022 by Gerben Voshol
 *                                                                          *
 * This file is part of FastqPuri.                                      *
 *                                                                          *
 *   FastqPuri is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as                *
 *   published by the Free Software Foundation, either version 3 of the     *
 *   License, or (at your option) any later version.                        *
 *                                                                          *
 *   FastqPuri is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with FastqPuri.                                              *
 *   If not, see <http://www.gnu.org/licenses/>.                            *
 ****************************************************************************/

/**
 * @file splitDS.c
 * @author Gerben Voshol <gpvoshol@gmail.com>
 * @date 10.01.2022
 * @brief splitDS main function
 *
 * This file contains the splitDS main function.
 * See README_splitDS.md for more details.
 *
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <omp.h>
#include "defines.h"
#include "config.h"
#include "str_manip.h"
#include "fopen_gen.h"
#include "trimDS.h"
#include "trim.h"
#include "tree.h"
#include "bloom.h"
#include "Lmer.h"
#include "adapters.h"
#include "fq_read.h"
#include "io_trimFilterDS.h"
#include "init_trimFilterDS.h"
#include "timing.h"

uint64_t alloc_mem = 0;    /**< global variable. Memory allocated in the heap.*/
Iparam_trimFilter par_TF;  /**< global variable: Input parameters of makeTree.*/

/**
 * @brief Function that prints trimFilterPE help dialog when called.
*/
void printHelpDialog_splitDS() {
  const char dialog[] =
   "Usage: splitPE --ifq <INPUT1.fq>:<INPUT2.fq>\n"
   "                  --output [O_PREFIX] --gzip [y|n]\n"
   "                  (--idx [<INDEX_FILE(s)>:<score>:<lmer_len>] |\n"
   "Reads in paired end fq files (gz, bz2, z formats also accepted) "
   "Options:\n"
   " -v, --version prints package version.\n"
   " -h, --help    prints help dialog.\n"
   " -f, --ifq     2 fastq input files [*fq|*fq.gz|*fq.bz2] separated by\n"
   "               colons, mandatory option.\n"
   " -l, --length  read length: length of the reads, mandatory option.\n"
   " -o, --output  output prefix (with path), optional (default ./out).\n"
   " -b, --best    output to filter with the best score\n"
   " -x, --idx     index input file. To be included with any methods to remove.\n"
   "               contaminations (TREE, BLOOM). Minimum 3 fields separated by colons: \n"
   "               <INDEX_FILE(s)>: output of makeTree, makeBloom,\n"
   "               <score>: score threshold to accept a match [0,1],\n"
   "               [lmer_len]: the length of the lmers to be \n"
   "                        looked for in the reads [1,READ_LENGTH].\n"
   " -A, --adapter adapter input. Four fields separated by colons:\n"
   "               <AD1.fa>: fasta file containing adapters,\n"
   "               <AD2.fa>: fasta file containing adapters,\n"
   "               <mismatches>: maximum mismatch count allowed,\n"
   "               <score>: score threshold  for the aligner.\n"
   " -Q, --trimQ   NO:       does nothing to low quality reads (default),\n"
   "               ALL:      removes all reads containing at least one low\n"
   "                         quality nucleotide.\n"
   "               ENDS:     trims the ends of the read if their quality is\n"
   "                         below the threshold -q,\n"
   "               FRAC:     discards a read if the fraction of bases with\n"
   "                         low quality scores (below -q) is over 5 percent \n"
   "                         or a user defined percentage (-p). \n"
   "               ENDSFRAC: trims the ends and then discards the read if \n"
   "                         there are more low quality nucleotides than \n"
   "                         allowed by the option -p.\n"
   "               GLOBAL:   removes n1 cycles on the left and n2 on the \n"
   "                         right, specified in -g.\n"
   "               All reads are discarded if they are shorter than MINL\n"
   "               (specified with -m or --minL).\n "   
   " -g, --trimG   <min_len>: enable trimming of polyG repeats with a minimum length of min_len\n"
   "                       at the end of read,\n"
   "               All reads are discarded if they are shorter than MINL\n"
   "               (specified with -m or --minL).\n "   
   " -m, --minL    minimum length allowed for a read before it is discarded\n"
   "               (default 36).\n";
  fprintf(stderr, "%s", dialog);
}

/**
 * Checks if string contains the suffix
 * 
 * */
int strsuff(const char *s, const char *suff) {
  size_t slen = strlen(s);
  size_t sufflen = strlen(suff);

  return slen >= sufflen && !memcmp(s + slen - sufflen, suff, sufflen);
}

// Trim 3prime (tail) of sequence if it contains a poly X of min_size length (default 10)
int trim_polyX(Fq_read *seq, const char X, const int min_size)
{
    const int one_mismatch_per = 8; // Allow one mismatch for each 8 base pairs
    const int max_mismatch = 5;     // Allow maximum of 5 mismatches

    int slen = seq->L - 1;
    int mismatch = 0;
    int i, allowed_mismatch, firstX = slen;
    for(i = slen; i >= 0; i--) {
        // Found another X update position
        if(seq->line2[i] == X) {
            firstX = i;
        // Not an X, add mismatch and check if we continue 
        } else {
          mismatch++;        
        }

        allowed_mismatch = (slen - i + 1) / one_mismatch_per;
        // NOTE: Leaving the minimum length requirement out, might lead to better quality trimming. More testing is required
        if (mismatch > max_mismatch || (mismatch > allowed_mismatch && (slen - i + 1) >= min_size)) {
            break;
        }
    }

    // If the poly X tail is long enough, trim the sequence
    if ((slen - firstX + 1) >= min_size) {
        seq->L = firstX + 1;
        seq->line2[firstX] = '\0';
        seq->line4[firstX] = '\0';
        return 1;
    }

    return 0;
}


/**
 * @brief contains splitDS main function. See README_splitDS.md
 *        for more details.
 *
 * */
int main(int argc, char *argv[]) {

  Timing t;

  timing_start(&t);

  // Read in command line arguments
  fprintf(stderr, "splitPE from FastqPuri\n");

  static struct option long_options[] = {
     {"version", no_argument, 0, 'v'},
     {"help", no_argument, 0, 'h'},
     {"best", no_argument, 0, 'b'},
     {"ifq", required_argument, 0, 'f'},
     {"idx", required_argument, 0, 'x'},
     {"adapter", required_argument, 0, 'A'},
     {"length", required_argument, 0, 'l'},
     {"minL", required_argument, 0, 'm'},
     {"trimG", required_argument, 0, 'g'},
     {"trimQ", required_argument, 0, 'Q'},
     {"percent", required_argument, 0, 'p'},
     {"output", required_argument, 0, 'o'}
  };

  char *Ifq = NULL, *Ifq2 = NULL;
  char *ad_fa = NULL, *ad2_fa = NULL;
  double score = 0.15;
  int kmersize = 25;
  int nrfilters = 0;
  char *prefix = "";
  int L = 151;
  int use_best = 0;
  int adapter_trim = 0;
  int mismatches = 2;
  double threshold = 8.0; 
  int option;
  int i;
  par_TF.minL = 36;
  int min_poly_G = 0;

  par_TF.trimQ = NO;
  par_TF.percent  = 5;
  int method_len = 20;
  Split index = { 0 }, in_fq = { 0 }, adapt = { 0 };
  while ((option = getopt_long(argc, argv, "hvbf:x:o:l:A:m:g:p:Q:", long_options, 0)) != -1) {
    switch (option) {
      case 'h':
        printHelpDialog_splitDS();
        exit(EXIT_SUCCESS);
        break;
      case 'v':
        printf("splitPE version %s \nWritten by Gerben Voshol\n", VERSION);
        exit(EXIT_SUCCESS);
        break;
      case 'b':
        use_best = 1;
        break;
      case 'm':
        par_TF.minL = atoi(optarg);
        break;
      case 'g':
        min_poly_G = atoi(optarg);
        break;
      case 'f':
         in_fq = strsplit(optarg, ':');
         if (in_fq.N != 2) {
            fprintf(stderr, "--ifq, -f: optionERR. You must pass two \n");
            fprintf(stderr, "  arguments separated by semicolons: \n");
            fprintf(stderr, "  <INPUT1.fq>:<INPUT2.fq>\n");
            fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
         }
         Ifq = (char*) malloc(MAX_FILENAME*sizeof(char));
         Ifq2 = (char*) malloc(MAX_FILENAME*sizeof(char));
         strncpy(Ifq, in_fq.s[0], MAX_FILENAME);
         strncpy(Ifq2, in_fq.s[1], MAX_FILENAME);
         break;
      case 'A':
         adapter_trim = true;
         adapt = strsplit(optarg, ':');
         if (adapt.N != 4) {
            fprintf(stderr, "--adapter,-A: optionERR. You must pass four \n");
            fprintf(stderr, "  arguments separated by semicolons: \n");
            fprintf(stderr, "   <AD1.fa>:<AD2.fa>:<mismatches>:<threshold>\n");
            fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
         }
         ad_fa = adapt.s[0];
         ad2_fa = adapt.s[1];
         par_TF.ad.mismatches = atoi(adapt.s[2]);
         par_TF.ad.threshold = atof(adapt.s[3]);
         mismatches = atoi(adapt.s[2]);
         threshold = atof(adapt.s[3]);
         break;
      case 'x':
         index = strsplit(optarg, ':');
         if (index.N < 3) {
            fprintf(stderr, "--idx,-x: optionERR. You must pass at least \n");
            fprintf(stderr, "  three arguments separated by semicolons: \n");
            fprintf(stderr, "  <INDEX_FILE.fa>:<score>:<lmer_len>\n");
            fprintf(stderr, "  and you passed %d\n", index.N);
            fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
         }
         nrfilters = index.N - 2;
         score = atof(index.s[index.N - 2]);
         par_TF.score = score; // Used by tim.c to score sequences
         kmersize = atoi(index.s[index.N - 1]);
         par_TF.kmersize = atoi(index.s[index.N - 1]);
         break;
      case 'o':
         prefix = optarg;
         break;
      case 'q':
         par_TF.minQ = atoi(optarg);
         break;
      case 'Q':
         par_TF.trimQ = (!strncmp(optarg, "NO", method_len)) ? NO :
            (!strncmp(optarg, "ALL", method_len)) ? ALL :
            (!strncmp(optarg, "FRAC", method_len)) ? FRAC :
            (!strncmp(optarg, "ENDS", method_len)) ? ENDS :
            (!strncmp(optarg, "ENDSFRAC", method_len)) ? ENDSFRAC :
            (!strncmp(optarg, "GLOBAL", method_len)) ? GLOBAL : ERROR;
         break;
      case 'p':
         par_TF.percent = atoi(optarg);
         break;
      case 'l':
         L = atoi(optarg);
         par_TF.L = atoi(optarg);
         break;
      default:
        fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                           argv[0], optopt);
        printHelpDialog_splitDS();
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        break;
    }
  }

  if (nrfilters / 2 > NFILES_DS) {
    fprintf(stderr, "Only %i filters supported. This includes undet and multi\n", NFILES_DS/2);
    return 1;
  }

  // Checking the input
  // Ifq is a mandatory argument
  if (Ifq == NULL || Ifq2 == NULL) {
    printHelpDialog_splitDS();
    fprintf(stderr, "Input *fq filenames were not properly initialized and \n");
    fprintf(stderr, "is a mandatory option. (--ifq <INPUT1.fq>:<INPUT2.fq>)\n");
    fprintf(stderr, "Exiting program\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  } else {
    fprintf(stderr, "- Fastq input file1: %s \n", Ifq);
    fprintf(stderr, "- Fastq input file2: %s \n", Ifq2);
  }

  if (nrfilters < 1) {
    printHelpDialog_splitDS();
    fprintf(stderr, "Specify at least 1 filter\n");
    return 1;
  }

  /**
   *  Loading the adapters file if the option is activated
   * 
   * */
  int Nad;
  DS_adap *adap_list = NULL;
  if (adapter_trim) {
    Fa_data *ad1 = malloc(sizeof(Fa_data));
    Fa_data *ad2 = malloc(sizeof(Fa_data));
    read_fasta(ad_fa, ad1);
    read_fasta(ad2_fa, ad2);
    Nad = ad1->nentries;
    par_TF.ad.Nad = ad1->nentries;
    adap_list = malloc(sizeof(DS_adap)*ad1->nentries);
    for (i = 0; i < ad1->nentries; i++) {
       adap_list[i] = init_DSadap(ad1->entry[i].seq, ad2->entry[i].seq,
                                  ad1->entry[i].N, ad2->entry[i].N);
    }
    init_alLUTs();
    init_map();
    free_fasta(ad1);
    free_fasta(ad2);
    fprintf(stderr, "- Adapters removal is activated!\n");
  }  // endif par_TF.is adapter

  /**
   * Load and initialize filters
   * 
   * */
  Tree **ptr_tree = calloc(sizeof(Tree *), nrfilters);
  Bfilter **ptr_bf = calloc(sizeof(Bfilter *), nrfilters);
  Bfkmer **ptr_bfkmer = calloc(sizeof(Bfkmer *), nrfilters);

  init_LUTs();
  init_map();
  char temp[MAX_FILENAME];
  for (i = 0; i < nrfilters; i++) {
    if (strsuff(index.s[i], ".bf")) {
      fprintf(stderr, "* DOING: Reading Bloom filter  from %s\n", index.s[i]);
      // Info filename
      snprintf(temp, MAX_FILENAME, "%s.txt", index.s[i]);
      ptr_bf[i] = read_Bfilter(index.s[i], temp);
      ptr_bfkmer[i] = init_Bfkmer(ptr_bf[i]->kmersize, ptr_bf[i]->hashNum);
    } else {
      fprintf(stderr, "* DOING: Reading tree structure from %s ... \n", index.s[i]);
      ptr_tree[i] = read_tree(index.s[i]);      
    }
  }

  /**
   * Keep stats of how many reads were written
   * 
   * */
  size_t *nrreads = calloc(sizeof(size_t), nrfilters + 2);

  /**
   * Output files (=filters+multi+undet)
   * 
   * */
  FILE **out_1 = calloc(sizeof(FILE *), nrfilters);
  FILE **out_2 = calloc(sizeof(FILE *), nrfilters);
  for (i = 0; i < nrfilters; i++) {
    char *suffix = NULL;
    if (!(suffix = strrchr(index.s[i], '/'))) {
      suffix = index.s[i];
    } else {
      suffix++;
    }
    snprintf(temp, MAX_FILENAME, "%s%s_1.fq.gz", prefix, suffix);
    out_1[i] = fopen_gen(temp, "w");
    snprintf(temp, MAX_FILENAME, "%s%s_2.fq.gz", prefix, suffix);
    out_2[i] = fopen_gen(temp, "w");
  }

  snprintf(temp, MAX_FILENAME, "%sundet_1.fq.gz", prefix);
  FILE *undet_1 = fopen_gen(temp, "w");
  snprintf(temp, MAX_FILENAME, "%sundet_2.fq.gz", prefix);
  FILE *undet_2 = fopen_gen(temp, "w");  

  snprintf(temp, MAX_FILENAME, "%smulti_1.fq.gz", prefix);
  FILE *multi_1 = fopen_gen(temp, "w");
  snprintf(temp, MAX_FILENAME, "%smulti_2.fq.gz", prefix);
  FILE *multi_2 = fopen_gen(temp, "w");  

  /**
   *  Opening fq file for reading
   * 
   * */
  FILE *fq_in1 = fopen_gen(Ifq, "r");
  FILE *fq_in2 = fopen_gen(Ifq2, "r");

  // Allocating memory for the fastq structure,
  Fq_read *seq1 = malloc(sizeof(Fq_read));
  Fq_read *seq2 = malloc(sizeof(Fq_read));

  int hit;
  int res;
  int hit_idx;
  int multi;
  int undet;
  int save_best = 0; // Save the best hit
  size_t counts = 0;
  char *char_seq1 = calloc(4*READ_MAXLEN, sizeof(char));  // string containing one fq read
  char *char_seq2 = calloc(4*READ_MAXLEN, sizeof(char));  // string containing one fq read
  int Nchar1, Nchar2;  // length of char_seq

  size_t nrundet = 0;
  size_t nrmulti = 0;
  size_t nrtrim = 0;
  size_t nrdisc = 0;
  size_t nrpolyG = 0;
  size_t nrlowQ = 0;

  int newl1 = 0, newl2 = 0;
  int offset1 = 0, offset2 = 0;
  int l1_i = 0, l1_f = 0, l2_i = 0, l2_f = 0;
  int j1 = 0, j2 = 0;
  int i_ad = 0;
  int nl1 = 0, nl2 = 0;
  int stop1 = 0, stop2 = 0;
  char *buffer1 = malloc(sizeof(char)*(B_LEN + 1));
  char *buffer2 = malloc(sizeof(char)*(B_LEN + 1));
  double score1;
  double score2;
  double curr_score;
  double best_score;
  do {
     newl1 = fread(buffer1+offset1, 1, B_LEN-offset1, fq_in1);
     newl2 = fread(buffer2+offset2, 1, B_LEN-offset2, fq_in2);
     newl1 += offset1;
     newl2 += offset2;
     buffer1[newl1] = '\0';
     buffer2[newl2] = '\0';
     j1 = 0; j2 = 0;
     while ((j1 <= newl1 && j2 <= newl2) && (newl1 ||  newl2)) {
        if ((buffer1[j1] == '\n') && (j1 < newl1)) {
          l1_f = j1;
          get_fqread(seq1, buffer1, l1_i, l1_f, nl1, L, 0);
          if ((nl1++ % 4 == 3)) {
             stop1 = 1;
             j1++;
          }
          l1_i = l1_f + 1;
        }
        if ((buffer2[j2] == '\n') && (j2 < newl2)) {
          l2_f = j2;
          get_fqread(seq2, buffer2, l2_i, l2_f, nl2, L, 0);
          if ((nl2++ % 4 == 3)) {
             stop2 = 1;
             j2++;
          }
          l2_i = l2_f + 1;
        }
        if (stop1 && !stop2) {
           j2++;
        } else if (!stop1 && stop2) {
           j1++;
        } else if (!stop1 && !stop2) {
           j1++;
           j2++;
        } else if (stop1 && stop2) {  // Do the stuff!!
          counts++;
          int trim = 0;
          int trim2 = 0;
          int discarded = 0;
          if (adapter_trim) {
            // omp_set_dynamic(0);     // Explicitly disable dynamic teams
            // omp_set_num_threads(2); // Use 4 threads for all consecutive parallel regions
            // #pragma omp parallel for private(position)
            for (i_ad=0; i_ad < Nad; i_ad++) {
              trim = trim_adapterDS(&adap_list[i_ad], seq1, seq2, DEFAULT_ZEROQ);
              // Too short or trimmed
              if (trim != 1)
                break;
            }
            // Read is too short
            if (trim == 0) {
              nrdisc++;
              // printf("disc\n");
              // Nchar1 = string_seq(seq1, char_seq1);
              // Nchar2 = string_seq(seq2, char_seq2);
              // printf("%s\n", char_seq1);
              discarded = 1;
            } else if (trim == 2) {
              nrtrim++;
              // printf("trim\n");
            }
          }
          if (min_poly_G && !discarded) {
            if (trim_polyX(seq1, 'G', min_poly_G) || trim_polyX(seq2, 'G', min_poly_G)) {
              nrpolyG++;
            }
            if (seq1->L < par_TF.minL || seq2->L < par_TF.minL) {
              nrdisc++;
              discarded = 1;
            }
          }
          if (par_TF.trimQ && !discarded) {
             trim = trim_sequenceQ(seq1);
             trim2 = trim_sequenceQ(seq2);
             discarded = (!trim) || (!trim2);
             if (discarded) {
                nrdisc++;
             } else if (trim == 2 || trim2 == 2) {
                nrlowQ++;
             }
           }
          if (!discarded) {
            hit = 0; // have a hit
            hit_idx = -1;
            undet = 1;
            multi = 0;
            for (i = 0; i < nrfilters; i++) {
              if (ptr_tree[i]) {
                res = (is_read_inTree(ptr_tree[i], seq1, &score1) || is_read_inTree(ptr_tree[i], seq2, &score2));
                if (res && !hit) {
                  hit = 1;
                  hit_idx = i;
                  undet = 0;
                  best_score = score1 + score2;
                } else if (res && hit) {
                  if (use_best) {
                    curr_score = score1 + score2;
                    if (curr_score > best_score) {
                      undet = 0;
                      hit_idx = i;
                      hit = 1;
                    }
                  } else {
                    multi = 1;
                    undet = 0;
                    hit = 0;
                    break;
                  }
                }
              } else if (ptr_bf[i]) {
                res =(is_read_inBloom(ptr_bf[i], seq1, ptr_bfkmer[i], &score1) || is_read_inBloom(ptr_bf[i], seq2, ptr_bfkmer[i], &score2));
                if (res && !hit) {
                  hit = 1;
                  hit_idx = i;
                  undet = 0;
                  best_score = score1 + score2;
                } else if (res && hit) {
                  if (use_best) {
                    curr_score = score1 + score2;
                    if (curr_score > best_score) {
                      undet = 0;
                      hit_idx = i;
                      hit = 1;
                    }
                  } else {
                    multi = 1;
                    undet = 0;
                    hit = 0;
                    break;
                  }
                }
              }
            }
            if (undet) {
              nrundet++;
              Nchar1 = string_seq(seq1, char_seq1);
              Nchar2 = string_seq(seq2, char_seq2);
              buffer_outputDS(undet_1, char_seq1, Nchar1, nrfilters*2);
              buffer_outputDS(undet_2, char_seq2, Nchar2, nrfilters*2+1);           
            } else if (multi) {
              nrmulti++;
              Nchar1 = string_seq(seq1, char_seq1);
              Nchar2 = string_seq(seq2, char_seq2);
              buffer_outputDS(multi_1, char_seq1, Nchar1, nrfilters*2+2);
              buffer_outputDS(multi_2, char_seq2, Nchar2, nrfilters*2+3);               
            } else if (hit) {
              nrreads[hit_idx]++;
              Nchar1 = string_seq(seq1, char_seq1);
              Nchar2 = string_seq(seq2, char_seq2);
              buffer_outputDS(out_1[hit_idx], char_seq1, Nchar1, hit_idx*2);
              buffer_outputDS(out_2[hit_idx], char_seq2, Nchar2, hit_idx*2+1);                
            }
          }
          if (counts % 1000000 == 0) {
            fprintf(stderr, "  %10ld reads have been read.\n", counts);
          }
          // reset stop
          stop1 = 0;
          stop2 = 0;
        }
     }  // end buffer loop
     offset1 = newl1 - l1_i;
     if (offset1 > -1)
       memmove(buffer1, buffer1 + l1_i, offset1);
     l1_f = -1;
     l1_i = 0;
     offset2 = newl2 - l2_i;
     if (offset2 > -1)
       memmove(buffer2, buffer2 + l2_i, offset2);
     l2_f = -1;
     l2_i = 0;
  }  while ((newl1 > offset1) ||  (newl2 > offset2));  // end read buffer

  // Check that the number of lines of both input files is the same
  if (nl1 != nl2) {
    fprintf(stderr, "ERROR: Input fq's contain different number of lines\n");
    fprintf(stderr, "%s contains %d lines \n", Ifq, nl1);
    fprintf(stderr, "%s contains %d lines \n", Ifq2, nl2);
    fprintf(stderr, "Exiting program\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    printf("stop1 %d, stop2 %d \n", stop1, stop2);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "- Number of lines in fq_files %d\n", nl1);
  // Printing the rest of the buffer outputs and closing file
  fprintf(stderr, "- Finished reading fq file.\n");
  fprintf(stderr, "- Closing files.\n");

  fclose(fq_in1);
  fclose(fq_in2);

  /**
   * Finish file output
   * 
   * */
  buffer_outputDS(undet_1, NULL, 0, nrfilters*2);
  fclose(undet_1);
  buffer_outputDS(undet_2, NULL, 0, nrfilters*2+1);
  fclose(undet_2);
  fprintf(stderr, "- Undetermined: %ld\n", nrundet);

  buffer_outputDS(multi_1, NULL, 0, nrfilters*2+2);
  fclose(multi_1);
  buffer_outputDS(multi_2, NULL, 0, nrfilters*2+3);
  fclose(multi_2);
  fprintf(stderr, "- Multimapped: %ld\n", nrmulti);

  if (adapter_trim) {
    fprintf(stderr, "- Too short: %ld\n", nrdisc);
    fprintf(stderr, "- Adapter trimmed: %ld\n", nrtrim);
  }

  for (i = 0; i < nrfilters; i++) {
    buffer_outputDS(out_1[i], NULL, 0, i*2);
    fclose(out_1[i]);
    buffer_outputDS(out_2[i], NULL, 0, i*2+1);
    fclose(out_2[i]);
    fprintf(stderr, "- %s: %ld\n", index.s[i], nrreads[i]);    
  }

  /**
   * Free file pointers
   * 
   * */
  free(out_1);
  free(out_2);

  /**
   * Free filters
   *
   * */
  for (i = 0; i < nrfilters; i++) {
    if (ptr_bf[i]) {
      free(ptr_bf[i]);
    }
    if (ptr_bfkmer[i]) {
      free(ptr_bfkmer[i]);
    }
    if (ptr_tree[i]) {
      free(ptr_tree[i]);
    }
  }
  free(ptr_bf);
  free(ptr_bfkmer);
  free(ptr_tree);

  free(adap_list);

  /**
   * Cleanup
   * 
   * */
  free(seq1);
  free(seq2);
  free(char_seq1);
  free(char_seq2);
  free(buffer1);
  free(buffer2);
  free(Ifq);
  free(Ifq2);

  timing_end(&t);
  
  // get to a pretty print version
  char* output = format_time_diff(&t);
  printf("Time: %s\n", output);
  free(output);

  return 0;
}
