# splitPE user manual

This program reads two paired end `fastq` files as input and splits/demultiplexes them 
based on matching to different contamination filters. It can optionally:
 - Trim/remove reads containing adapter remnants.
 - Trim low quality reads.
 - Trim reads containing polyG tails (common in NovaSeq data).

Reads are classified into separate output files based on which filter (index) they match. 
If a read matches multiple filters, it goes into a "multi" file. If it matches no filters, 
it goes into an "undetermined" file.

## Running the program

Usage `C` executable (in folder `bin`):

```
Usage: splitPE --ifq <INPUT1.fq>:<INPUT2.fq> --idx <INDEX_FILE>:[<INDEX_FILE2>:...]:<score>:<lmer_len>
                  --length <READ_LENGTH> --output [O_PREFIX] [options]

Reads in paired end fq files (gz, bz2, z formats also accepted) and splits them
based on contamination filters.

Options:
 -v, --version prints package version.
 -h, --help    prints help dialog.
 -f, --ifq     2 fastq input files [*fq|*fq.gz|*fq.bz2] separated by
               colons, mandatory option.
 -l, --length  read length: length of the reads, mandatory option.
 -o, --output  output prefix (with path), optional (default ./out).
 -b, --best    output to filter with the best score
 -x, --idx     index input file(s). One or more index files created by makeTree or
               makeBloom, followed by score and kmer length, all separated by colons.
               Format: <INDEX_FILE1>[:<INDEX_FILE2>:...]:<score>:<lmer_len>
               <INDEX_FILE>: output of makeTree or makeBloom (one or more files),
               <score>: score threshold to accept a match [0,1],
               <lmer_len>: the length of the lmers to be 
                        looked for in the reads [1,READ_LENGTH].
               Multiple index files allow splitting reads into different categories.
 -A, --adapter adapter input. Four fields separated by colons:
               <AD1.fa>: fasta file containing adapters,
               <AD2.fa>: fasta file containing adapters,
               <mismatches>: maximum mismatch count allowed,
               <score>: score threshold  for the aligner.
 -Q, --trimQ   NO:       does nothing to low quality reads (default),
               ALL:      removes all reads containing at least one low
                         quality nucleotide.
               ENDS:     trims the ends of the read if their quality is
                         below the threshold -q,
               FRAC:     discards a read if the fraction of bases with
                         low quality scores (below -q) is over 5 percent 
                         or a user defined percentage (-p). 
               ENDSFRAC: trims the ends and then discards the read if 
                         there are more low quality nucleotides than 
                         allowed by the option -p.
               GLOBAL:   removes n1 cycles on the left and n2 on the 
                         right, specified in -g.
               All reads are discarded if they are shorter than MINL
               (specified with -m or --minL).
 -g, --trimG   <min_len>: enable trimming of polyG repeats with a minimum length of min_len
                       at the end of read,
               All reads are discarded if they are shorter than MINL
               (specified with -m or --minL).
 -m, --minL    minimum length allowed for a read before it is discarded
               (default 36).
 -s, --stats   only report stats (no fastq output)
 -r, --reads   only test x reads (default ALL)
```

NOTE: the parameters -l or --length are meant to identify the length
of the reads in the input data. Actually, `splitPE` also copes with
data holding reads with different lengths. The length parameter must
hold the length of the longest read in the dataset.

## Output description

For each index file provided, splitPE creates two output files (one for each read):
- `[O_PREFIX]<INDEX_FILENAME>_1.fq.gz`: contains reads from read1 that matched this filter.
- `[O_PREFIX]<INDEX_FILENAME>_2.fq.gz`: contains reads from read2 that matched this filter.

Additionally, two special categories are created:
- `[O_PREFIX]undet_1.fq.gz` and `[O_PREFIX]undet_2.fq.gz`: contains reads that did not match any filter (undetermined).
- `[O_PREFIX]multi_1.fq.gz` and `[O_PREFIX]multi_2.fq.gz`: contains reads that matched multiple filters (unless `--best` is used).

When the `--best` option is enabled, reads matching multiple filters will be assigned to the filter with the highest score instead of being placed in the "multi" files.

## Filters

### Index Files (Contamination Detection)

Contamination detection is performed using index files created by either `makeTree` or `makeBloom`.
Multiple index files can be provided separated by colons to split reads into different categories.

The tool supports two methods:
- **TREE**: Index files created by `makeTree` (suffix varies). Suitable for small reference sequences.
- **BLOOM**: Index files created by `makeBloom` (suffix `.bf`). Suitable for large reference sequences.

Each read pair is checked against all provided index files. The score for matching is computed as 
the proportion of k-mers found in the reference. If either read in the pair matches a filter 
(score above threshold), the pair is assigned to that filter's output.

See `README_trimFilter.md` for more details on how contaminations are detected using TREE and BLOOM methods.

### Adapters

Technical sequences within the reads are detected if the option
`--adapter <ADAPTERS1.fa>:<ADAPTERS2.fa>:<mismatches>:<score>` is given.

The adapter sequences are read from the fasta files and prepended to their respective reads. 
A 'seed and extend' approach is used to look for overlaps. Reads are trimmed to remove 
adapter sequences. If the remaining read is shorter than `minL`, the read pair is discarded.

See `README_trimFilterPE.md` for detailed information on how adapters are detected and trimmed 
in paired-end mode.

### Low quality trimming

Low quality trimming follows the same procedure as for `trimFilterPE`. 
We list the options below; see `README_trimFilter.md` for more details.

- `--trimQ NO` or flag absent: nothing is done to reads with low quality.
- `--trimQ ALL`: all reads containing at least one low quality nucleotide are discarded.
- `--trimQ ENDS`: look for low quality (below minQ) base callings at the
  beginning and at the end of the read. Trim them at both ends until the
  quality is above the threshold. Discard the read if the remaining length is smaller than `minL`.
- `--trimQ FRAC [--percent p]`: discard the read if there are more than `p%` 
  nucleotides whose quality lies below the threshold (default `-p 5`).
- `--trimQ ENDSFRAC --percent p`: first trim the ends as in the `ENDS` option,
  then discard if the number of low quality nucleotides exceeds `p%` (default `-p 5`).
- `--trimQ GLOBAL --global n1:n2`: cut all reads globally `n1` nucleotides from
  the left and `n2` from the right.

**Note:** qualities are evaluated assuming the reads follow the
L - Illumina 1.8+ Phred+33 convention, see [Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format#Encoding).
Adjust the values for a different convention.

### PolyG trimming

NovaSeq sequencers can produce artificial polyG tails when there is insufficient signal. 
The `--trimG` option allows trimming of these polyG sequences from the 3' end of reads.

Usage: `--trimG <min_len>`

Where `min_len` is the minimum length of a polyG sequence to trigger trimming. 
The algorithm allows for some mismatches within the polyG tail (1 mismatch per 8 bp, 
maximum 5 mismatches total) to handle sequencing errors within the polyG region.

If trimming results in reads shorter than `minL`, the read pair is discarded.

## Statistics

The program outputs statistics to stderr showing:
- Proportion of undetermined reads (did not match any filter)
- Proportion of multimapped reads (matched multiple filters, if `--best` is not used)
- Proportion of reads that were too short after trimming
- Proportion of reads with adapter trimming (if enabled)
- Proportion of reads with polyG trimming (if enabled)
- Proportion of reads with quality trimming (if enabled)
- Proportion of reads matching each filter

When the `--stats` option is used, only statistics are reported and no output fastq files are created.

## Example Usage

### Basic splitting with contamination filters

Split paired-end reads based on two reference sequences (e.g., different species):

```bash
# Create bloom filters for two organisms
makeBloom --fasta organism1.fa --output org1 --kmersize 25 --fal_pos_rate 0.01
makeBloom --fasta organism2.fa --output org2 --kmersize 25 --fal_pos_rate 0.01

# Split reads
splitPE --ifq sample_R1.fq.gz:sample_R2.fq.gz \
        --idx org1.bf:org2.bf:0.15:25 \
        --length 150 \
        --output output_prefix_
```

This will create:
- `output_prefix_org1.bf_1.fq.gz` and `output_prefix_org1.bf_2.fq.gz` (organism 1 reads)
- `output_prefix_org2.bf_1.fq.gz` and `output_prefix_org2.bf_2.fq.gz` (organism 2 reads)
- `output_prefix_undet_1.fq.gz` and `output_prefix_undet_2.fq.gz` (undetermined reads)
- `output_prefix_multi_1.fq.gz` and `output_prefix_multi_2.fq.gz` (reads matching both)

### With adapter and quality trimming

```bash
splitPE --ifq sample_R1.fq.gz:sample_R2.fq.gz \
        --idx org1.bf:org2.bf:0.15:25 \
        --length 150 \
        --adapter adapters_R1.fa:adapters_R2.fa:2:8.0 \
        --trimQ ENDSFRAC \
        --minL 36 \
        --output output_prefix_
```

### With polyG trimming (NovaSeq data)

```bash
splitPE --ifq sample_R1.fq.gz:sample_R2.fq.gz \
        --idx org1.bf:org2.bf:0.15:25 \
        --length 150 \
        --trimG 10 \
        --minL 36 \
        --output output_prefix_
```

### Best match mode

Use `--best` to assign reads to the filter with the highest score when multiple matches occur:

```bash
splitPE --ifq sample_R1.fq.gz:sample_R2.fq.gz \
        --idx org1.bf:org2.bf:0.15:25 \
        --length 150 \
        --best \
        --output output_prefix_
```

In this mode, the `multi` files will be empty as reads are assigned to their best match.

### Statistics only

To check classification statistics without creating output files:

```bash
splitPE --ifq sample_R1.fq.gz:sample_R2.fq.gz \
        --idx org1.bf:org2.bf:0.15:25 \
        --length 150 \
        --stats
```

## Contributors

Gerben Voshol, Paula Pérez Rubio

## License

GPL v3 (see LICENSE.txt)
