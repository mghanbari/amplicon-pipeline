---
kernelspec:
  display_name: R
  language: R
  name: ir
language_info:
  codemirror_mode: r
  file_extension: .r
  mimetype: text/x-r-source
  name: R
  pygments_lexer: r
  version: 3.2.4
---

# Using DADA2 for QIIME using the Moving Pictures of the Human Microbiome
This tutorial covers processing raw sequence data with [DADA2](https://github.com/benjjneb/dada2) for use with QIIME using Illumina sequencing data. This tutorial is intended to be a quick to run, and as such, uses only a subset of a full Illumina Genome Analyzer II (GAIIx) run. This tutorial is only intended to illustrate the workflow necessary to generate an analagous table to the OTU table typically used in QIIME analyses. For a more in depth use of QIIME please refer to the [other tutorials](http://qiime.org/tutorials/#). Additionally it is recommended that you first follow the official [DADA2 tutorial](http://benjjneb.github.io/dada2/tutorial.html) prior to running this tutorial as this tutorial is intended to provide only the minimum workflow commands to create a biom table for use with QIIME.

The data used in this tutorial are derived from the Moving Pictures of the Human Microbiome study, where two human subjects collected daily samples from four body sites: the tongue, the palm of the left hand, the palm of the right hand, and the gut (via fecal samples obtained by swapping used toilet paper). These data were sequenced using the barcoded amplicon sequencing protocol described in Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. A more recent version of this protocol that can be used with the Illumina HiSeq 2000 and MiSeq can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/22402401).

This tutorial is presented as a Jupyter Notebook. You can find more information on the Jupyter Notebook [here](http://jupyter.org/).

## Getting started
We'll begin by downloading the tutorial data. Commands like the following must be run from the bash terminal.

```bash
wget ftp://ftp.microbio.me/qiime/tutorial_files/moving_pictures_tutorial-1.9.0.tgz || curl -O ftp://ftp.microbio.me/qiime/tutorial_files/moving_pictures_tutorial-1.9.0.tgz  

tar -xzf moving_pictures_tutorial-1.9.0.tgz
```

### Demultiplexing and quality filtering sequences
We next need to demultiplex (i.e. assign barcoded reads to the samples they are derived from). In general, you will get separate fastq files for your sequence and barcode reads. Note that we pass these files while still gzipped. split_libraries_fastq.py can handle gzipped or unzipped fastq files.  
  
The difference in processing the fastq files for use with DADA2 as opposed to the OTU picking pipeline in QIIME is that we will not be quality filtering the files here, but rather later in the DADA2 pipeline. The default strategy in QIIME for quality filtering of Illumina data is described in [Bokulich et al (2013)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531572/).

We will run split libraries with the quality filtering effectively turned off. Additionally we will pass the `--store_demultiplexed_fastq` flag as we want to retain the quality information stored in the fastq files for later use.

Note: If you received the sequence data where each sample is represented in a unique fastq file you should skip to the DADA2 pipeline commands

```bash
cd moving_pictures_tutorial-1.9.0  

split_libraries_fastq.py -o slout/ -i illumina/forward_reads.fastq.gz -b illumina/barcodes.fastq.gz -m illumina/map.tsv --store_demultiplexed_fastq -r 999 -n 999 -q 0 -p 0.000001
```

### Split sequence files on a per sample basis
DADA2 requires a unique fastq file for each sample, the following command will split the fulll fastq file on a per sample basis

```bash
split_sequence_file_on_sample_ids.py -i slout/seqs.fastq -o sample_fastqs --file_type fastq 
```

### Start the DADA2 pipeline
At this point you are ready to start the DADA2 pipeline. DADA2 is available only with an R api. This can be run from the [R](https://www.r-project.org/)m terminal, [R-studio](https://www.rstudio.com/First), or from with within the Jupyter notebook environment with [IRkernel](http://irkernel.github.io/). More details can be found in the official [DADA2 tutorial](http://benjjneb.github.io/dada2/tutorial.html)

### Load the necessary libraries. 
If you don’t already have the dada2 package, see the [dada2 installation instructions](http://benjjneb.github.io/dada2/dada-installation.html). The ShortRead package is available from Bioconductor, and ggplot2 from CRAN or Bioconductor:

```python
>>> suppressPackageStartupMessages(library(dada2)); packageVersion("dada2")
>>> suppressPackageStartupMessages(library(ShortRead)); packageVersion("ShortRead")
>>> suppressPackageStartupMessages(library(ggplot2)); packageVersion("ggplot2")
Warning message:
: package ‘Rcpp’ was built under R version 3.2.5Creating a generic function for ‘nchar’ from package ‘base’ in package ‘S4Vectors’
[1] ‘1.1.1’
[1] ‘1.26.0’
[1] ‘2.1.0’
```

### Filtering and Trimming
First we read in the file names for all the fastq files and do a little string manipulation to get a list of the fastq files

```python
>>> path <- "/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/"#replace this with the path to your sample fastqs
>>> fns <- list.files(path)
>>> sample.names <- sapply(strsplit(fns, ".f"), `[`, 1)#get list of sample names
>>> fastqs <- paste0(path, fns)#create list of filepaths to each fastq
```

### Examine quality profiles of forward and reverse reads
It is always important to look at your data. We start by visualizing the quality profiles along the sequencing reads.

Visualize the quality profile of the forward reads:

```python
>>> plotQualityProfile(fastqs[[1]])

plot without title
```

### Perform filtering and trimming
Based on the quality score plots we will trim the reads at 130 where the read quality really starts to drop off

```python
>>> # Make filenames for the outputed filtered fastq files
... filts <- paste0(path, sample.names, "_filt.fastq.gz")
...
>>> for(i in seq_along(fastqs)) {
...   fastqFilter(fastqs[i], filts[i],
...               trimLeft=10, truncLen=c(130),#trim the first 10 and last 20 base pairs from each sequence
...               maxN=0, maxEE=2, truncQ=2,
...               compress=TRUE, verbose=TRUE)
>>> }
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S105_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 6843, output 0 filtered sequences.
Read in 11340, output 6843 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S140_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5616, output 0 filtered sequences.
Read in 9736, output 5616 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S208_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 7117, output 0 filtered sequences.
Read in 11335, output 7117 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S257_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5132, output 0 filtered sequences.
Read in 8216, output 5132 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S281_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5370, output 0 filtered sequences.
Read in 8907, output 5370 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S57_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 7120, output 0 filtered sequences.
Read in 11752, output 7120 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S76_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 6493, output 0 filtered sequences.
Read in 10100, output 6493 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S8_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5700, output 0 filtered sequences.
Read in 12388, output 5700 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S155_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 2327, output 0 filtered sequences.
Read in 9262, output 2327 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S175_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 2698, output 0 filtered sequences.
Read in 10692, output 2698 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S204_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1914, output 0 filtered sequences.
Read in 7297, output 1914 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S222_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 2001, output 0 filtered sequences.
Read in 8386, output 2001 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S240_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 3213, output 0 filtered sequences.
Read in 11986, output 3213 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S309_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 935, output 0 filtered sequences.
Read in 3836, output 935 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S357_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1502, output 0 filtered sequences.
Read in 6254, output 1502 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S382_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 2700, output 0 filtered sequences.
Read in 9301, output 2700 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S242_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 651, output 0 filtered sequences.
Read in 2055, output 651 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S294_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 846, output 0 filtered sequences.
Read in 2318, output 846 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S313_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 810, output 0 filtered sequences.
Read in 1996, output 810 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S341_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 689, output 0 filtered sequences.
Read in 1854, output 689 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S360_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 685, output 0 filtered sequences.
Read in 2085, output 685 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S378_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 856, output 0 filtered sequences.
Read in 2147, output 856 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S112_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5559, output 0 filtered sequences.
Read in 16266, output 5559 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S137_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 6449, output 0 filtered sequences.
Read in 18787, output 6449 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S63_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 6952, output 0 filtered sequences.
Read in 17168, output 6952 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S104_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1698, output 0 filtered sequences.
Read in 3460, output 1698 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S155_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1489, output 0 filtered sequences.
Read in 2492, output 1489 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S174_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1439, output 0 filtered sequences.
Read in 2792, output 1439 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S203_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1722, output 0 filtered sequences.
Read in 3037, output 1722 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S222_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 2062, output 0 filtered sequences.
Read in 3394, output 2062 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S240_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 1523, output 0 filtered sequences.
Read in 2667, output 1523 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S20_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5768, output 0 filtered sequences.
Read in 9776, output 5768 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S68_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 4792, output 0 filtered sequences.
Read in 9557, output 4792 filtered sequences.
Overwriting file:/Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S93_filt.fastq.gz
Warning message:
In min(qs, na.rm = TRUE): no non-missing arguments to min; returning InfRead in 5434, output 0 filtered sequences.
Read in 11270, output 5434 filtered sequences.
```

### Dereplication
In the dereplication step, all reads with identical sequences are combined into “unique sequences” with a corresponding abundance, i.e. the number of reads with that same sequence.

DADA2 retains a summary of the quality information associated with each unique sequence. DADA2 constructs a “consensus” quality profile for each unique sequence by averaging the positional qualities from the dereplicated reads. These consensus quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.

```python
>>> dereps <- derepFastq(filts, verbose=TRUE)
>>> # Name the derep-class objects by the sample names
... names(dereps) <- sample.names
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S105_filt.fastq.gz
Encountered 1689 unique sequences from 6843 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S105_filt.fastq.gz
Encountered 1689 unique sequences from 6843 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S140_filt.fastq.gz
Encountered 1330 unique sequences from 5616 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S140_filt.fastq.gz
Encountered 1330 unique sequences from 5616 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S208_filt.fastq.gz
Encountered 1971 unique sequences from 7117 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S208_filt.fastq.gz
Encountered 1971 unique sequences from 7117 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S257_filt.fastq.gz
Encountered 1528 unique sequences from 5132 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S257_filt.fastq.gz
Encountered 1528 unique sequences from 5132 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S281_filt.fastq.gz
Encountered 1624 unique sequences from 5370 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S281_filt.fastq.gz
Encountered 1624 unique sequences from 5370 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S57_filt.fastq.gz
Encountered 1767 unique sequences from 7120 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S57_filt.fastq.gz
Encountered 1767 unique sequences from 7120 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S76_filt.fastq.gz
Encountered 1501 unique sequences from 6493 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S76_filt.fastq.gz
Encountered 1501 unique sequences from 6493 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S8_filt.fastq.gz
Encountered 1230 unique sequences from 5700 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L1S8_filt.fastq.gz
Encountered 1230 unique sequences from 5700 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S155_filt.fastq.gz
Encountered 813 unique sequences from 2327 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S155_filt.fastq.gz
Encountered 813 unique sequences from 2327 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S175_filt.fastq.gz
Encountered 858 unique sequences from 2698 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S175_filt.fastq.gz
Encountered 858 unique sequences from 2698 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S204_filt.fastq.gz
Encountered 751 unique sequences from 1914 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S204_filt.fastq.gz
Encountered 751 unique sequences from 1914 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S222_filt.fastq.gz
Encountered 992 unique sequences from 2001 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S222_filt.fastq.gz
Encountered 992 unique sequences from 2001 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S240_filt.fastq.gz
Encountered 735 unique sequences from 3213 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S240_filt.fastq.gz
Encountered 735 unique sequences from 3213 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S309_filt.fastq.gz
Encountered 407 unique sequences from 935 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S309_filt.fastq.gz
Encountered 407 unique sequences from 935 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S357_filt.fastq.gz
Encountered 580 unique sequences from 1502 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S357_filt.fastq.gz
Encountered 580 unique sequences from 1502 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S382_filt.fastq.gz
Encountered 887 unique sequences from 2700 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L2S382_filt.fastq.gz
Encountered 887 unique sequences from 2700 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S242_filt.fastq.gz
Encountered 171 unique sequences from 651 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S242_filt.fastq.gz
Encountered 171 unique sequences from 651 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S294_filt.fastq.gz
Encountered 349 unique sequences from 846 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S294_filt.fastq.gz
Encountered 349 unique sequences from 846 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S313_filt.fastq.gz
Encountered 343 unique sequences from 810 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S313_filt.fastq.gz
Encountered 343 unique sequences from 810 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S341_filt.fastq.gz
Encountered 327 unique sequences from 689 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S341_filt.fastq.gz
Encountered 327 unique sequences from 689 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S360_filt.fastq.gz
Encountered 413 unique sequences from 685 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S360_filt.fastq.gz
Encountered 413 unique sequences from 685 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S378_filt.fastq.gz
Encountered 255 unique sequences from 856 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L3S378_filt.fastq.gz
Encountered 255 unique sequences from 856 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S112_filt.fastq.gz
Encountered 1542 unique sequences from 5559 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S112_filt.fastq.gz
Encountered 1542 unique sequences from 5559 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S137_filt.fastq.gz
Encountered 1503 unique sequences from 6449 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S137_filt.fastq.gz
Encountered 1503 unique sequences from 6449 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S63_filt.fastq.gz
Encountered 2126 unique sequences from 6952 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L4S63_filt.fastq.gz
Encountered 2126 unique sequences from 6952 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S104_filt.fastq.gz
Encountered 369 unique sequences from 1698 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S104_filt.fastq.gz
Encountered 369 unique sequences from 1698 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S155_filt.fastq.gz
Encountered 320 unique sequences from 1489 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S155_filt.fastq.gz
Encountered 320 unique sequences from 1489 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S174_filt.fastq.gz
Encountered 309 unique sequences from 1439 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S174_filt.fastq.gz
Encountered 309 unique sequences from 1439 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S203_filt.fastq.gz
Encountered 385 unique sequences from 1722 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S203_filt.fastq.gz
Encountered 385 unique sequences from 1722 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S222_filt.fastq.gz
Encountered 397 unique sequences from 2062 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S222_filt.fastq.gz
Encountered 397 unique sequences from 2062 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S240_filt.fastq.gz
Encountered 337 unique sequences from 1523 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L5S240_filt.fastq.gz
Encountered 337 unique sequences from 1523 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S20_filt.fastq.gz
Encountered 814 unique sequences from 5768 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S20_filt.fastq.gz
Encountered 814 unique sequences from 5768 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S68_filt.fastq.gz
Encountered 791 unique sequences from 4792 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S68_filt.fastq.gz
Encountered 791 unique sequences from 4792 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S93_filt.fastq.gz
Encountered 990 unique sequences from 5434 total sequences read.
Dereplicating sequence entries in Fastq file: /Users/jc33/dev/dada2_tutorial/moving_pictures_tutorial-1.9.0/sample_fastqs/L6S93_filt.fastq.gz
Encountered 990 unique sequences from 5434 total sequences read.
```

### Sample Inference
We are now ready to apply DADA2’s core sample inference algorithm to the dereplicated sequences.

```python
>>> dada_seqs <- dada(dereps, err=inflateErr(tperr1,3), selfConsist = TRUE)
Sample 1 - 6843 reads in 1689 unique sequences.
Sample 2 - 6843 reads in 1689 unique sequences.
Sample 3 - 5616 reads in 1330 unique sequences.
Sample 4 - 5616 reads in 1330 unique sequences.
Sample 5 - 7117 reads in 1971 unique sequences.
Sample 6 - 7117 reads in 1971 unique sequences.
Sample 7 - 5132 reads in 1528 unique sequences.
Sample 8 - 5132 reads in 1528 unique sequences.
Sample 9 - 5370 reads in 1624 unique sequences.
Sample 10 - 5370 reads in 1624 unique sequences.
Sample 11 - 7120 reads in 1767 unique sequences.
Sample 12 - 7120 reads in 1767 unique sequences.
Sample 13 - 6493 reads in 1501 unique sequences.
Sample 14 - 6493 reads in 1501 unique sequences.
Sample 15 - 5700 reads in 1230 unique sequences.
Sample 16 - 5700 reads in 1230 unique sequences.
Sample 17 - 2327 reads in 813 unique sequences.
Sample 18 - 2327 reads in 813 unique sequences.
Sample 19 - 2698 reads in 858 unique sequences.
Sample 20 - 2698 reads in 858 unique sequences.
Sample 21 - 1914 reads in 751 unique sequences.
Sample 22 - 1914 reads in 751 unique sequences.
Sample 23 - 2001 reads in 992 unique sequences.
Sample 24 - 2001 reads in 992 unique sequences.
Sample 25 - 3213 reads in 735 unique sequences.
Sample 26 - 3213 reads in 735 unique sequences.
Sample 27 - 935 reads in 407 unique sequences.
Sample 28 - 935 reads in 407 unique sequences.
Sample 29 - 1502 reads in 580 unique sequences.
Sample 30 - 1502 reads in 580 unique sequences.
Sample 31 - 2700 reads in 887 unique sequences.
Sample 32 - 2700 reads in 887 unique sequences.
Sample 33 - 651 reads in 171 unique sequences.
Sample 34 - 651 reads in 171 unique sequences.
Sample 35 - 846 reads in 349 unique sequences.
Sample 36 - 846 reads in 349 unique sequences.
Sample 37 - 810 reads in 343 unique sequences.
Sample 38 - 810 reads in 343 unique sequences.
Sample 39 - 689 reads in 327 unique sequences.
Sample 40 - 689 reads in 327 unique sequences.
Sample 41 - 685 reads in 413 unique sequences.
Sample 42 - 685 reads in 413 unique sequences.
Sample 43 - 856 reads in 255 unique sequences.
Sample 44 - 856 reads in 255 unique sequences.
Sample 45 - 5559 reads in 1542 unique sequences.
Sample 46 - 5559 reads in 1542 unique sequences.
Sample 47 - 6449 reads in 1503 unique sequences.
Sample 48 - 6449 reads in 1503 unique sequences.
Sample 49 - 6952 reads in 2126 unique sequences.
Sample 50 - 6952 reads in 2126 unique sequences.
Sample 51 - 1698 reads in 369 unique sequences.
Sample 52 - 1698 reads in 369 unique sequences.
Sample 53 - 1489 reads in 320 unique sequences.
Sample 54 - 1489 reads in 320 unique sequences.
Sample 55 - 1439 reads in 309 unique sequences.
Sample 56 - 1439 reads in 309 unique sequences.
Sample 57 - 1722 reads in 385 unique sequences.
Sample 58 - 1722 reads in 385 unique sequences.
Sample 59 - 2062 reads in 397 unique sequences.
Sample 60 - 2062 reads in 397 unique sequences.
Sample 61 - 1523 reads in 337 unique sequences.
Sample 62 - 1523 reads in 337 unique sequences.
Sample 63 - 5768 reads in 814 unique sequences.
Sample 64 - 5768 reads in 814 unique sequences.
Sample 65 - 4792 reads in 791 unique sequences.
Sample 66 - 4792 reads in 791 unique sequences.
Sample 67 - 5434 reads in 990 unique sequences.
Sample 68 - 5434 reads in 990 unique sequences.
   selfConsist step 2 
   selfConsist step 3 
   selfConsist step 4 


Convergence after  4  rounds.
```

### Create a table of sequences per sample

```python
>>> seqtab <- makeSequenceTable(dada_seqs)
```

### Create a fasta file of the sequences present in the table
This is necessary for many qiime scripts and is analogous to the `rep_set.fna` output by QIIME's OTU picking pipelines

```python
>>> uniquesToFasta(dada_seqs, "moving_pictures_tutorial-1.9.0/illumina/slout/rep_set.fna")
```

### Write the resulting table to file
In order to convert the resulting table to a biom table for use with QIIME, we need to format the output so that [biom](http://biom-format.org/) is able to convert it from text to biom format

```python
>>> seqtab <- t(seqtab)#transpose the table
>>> seqtab <- cbind('#OTUID' = rownames(seqtab), seqtab)#Add '#OTUID' to the header (required by biom)
>>> write.table(seqtab, "moving_pictures_tutorial-1.9.0/dada2_seq_table.txt", sep='\t', row.names=FALSE, quote=FALSE)
```

### Convert the biom table to hdf5 format

```bash
biom convert -i dada2_seq_table.txt -o seq_table.biom --table-type "OTU table" --to-hdf5
```

### Find addtional tutorials [here](http://qiime.org/tutorials/#)
At this point the core files for QIIME's diversity analyses scripts are ready to go.
