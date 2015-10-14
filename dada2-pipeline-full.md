Workflow to create an allele table for samples
=========================

```python
>>> library(dada2); packageVersion("dada2")
>>> library(ShortRead); packageVersion("ShortRead")
>>> library(ggplot2); packageVersion("ggplot2")
Loading required package: Rcpp
Creating a generic function for ‘nchar’ from package ‘base’ in package ‘S4Vectors’

Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: ‘BiocGenerics’

The following objects are masked from ‘package:parallel’:

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following object is masked from ‘package:stats’:

    xtabs

The following objects are masked from ‘package:base’:

    anyDuplicated, append, as.data.frame, as.vector, cbind, colnames,
    do.call, duplicated, eval, evalq, Filter, Find, get, intersect,
    is.unsorted, lapply, Map, mapply, match, mget, order, paste, pmax,
    pmax.int, pmin, pmin.int, Position, rank, rbind, Reduce, rep.int,
    rownames, sapply, setdiff, sort, table, tapply, union, unique,
    unlist, unsplit

Loading required package: BiocParallel
Loading required package: Biostrings
Loading required package: S4Vectors
Loading required package: stats4
Loading required package: IRanges
Loading required package: XVector
Loading required package: Rsamtools
Loading required package: GenomeInfoDb
Loading required package: GenomicRanges
Loading required package: GenomicAlignments
[1] ‘0.9.5’

[1] ‘1.26.0’
[1] ‘1.0.1’
```

```python
>>> path <- "AmpliconSeq/"
>>> amplicons <- list.files(path)
```

```python
>>> length(amplicons)
[1] 170
```

```python
>>> fastqs <- amplicons[grepl(".fastq", amplicons)]
>>> fastqs <- sort(fastqs) # Sort should keep them paired in order
>>> fnFs <- fastqs[grepl("_R1", fastqs)]
>>> fnRs <- fastqs[grepl("_R2", fastqs)]
```

```python
>>> filtFs <- paste0(path, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq.gz")
>>> filtRs <- paste0(path, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq.gz")
>>> for(i in seq_along(fnFs)) {
...     print(i)
...     print(fnFs[i])
...     fastqPairedFilter(paste0(path, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]), maxN=0, maxEE=2, truncQ=2, trimLeft=c(0, 22), truncLen=c(240,160), compress=TRUE, verbose=TRUE)
>>> }
[1] 1
[1] "C-gunnisoni-BLS01-PDMHC-Amplicon_S77_L001_R1_001.fastq.gz"
Read in 1950 paired-sequences, outputted 1273 filtered paired-sequences.
[1] 2
[1] "C-gunnisoni-BLS02-PDMHC-Amplicon_S76_L001_R1_001.fastq.gz"
Read in 1541 paired-sequences, outputted 947 filtered paired-sequences.
[1] 3
[1] "C-gunnisoni-BLS03-PDMHC-Amplicon_S75_L001_R1_001.fastq.gz"
Read in 3159 paired-sequences, outputted 2039 filtered paired-sequences.
[1] 4
[1] "C-gunnisoni-BLS04-PDMHC-Amplicon_S74_L001_R1_001.fastq.gz"
Read in 3201 paired-sequences, outputted 2051 filtered paired-sequences.
[1] 5
[1] "C-gunnisoni-BLS05-PDMHC-Amplicon_S73_L001_R1_001.fastq.gz"
Read in 5634 paired-sequences, outputted 3347 filtered paired-sequences.
[1] 6
[1] "C-gunnisoni-BLS06-PDMHC-Amplicon_S72_L001_R1_001.fastq.gz"
Read in 947 paired-sequences, outputted 585 filtered paired-sequences.
[1] 7
[1] "C-gunnisoni-CRL01-PDMHC-Amplicon_S54_L001_R1_001.fastq.gz"
Read in 6251 paired-sequences, outputted 5292 filtered paired-sequences.
[1] 8
[1] "C-gunnisoni-CRL02-PDMHC-Amplicon_S53_L001_R1_001.fastq.gz"
Read in 1530 paired-sequences, outputted 1200 filtered paired-sequences.
[1] 9
[1] "C-gunnisoni-CRL03-PDMHC-Amplicon_S52_L001_R1_001.fastq.gz"
Read in 2491 paired-sequences, outputted 2020 filtered paired-sequences.
[1] 10
[1] "C-gunnisoni-CRL04-PDMHC-Amplicon_S51_L001_R1_001.fastq.gz"
Read in 1974 paired-sequences, outputted 1699 filtered paired-sequences.
[1] 11
[1] "C-gunnisoni-CRL05-PDMHC-Amplicon_S50_L001_R1_001.fastq.gz"
Read in 2785 paired-sequences, outputted 2413 filtered paired-sequences.
[1] 12
[1] "C-gunnisoni-CRL06-PDMHC-Amplicon_S49_L001_R1_001.fastq.gz"
Read in 2472 paired-sequences, outputted 2140 filtered paired-sequences.
[1] 13
[1] "C-gunnisoni-CRL07-PDMHC-Amplicon_S48_L001_R1_001.fastq.gz"
Read in 2621 paired-sequences, outputted 2291 filtered paired-sequences.
[1] 14
[1] "C-gunnisoni-CRL08-PDMHC-Amplicon_S47_L001_R1_001.fastq.gz"
Read in 3045 paired-sequences, outputted 2651 filtered paired-sequences.
[1] 15
[1] "C-gunnisoni-CRL09-PDMHC-Amplicon_S62_L001_R1_001.fastq.gz"
Read in 4009 paired-sequences, outputted 3407 filtered paired-sequences.
[1] 16
[1] "C-gunnisoni-CRL10-PDMHC-Amplicon_S61_L001_R1_001.fastq.gz"
Read in 3118 paired-sequences, outputted 2634 filtered paired-sequences.
[1] 17
[1] "C-gunnisoni-CRL11-PDMHC-Amplicon_S60_L001_R1_001.fastq.gz"
Read in 3718 paired-sequences, outputted 3159 filtered paired-sequences.
[1] 18
[1] "C-gunnisoni-CRL12-PDMHC-Amplicon_S59_L001_R1_001.fastq.gz"
Read in 4220 paired-sequences, outputted 3530 filtered paired-sequences.
[1] 19
[1] "C-gunnisoni-CRL13-PDMHC-Amplicon_S58_L001_R1_001.fastq.gz"
Read in 6537 paired-sequences, outputted 5693 filtered paired-sequences.
[1] 20
[1] "C-gunnisoni-CRL14-PDMHC-Amplicon_S57_L001_R1_001.fastq.gz"
Read in 1483 paired-sequences, outputted 1233 filtered paired-sequences.
[1] 21
[1] "C-gunnisoni-CRL15-PDMHC-Amplicon_S56_L001_R1_001.fastq.gz"
Read in 4299 paired-sequences, outputted 3728 filtered paired-sequences.
[1] 22
[1] "C-gunnisoni-CRL16-PDMHC-Amplicon_S55_L001_R1_001.fastq.gz"
Read in 6097 paired-sequences, outputted 5412 filtered paired-sequences.
[1] 23
[1] "C-gunnisoni-CRL17-PDMHC-Amplicon_S70_L001_R1_001.fastq.gz"
Read in 4340 paired-sequences, outputted 3671 filtered paired-sequences.
[1] 24
[1] "C-gunnisoni-CRL18-PDMHC-Amplicon_S69_L001_R1_001.fastq.gz"
Read in 974 paired-sequences, outputted 836 filtered paired-sequences.
[1] 25
[1] "C-gunnisoni-CRL19-PDMHC-Amplicon_S68_L001_R1_001.fastq.gz"
Read in 3472 paired-sequences, outputted 2997 filtered paired-sequences.
[1] 26
[1] "C-gunnisoni-CRL20-PDMHC-Amplicon_S67_L001_R1_001.fastq.gz"
Read in 4968 paired-sequences, outputted 4283 filtered paired-sequences.
[1] 27
[1] "C-gunnisoni-CRL21-PDMHC-Amplicon_S66_L001_R1_001.fastq.gz"
Read in 1489 paired-sequences, outputted 1238 filtered paired-sequences.
[1] 28
[1] "C-gunnisoni-CRL22-PDMHC-Amplicon_S65_L001_R1_001.fastq.gz"
Read in 3221 paired-sequences, outputted 2821 filtered paired-sequences.
[1] 29
[1] "C-gunnisoni-CRL23-PDMHC-Amplicon_S64_L001_R1_001.fastq.gz"
Read in 2846 paired-sequences, outputted 2495 filtered paired-sequences.
[1] 30
[1] "C-gunnisoni-CRL24-PDMHC-Amplicon_S63_L001_R1_001.fastq.gz"
Read in 2103 paired-sequences, outputted 1808 filtered paired-sequences.
[1] 31
[1] "C-gunnisoni-ELMA01-PDMHC-Amplicon_S71_L001_R1_001.fastq.gz"
Read in 2207 paired-sequences, outputted 1364 filtered paired-sequences.
[1] 32
[1] "C-gunnisoni-ELMA02-PDMHC-Amplicon_S70_L001_R1_001.fastq.gz"
Read in 3136 paired-sequences, outputted 2202 filtered paired-sequences.
[1] 33
[1] "C-gunnisoni-ELMA04-PDMHC-Amplicon_S84_L001_R1_001.fastq.gz"
Read in 2932 paired-sequences, outputted 2007 filtered paired-sequences.
[1] 34
[1] "C-gunnisoni-ELMA05-PDMHC-Amplicon_S83_L001_R1_001.fastq.gz"
Read in 4455 paired-sequences, outputted 2882 filtered paired-sequences.
[1] 35
[1] "C-gunnisoni-ELMA06-PDMHC-Amplicon_S82_L001_R1_001.fastq.gz"
Read in 5414 paired-sequences, outputted 3136 filtered paired-sequences.
[1] 36
[1] "C-gunnisoni-ELMA07-PDMHC-Amplicon_S81_L001_R1_001.fastq.gz"
Read in 1454 paired-sequences, outputted 691 filtered paired-sequences.
[1] 37
[1] "C-gunnisoni-ELMA08-PDMHC-Amplicon_S80_L001_R1_001.fastq.gz"
Read in 5336 paired-sequences, outputted 3389 filtered paired-sequences.
[1] 38
[1] "C-gunnisoni-ELMA09-PDMHC-Amplicon_S79_L001_R1_001.fastq.gz"
Read in 1903 paired-sequences, outputted 1113 filtered paired-sequences.
[1] 39
[1] "C-gunnisoni-ELMA10-PDMHC-Amplicon_S78_L001_R1_001.fastq.gz"
Read in 2448 paired-sequences, outputted 1786 filtered paired-sequences.
[1] 40
[1] "C-gunnisoni-HOVE01-PDMHC-Amplicon_S78_L001_R1_001.fastq.gz"
Read in 4282 paired-sequences, outputted 3751 filtered paired-sequences.
[1] 41
[1] "C-gunnisoni-HOVE02-PDMHC-Amplicon_S77_L001_R1_001.fastq.gz"
Read in 5363 paired-sequences, outputted 4777 filtered paired-sequences.
[1] 42
[1] "C-gunnisoni-HOVE03-PDMHC-Amplicon_S76_L001_R1_001.fastq.gz"
Read in 5307 paired-sequences, outputted 4764 filtered paired-sequences.
[1] 43
[1] "C-gunnisoni-HOVE04-PDMHC-Amplicon_S75_L001_R1_001.fastq.gz"
Read in 2608 paired-sequences, outputted 2354 filtered paired-sequences.
[1] 44
[1] "C-gunnisoni-HOVE05-PDMHC-Amplicon_S74_L001_R1_001.fastq.gz"
Read in 5072 paired-sequences, outputted 4506 filtered paired-sequences.
[1] 45
[1] "C-gunnisoni-HOVE06-PDMHC-Amplicon_S73_L001_R1_001.fastq.gz"
Read in 2035 paired-sequences, outputted 1802 filtered paired-sequences.
[1] 46
[1] "C-gunnisoni-HOVE07-PDMHC-Amplicon_S72_L001_R1_001.fastq.gz"
Read in 3845 paired-sequences, outputted 3429 filtered paired-sequences.
[1] 47
[1] "C-gunnisoni-HOVE08-PDMHC-Amplicon_S71_L001_R1_001.fastq.gz"
Read in 3829 paired-sequences, outputted 3344 filtered paired-sequences.
[1] 48
[1] "C-gunnisoni-HOVE09-PDMHC-Amplicon_S86_L001_R1_001.fastq.gz"
Read in 7144 paired-sequences, outputted 6310 filtered paired-sequences.
[1] 49
[1] "C-gunnisoni-HOVE10-PDMHC-Amplicon_S85_L001_R1_001.fastq.gz"
Read in 4077 paired-sequences, outputted 3613 filtered paired-sequences.
[1] 50
[1] "C-gunnisoni-HOVE11-PDMHC-Amplicon_S84_L001_R1_001.fastq.gz"
Read in 2388 paired-sequences, outputted 2134 filtered paired-sequences.
[1] 51
[1] "C-gunnisoni-HOVE12-PDMHC-Amplicon_S83_L001_R1_001.fastq.gz"
Read in 4639 paired-sequences, outputted 4063 filtered paired-sequences.
[1] 52
[1] "C-gunnisoni-HOVE13-PDMHC-Amplicon_S82_L001_R1_001.fastq.gz"
Read in 5118 paired-sequences, outputted 4510 filtered paired-sequences.
[1] 53
[1] "C-gunnisoni-HOVE14-PDMHC-Amplicon_S81_L001_R1_001.fastq.gz"
Read in 3668 paired-sequences, outputted 3246 filtered paired-sequences.
[1] 54
[1] "C-gunnisoni-HOVE15-PDMHC-Amplicon_S80_L001_R1_001.fastq.gz"
Read in 4138 paired-sequences, outputted 3785 filtered paired-sequences.
[1] 55
[1] "C-gunnisoni-HOVE16-PDMHC-Amplicon_S79_L001_R1_001.fastq.gz"
Read in 2369 paired-sequences, outputted 2096 filtered paired-sequences.
[1] 56
[1] "C-gunnisoni-HOVE18-PDMHC-Amplicon_S50_L001_R1_001.fastq.gz"
Read in 2091 paired-sequences, outputted 1133 filtered paired-sequences.
[1] 57
[1] "C-gunnisoni-HOVE19-PDMHC-Amplicon_S49_L001_R1_001.fastq.gz"
Read in 4080 paired-sequences, outputted 2788 filtered paired-sequences.
[1] 58
[1] "C-gunnisoni-HOVE20-PDMHC-Amplicon_S48_L001_R1_001.fastq.gz"
Read in 1749 paired-sequences, outputted 966 filtered paired-sequences.
[1] 59
[1] "C-gunnisoni-HOVE21-PDMHC-Amplicon_S89_L001_R1_001.fastq.gz"
Read in 5552 paired-sequences, outputted 4889 filtered paired-sequences.
[1] 60
[1] "C-gunnisoni-HUTR01-PDMHC-Amplicon_S88_L001_R1_001.fastq.gz"
Read in 1901 paired-sequences, outputted 1668 filtered paired-sequences.
[1] 61
[1] "C-gunnisoni-HUTR02-PDMHC-Amplicon_S87_L001_R1_001.fastq.gz"
Read in 2321 paired-sequences, outputted 2047 filtered paired-sequences.
[1] 62
[1] "C-gunnisoni-HUTR04-PDMHC-Amplicon_S57_L001_R1_001.fastq.gz"
Read in 1730 paired-sequences, outputted 1070 filtered paired-sequences.
[1] 63
[1] "C-gunnisoni-HUTR05-PDMHC-Amplicon_S56_L001_R1_001.fastq.gz"
Read in 2647 paired-sequences, outputted 1368 filtered paired-sequences.
[1] 64
[1] "C-gunnisoni-HUTR06-PDMHC-Amplicon_S55_L001_R1_001.fastq.gz"
Read in 1405 paired-sequences, outputted 682 filtered paired-sequences.
[1] 65
[1] "C-gunnisoni-HUTR07-PDMHC-Amplicon_S54_L001_R1_001.fastq.gz"
Read in 841 paired-sequences, outputted 506 filtered paired-sequences.
[1] 66
[1] "C-gunnisoni-HUTR08-PDMHC-Amplicon_S53_L001_R1_001.fastq.gz"
Read in 777 paired-sequences, outputted 473 filtered paired-sequences.
[1] 67
[1] "C-gunnisoni-HUTR09-PDMHC-Amplicon_S52_L001_R1_001.fastq.gz"
Read in 1961 paired-sequences, outputted 1279 filtered paired-sequences.
[1] 68
[1] "C-gunnisoni-HUTR11-PDMHC-Amplicon_S51_L001_R1_001.fastq.gz"
Read in 2402 paired-sequences, outputted 1481 filtered paired-sequences.
[1] 69
[1] "C-gunnisoni-HUTR12-PDMHC-Amplicon_S62_L001_R1_001.fastq.gz"
Read in 1676 paired-sequences, outputted 1025 filtered paired-sequences.
[1] 70
[1] "C-gunnisoni-HUTR14-PDMHC-Amplicon_S61_L001_R1_001.fastq.gz"
Read in 3129 paired-sequences, outputted 2005 filtered paired-sequences.
[1] 71
[1] "C-gunnisoni-HUTR15-PDMHC-Amplicon_S60_L001_R1_001.fastq.gz"
Read in 4418 paired-sequences, outputted 2722 filtered paired-sequences.
[1] 72
[1] "C-gunnisoni-HUTR16-PDMHC-Amplicon_S59_L001_R1_001.fastq.gz"
Read in 3630 paired-sequences, outputted 2398 filtered paired-sequences.
[1] 73
[1] "C-gunnisoni-HUTR17-PDMHC-Amplicon_S58_L001_R1_001.fastq.gz"
Read in 3196 paired-sequences, outputted 1896 filtered paired-sequences.
[1] 74
[1] "C-gunnisoni-HUTR20-PDMHC-Amplicon_S69_L001_R1_001.fastq.gz"
Read in 2255 paired-sequences, outputted 1413 filtered paired-sequences.
[1] 75
[1] "C-gunnisoni-HUTR21-PDMHC-Amplicon_S68_L001_R1_001.fastq.gz"
Read in 4599 paired-sequences, outputted 3275 filtered paired-sequences.
[1] 76
[1] "C-gunnisoni-HUTR22-PDMHC-Amplicon_S67_L001_R1_001.fastq.gz"
Read in 2437 paired-sequences, outputted 1439 filtered paired-sequences.
[1] 77
[1] "C-gunnisoni-HUTR24-PDMHC-Amplicon_S66_L001_R1_001.fastq.gz"
Read in 1646 paired-sequences, outputted 1024 filtered paired-sequences.
[1] 78
[1] "C-gunnisoni-RSF01-PDMHC-Amplicon_S65_L001_R1_001.fastq.gz"
Read in 4445 paired-sequences, outputted 2705 filtered paired-sequences.
[1] 79
[1] "C-gunnisoni-RSF02-PDMHC-Amplicon_S64_L001_R1_001.fastq.gz"
Read in 1825 paired-sequences, outputted 949 filtered paired-sequences.
[1] 80
[1] "C-gunnisoni-RSF03-PDMHC-Amplicon_S63_L001_R1_001.fastq.gz"
Read in 3455 paired-sequences, outputted 2125 filtered paired-sequences.
[1] 81
[1] "Cynomys-gunnisoni-105-55_S8_L001_R1_001.fastq.gz"
Read in 24892 paired-sequences, outputted 15327 filtered paired-sequences.
[1] 82
[1] "Cynomys-gunnisoni-61-30_S5_L001_R1_001.fastq.gz"
Read in 12756 paired-sequences, outputted 8343 filtered paired-sequences.
[1] 83
[1] "Cynomys-gunnisoni-74-35_S6_L001_R1_001.fastq.gz"
Read in 9990 paired-sequences, outputted 6575 filtered paired-sequences.
[1] 84
[1] "Cynomys-gunnisoni-77-1_S4_L001_R1_001.fastq.gz"
Read in 11668 paired-sequences, outputted 7871 filtered paired-sequences.
[1] 85
[1] "Cynomys-gunnisoni-80-34_S7_L001_R1_001.fastq.gz"
Read in 17542 paired-sequences, outputted 11515 filtered paired-sequences.
```

```python
>>> derepFs <- lapply(filtFs, derepFastq, verbose=TRUE)
>>> derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS01-PDMHC-Amplicon_S77_L001_R1_001_filt.fastq.gz
Encountered 727 unique sequences from 1273 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS02-PDMHC-Amplicon_S76_L001_R1_001_filt.fastq.gz
Encountered 502 unique sequences from 947 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS03-PDMHC-Amplicon_S75_L001_R1_001_filt.fastq.gz
Encountered 1079 unique sequences from 2039 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS04-PDMHC-Amplicon_S74_L001_R1_001_filt.fastq.gz
Encountered 1117 unique sequences from 2051 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS05-PDMHC-Amplicon_S73_L001_R1_001_filt.fastq.gz
Encountered 1743 unique sequences from 3347 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS06-PDMHC-Amplicon_S72_L001_R1_001_filt.fastq.gz
Encountered 370 unique sequences from 585 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL01-PDMHC-Amplicon_S54_L001_R1_001_filt.fastq.gz
Encountered 1413 unique sequences from 5292 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL02-PDMHC-Amplicon_S53_L001_R1_001_filt.fastq.gz
Encountered 474 unique sequences from 1200 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL03-PDMHC-Amplicon_S52_L001_R1_001_filt.fastq.gz
Encountered 745 unique sequences from 2020 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL04-PDMHC-Amplicon_S51_L001_R1_001_filt.fastq.gz
Encountered 592 unique sequences from 1699 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL05-PDMHC-Amplicon_S50_L001_R1_001_filt.fastq.gz
Encountered 728 unique sequences from 2413 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL06-PDMHC-Amplicon_S49_L001_R1_001_filt.fastq.gz
Encountered 708 unique sequences from 2140 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL07-PDMHC-Amplicon_S48_L001_R1_001_filt.fastq.gz
Encountered 704 unique sequences from 2291 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL08-PDMHC-Amplicon_S47_L001_R1_001_filt.fastq.gz
Encountered 826 unique sequences from 2651 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL09-PDMHC-Amplicon_S62_L001_R1_001_filt.fastq.gz
Encountered 1016 unique sequences from 3407 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL10-PDMHC-Amplicon_S61_L001_R1_001_filt.fastq.gz
Encountered 879 unique sequences from 2634 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL11-PDMHC-Amplicon_S60_L001_R1_001_filt.fastq.gz
Encountered 996 unique sequences from 3159 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL12-PDMHC-Amplicon_S59_L001_R1_001_filt.fastq.gz
Encountered 1065 unique sequences from 3530 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL13-PDMHC-Amplicon_S58_L001_R1_001_filt.fastq.gz
Encountered 1353 unique sequences from 5693 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL14-PDMHC-Amplicon_S57_L001_R1_001_filt.fastq.gz
Encountered 512 unique sequences from 1233 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL15-PDMHC-Amplicon_S56_L001_R1_001_filt.fastq.gz
Encountered 1063 unique sequences from 3728 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL16-PDMHC-Amplicon_S55_L001_R1_001_filt.fastq.gz
Encountered 1690 unique sequences from 5412 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL17-PDMHC-Amplicon_S70_L001_R1_001_filt.fastq.gz
Encountered 1261 unique sequences from 3671 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL18-PDMHC-Amplicon_S69_L001_R1_001_filt.fastq.gz
Encountered 339 unique sequences from 836 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL19-PDMHC-Amplicon_S68_L001_R1_001_filt.fastq.gz
Encountered 940 unique sequences from 2997 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL20-PDMHC-Amplicon_S67_L001_R1_001_filt.fastq.gz
Encountered 1141 unique sequences from 4283 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL21-PDMHC-Amplicon_S66_L001_R1_001_filt.fastq.gz
Encountered 437 unique sequences from 1238 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL22-PDMHC-Amplicon_S65_L001_R1_001_filt.fastq.gz
Encountered 818 unique sequences from 2821 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL23-PDMHC-Amplicon_S64_L001_R1_001_filt.fastq.gz
Encountered 760 unique sequences from 2495 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL24-PDMHC-Amplicon_S63_L001_R1_001_filt.fastq.gz
Encountered 601 unique sequences from 1808 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA01-PDMHC-Amplicon_S71_L001_R1_001_filt.fastq.gz
Encountered 717 unique sequences from 1364 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA02-PDMHC-Amplicon_S70_L001_R1_001_filt.fastq.gz
Encountered 1026 unique sequences from 2202 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA04-PDMHC-Amplicon_S84_L001_R1_001_filt.fastq.gz
Encountered 1027 unique sequences from 2007 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA05-PDMHC-Amplicon_S83_L001_R1_001_filt.fastq.gz
Encountered 1366 unique sequences from 2882 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA06-PDMHC-Amplicon_S82_L001_R1_001_filt.fastq.gz
Encountered 1576 unique sequences from 3136 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA07-PDMHC-Amplicon_S81_L001_R1_001_filt.fastq.gz
Encountered 446 unique sequences from 691 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA08-PDMHC-Amplicon_S80_L001_R1_001_filt.fastq.gz
Encountered 1531 unique sequences from 3389 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA09-PDMHC-Amplicon_S79_L001_R1_001_filt.fastq.gz
Encountered 671 unique sequences from 1113 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA10-PDMHC-Amplicon_S78_L001_R1_001_filt.fastq.gz
Encountered 817 unique sequences from 1786 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE01-PDMHC-Amplicon_S78_L001_R1_001_filt.fastq.gz
Encountered 1225 unique sequences from 3751 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE02-PDMHC-Amplicon_S77_L001_R1_001_filt.fastq.gz
Encountered 1460 unique sequences from 4777 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE03-PDMHC-Amplicon_S76_L001_R1_001_filt.fastq.gz
Encountered 1335 unique sequences from 4764 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE04-PDMHC-Amplicon_S75_L001_R1_001_filt.fastq.gz
Encountered 791 unique sequences from 2354 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE05-PDMHC-Amplicon_S74_L001_R1_001_filt.fastq.gz
Encountered 1337 unique sequences from 4506 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE06-PDMHC-Amplicon_S73_L001_R1_001_filt.fastq.gz
Encountered 622 unique sequences from 1802 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE07-PDMHC-Amplicon_S72_L001_R1_001_filt.fastq.gz
Encountered 1063 unique sequences from 3429 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE08-PDMHC-Amplicon_S71_L001_R1_001_filt.fastq.gz
Encountered 1251 unique sequences from 3344 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE09-PDMHC-Amplicon_S86_L001_R1_001_filt.fastq.gz
Encountered 1688 unique sequences from 6310 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE10-PDMHC-Amplicon_S85_L001_R1_001_filt.fastq.gz
Encountered 1208 unique sequences from 3613 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE11-PDMHC-Amplicon_S84_L001_R1_001_filt.fastq.gz
Encountered 772 unique sequences from 2134 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE12-PDMHC-Amplicon_S83_L001_R1_001_filt.fastq.gz
Encountered 1280 unique sequences from 4063 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE13-PDMHC-Amplicon_S82_L001_R1_001_filt.fastq.gz
Encountered 1420 unique sequences from 4510 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE14-PDMHC-Amplicon_S81_L001_R1_001_filt.fastq.gz
Encountered 1054 unique sequences from 3246 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE15-PDMHC-Amplicon_S80_L001_R1_001_filt.fastq.gz
Encountered 908 unique sequences from 3785 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE16-PDMHC-Amplicon_S79_L001_R1_001_filt.fastq.gz
Encountered 795 unique sequences from 2096 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE18-PDMHC-Amplicon_S50_L001_R1_001_filt.fastq.gz
Encountered 627 unique sequences from 1133 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE19-PDMHC-Amplicon_S49_L001_R1_001_filt.fastq.gz
Encountered 1353 unique sequences from 2788 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE20-PDMHC-Amplicon_S48_L001_R1_001_filt.fastq.gz
Encountered 559 unique sequences from 966 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE21-PDMHC-Amplicon_S89_L001_R1_001_filt.fastq.gz
Encountered 1335 unique sequences from 4889 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR01-PDMHC-Amplicon_S88_L001_R1_001_filt.fastq.gz
Encountered 608 unique sequences from 1668 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR02-PDMHC-Amplicon_S87_L001_R1_001_filt.fastq.gz
Encountered 778 unique sequences from 2047 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR04-PDMHC-Amplicon_S57_L001_R1_001_filt.fastq.gz
Encountered 609 unique sequences from 1070 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR05-PDMHC-Amplicon_S56_L001_R1_001_filt.fastq.gz
Encountered 837 unique sequences from 1368 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR06-PDMHC-Amplicon_S55_L001_R1_001_filt.fastq.gz
Encountered 462 unique sequences from 682 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR07-PDMHC-Amplicon_S54_L001_R1_001_filt.fastq.gz
Encountered 289 unique sequences from 506 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR08-PDMHC-Amplicon_S53_L001_R1_001_filt.fastq.gz
Encountered 304 unique sequences from 473 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR09-PDMHC-Amplicon_S52_L001_R1_001_filt.fastq.gz
Encountered 704 unique sequences from 1279 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR11-PDMHC-Amplicon_S51_L001_R1_001_filt.fastq.gz
Encountered 858 unique sequences from 1481 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR12-PDMHC-Amplicon_S62_L001_R1_001_filt.fastq.gz
Encountered 584 unique sequences from 1025 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR14-PDMHC-Amplicon_S61_L001_R1_001_filt.fastq.gz
Encountered 975 unique sequences from 2005 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR15-PDMHC-Amplicon_S60_L001_R1_001_filt.fastq.gz
Encountered 1283 unique sequences from 2722 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR16-PDMHC-Amplicon_S59_L001_R1_001_filt.fastq.gz
Encountered 1120 unique sequences from 2398 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR17-PDMHC-Amplicon_S58_L001_R1_001_filt.fastq.gz
Encountered 946 unique sequences from 1896 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR20-PDMHC-Amplicon_S69_L001_R1_001_filt.fastq.gz
Encountered 772 unique sequences from 1413 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR21-PDMHC-Amplicon_S68_L001_R1_001_filt.fastq.gz
Encountered 1305 unique sequences from 3275 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR22-PDMHC-Amplicon_S67_L001_R1_001_filt.fastq.gz
Encountered 816 unique sequences from 1439 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR24-PDMHC-Amplicon_S66_L001_R1_001_filt.fastq.gz
Encountered 565 unique sequences from 1024 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-RSF01-PDMHC-Amplicon_S65_L001_R1_001_filt.fastq.gz
Encountered 1396 unique sequences from 2705 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-RSF02-PDMHC-Amplicon_S64_L001_R1_001_filt.fastq.gz
Encountered 597 unique sequences from 949 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-RSF03-PDMHC-Amplicon_S63_L001_R1_001_filt.fastq.gz
Encountered 1124 unique sequences from 2125 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-105-55_S8_L001_R1_001_filt.fastq.gz
Encountered 6637 unique sequences from 15327 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-61-30_S5_L001_R1_001_filt.fastq.gz
Encountered 3301 unique sequences from 8343 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-74-35_S6_L001_R1_001_filt.fastq.gz
Encountered 2927 unique sequences from 6575 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-77-1_S4_L001_R1_001_filt.fastq.gz
Encountered 3459 unique sequences from 7871 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-80-34_S7_L001_R1_001_filt.fastq.gz
Encountered 4788 unique sequences from 11515 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS01-PDMHC-Amplicon_S77_L001_R2_001_filt.fastq.gz
Encountered 303 unique sequences from 1273 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS02-PDMHC-Amplicon_S76_L001_R2_001_filt.fastq.gz
Encountered 224 unique sequences from 947 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS03-PDMHC-Amplicon_S75_L001_R2_001_filt.fastq.gz
Encountered 492 unique sequences from 2039 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS04-PDMHC-Amplicon_S74_L001_R2_001_filt.fastq.gz
Encountered 483 unique sequences from 2051 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS05-PDMHC-Amplicon_S73_L001_R2_001_filt.fastq.gz
Encountered 701 unique sequences from 3347 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-BLS06-PDMHC-Amplicon_S72_L001_R2_001_filt.fastq.gz
Encountered 189 unique sequences from 585 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL01-PDMHC-Amplicon_S54_L001_R2_001_filt.fastq.gz
Encountered 686 unique sequences from 5292 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL02-PDMHC-Amplicon_S53_L001_R2_001_filt.fastq.gz
Encountered 297 unique sequences from 1200 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL03-PDMHC-Amplicon_S52_L001_R2_001_filt.fastq.gz
Encountered 374 unique sequences from 2020 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL04-PDMHC-Amplicon_S51_L001_R2_001_filt.fastq.gz
Encountered 334 unique sequences from 1699 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL05-PDMHC-Amplicon_S50_L001_R2_001_filt.fastq.gz
Encountered 450 unique sequences from 2413 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL06-PDMHC-Amplicon_S49_L001_R2_001_filt.fastq.gz
Encountered 405 unique sequences from 2140 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL07-PDMHC-Amplicon_S48_L001_R2_001_filt.fastq.gz
Encountered 356 unique sequences from 2291 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL08-PDMHC-Amplicon_S47_L001_R2_001_filt.fastq.gz
Encountered 376 unique sequences from 2651 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL09-PDMHC-Amplicon_S62_L001_R2_001_filt.fastq.gz
Encountered 448 unique sequences from 3407 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL10-PDMHC-Amplicon_S61_L001_R2_001_filt.fastq.gz
Encountered 495 unique sequences from 2634 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL11-PDMHC-Amplicon_S60_L001_R2_001_filt.fastq.gz
Encountered 408 unique sequences from 3159 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL12-PDMHC-Amplicon_S59_L001_R2_001_filt.fastq.gz
Encountered 564 unique sequences from 3530 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL13-PDMHC-Amplicon_S58_L001_R2_001_filt.fastq.gz
Encountered 740 unique sequences from 5693 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL14-PDMHC-Amplicon_S57_L001_R2_001_filt.fastq.gz
Encountered 235 unique sequences from 1233 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL15-PDMHC-Amplicon_S56_L001_R2_001_filt.fastq.gz
Encountered 470 unique sequences from 3728 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL16-PDMHC-Amplicon_S55_L001_R2_001_filt.fastq.gz
Encountered 783 unique sequences from 5412 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL17-PDMHC-Amplicon_S70_L001_R2_001_filt.fastq.gz
Encountered 616 unique sequences from 3671 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL18-PDMHC-Amplicon_S69_L001_R2_001_filt.fastq.gz
Encountered 160 unique sequences from 836 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL19-PDMHC-Amplicon_S68_L001_R2_001_filt.fastq.gz
Encountered 479 unique sequences from 2997 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL20-PDMHC-Amplicon_S67_L001_R2_001_filt.fastq.gz
Encountered 604 unique sequences from 4283 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL21-PDMHC-Amplicon_S66_L001_R2_001_filt.fastq.gz
Encountered 245 unique sequences from 1238 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL22-PDMHC-Amplicon_S65_L001_R2_001_filt.fastq.gz
Encountered 429 unique sequences from 2821 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL23-PDMHC-Amplicon_S64_L001_R2_001_filt.fastq.gz
Encountered 429 unique sequences from 2495 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-CRL24-PDMHC-Amplicon_S63_L001_R2_001_filt.fastq.gz
Encountered 264 unique sequences from 1808 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA01-PDMHC-Amplicon_S71_L001_R2_001_filt.fastq.gz
Encountered 280 unique sequences from 1364 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA02-PDMHC-Amplicon_S70_L001_R2_001_filt.fastq.gz
Encountered 495 unique sequences from 2202 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA04-PDMHC-Amplicon_S84_L001_R2_001_filt.fastq.gz
Encountered 463 unique sequences from 2007 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA05-PDMHC-Amplicon_S83_L001_R2_001_filt.fastq.gz
Encountered 538 unique sequences from 2882 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA06-PDMHC-Amplicon_S82_L001_R2_001_filt.fastq.gz
Encountered 614 unique sequences from 3136 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA07-PDMHC-Amplicon_S81_L001_R2_001_filt.fastq.gz
Encountered 179 unique sequences from 691 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA08-PDMHC-Amplicon_S80_L001_R2_001_filt.fastq.gz
Encountered 617 unique sequences from 3389 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA09-PDMHC-Amplicon_S79_L001_R2_001_filt.fastq.gz
Encountered 261 unique sequences from 1113 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-ELMA10-PDMHC-Amplicon_S78_L001_R2_001_filt.fastq.gz
Encountered 373 unique sequences from 1786 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE01-PDMHC-Amplicon_S78_L001_R2_001_filt.fastq.gz
Encountered 589 unique sequences from 3751 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE02-PDMHC-Amplicon_S77_L001_R2_001_filt.fastq.gz
Encountered 680 unique sequences from 4777 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE03-PDMHC-Amplicon_S76_L001_R2_001_filt.fastq.gz
Encountered 715 unique sequences from 4764 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE04-PDMHC-Amplicon_S75_L001_R2_001_filt.fastq.gz
Encountered 434 unique sequences from 2354 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE05-PDMHC-Amplicon_S74_L001_R2_001_filt.fastq.gz
Encountered 606 unique sequences from 4506 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE06-PDMHC-Amplicon_S73_L001_R2_001_filt.fastq.gz
Encountered 290 unique sequences from 1802 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE07-PDMHC-Amplicon_S72_L001_R2_001_filt.fastq.gz
Encountered 539 unique sequences from 3429 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE08-PDMHC-Amplicon_S71_L001_R2_001_filt.fastq.gz
Encountered 647 unique sequences from 3344 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE09-PDMHC-Amplicon_S86_L001_R2_001_filt.fastq.gz
Encountered 809 unique sequences from 6310 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE10-PDMHC-Amplicon_S85_L001_R2_001_filt.fastq.gz
Encountered 605 unique sequences from 3613 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE11-PDMHC-Amplicon_S84_L001_R2_001_filt.fastq.gz
Encountered 381 unique sequences from 2134 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE12-PDMHC-Amplicon_S83_L001_R2_001_filt.fastq.gz
Encountered 648 unique sequences from 4063 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE13-PDMHC-Amplicon_S82_L001_R2_001_filt.fastq.gz
Encountered 598 unique sequences from 4510 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE14-PDMHC-Amplicon_S81_L001_R2_001_filt.fastq.gz
Encountered 496 unique sequences from 3246 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE15-PDMHC-Amplicon_S80_L001_R2_001_filt.fastq.gz
Encountered 448 unique sequences from 3785 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE16-PDMHC-Amplicon_S79_L001_R2_001_filt.fastq.gz
Encountered 395 unique sequences from 2096 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE18-PDMHC-Amplicon_S50_L001_R2_001_filt.fastq.gz
Encountered 250 unique sequences from 1133 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE19-PDMHC-Amplicon_S49_L001_R2_001_filt.fastq.gz
Encountered 553 unique sequences from 2788 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE20-PDMHC-Amplicon_S48_L001_R2_001_filt.fastq.gz
Encountered 254 unique sequences from 966 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HOVE21-PDMHC-Amplicon_S89_L001_R2_001_filt.fastq.gz
Encountered 613 unique sequences from 4889 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR01-PDMHC-Amplicon_S88_L001_R2_001_filt.fastq.gz
Encountered 315 unique sequences from 1668 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR02-PDMHC-Amplicon_S87_L001_R2_001_filt.fastq.gz
Encountered 370 unique sequences from 2047 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR04-PDMHC-Amplicon_S57_L001_R2_001_filt.fastq.gz
Encountered 263 unique sequences from 1070 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR05-PDMHC-Amplicon_S56_L001_R2_001_filt.fastq.gz
Encountered 329 unique sequences from 1368 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR06-PDMHC-Amplicon_S55_L001_R2_001_filt.fastq.gz
Encountered 185 unique sequences from 682 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR07-PDMHC-Amplicon_S54_L001_R2_001_filt.fastq.gz
Encountered 124 unique sequences from 506 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR08-PDMHC-Amplicon_S53_L001_R2_001_filt.fastq.gz
Encountered 150 unique sequences from 473 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR09-PDMHC-Amplicon_S52_L001_R2_001_filt.fastq.gz
Encountered 343 unique sequences from 1279 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR11-PDMHC-Amplicon_S51_L001_R2_001_filt.fastq.gz
Encountered 369 unique sequences from 1481 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR12-PDMHC-Amplicon_S62_L001_R2_001_filt.fastq.gz
Encountered 246 unique sequences from 1025 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR14-PDMHC-Amplicon_S61_L001_R2_001_filt.fastq.gz
Encountered 390 unique sequences from 2005 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR15-PDMHC-Amplicon_S60_L001_R2_001_filt.fastq.gz
Encountered 553 unique sequences from 2722 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR16-PDMHC-Amplicon_S59_L001_R2_001_filt.fastq.gz
Encountered 400 unique sequences from 2398 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR17-PDMHC-Amplicon_S58_L001_R2_001_filt.fastq.gz
Encountered 375 unique sequences from 1896 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR20-PDMHC-Amplicon_S69_L001_R2_001_filt.fastq.gz
Encountered 369 unique sequences from 1413 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR21-PDMHC-Amplicon_S68_L001_R2_001_filt.fastq.gz
Encountered 556 unique sequences from 3275 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR22-PDMHC-Amplicon_S67_L001_R2_001_filt.fastq.gz
Encountered 346 unique sequences from 1439 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-HUTR24-PDMHC-Amplicon_S66_L001_R2_001_filt.fastq.gz
Encountered 281 unique sequences from 1024 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-RSF01-PDMHC-Amplicon_S65_L001_R2_001_filt.fastq.gz
Encountered 589 unique sequences from 2705 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-RSF02-PDMHC-Amplicon_S64_L001_R2_001_filt.fastq.gz
Encountered 250 unique sequences from 949 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/C-gunnisoni-RSF03-PDMHC-Amplicon_S63_L001_R2_001_filt.fastq.gz
Encountered 480 unique sequences from 2125 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-105-55_S8_L001_R2_001_filt.fastq.gz
Encountered 2165 unique sequences from 15327 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-61-30_S5_L001_R2_001_filt.fastq.gz
Encountered 1408 unique sequences from 8343 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-74-35_S6_L001_R2_001_filt.fastq.gz
Encountered 1207 unique sequences from 6575 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-77-1_S4_L001_R2_001_filt.fastq.gz
Encountered 1655 unique sequences from 7871 total sequences read.
Dereplicating sequence entries in Fastq file: AmpliconSeq/Cynomys-gunnisoni-80-34_S7_L001_R2_001_filt.fastq.gz
Encountered 2122 unique sequences from 11515 total sequences read.
```

```python
>>> # Name the derep-class objects by the sample names
... sam_names <- sapply(strsplit(fnFs, "/"), tail, n=1)
>>> sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)
>>> names(derepFs) <- sam_names
>>> names(derepRs) <- sam_names
```

```python
>>> dadaFs <- dada(derepFs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
Sample 1 - 1273 reads in 727 unique sequences.
Sample 2 - 947 reads in 502 unique sequences.
Sample 3 - 2039 reads in 1079 unique sequences.
Sample 4 - 2051 reads in 1117 unique sequences.
Sample 5 - 3347 reads in 1743 unique sequences.
Sample 6 - 585 reads in 370 unique sequences.
Sample 7 - 5292 reads in 1413 unique sequences.
Sample 8 - 1200 reads in 474 unique sequences.
Sample 9 - 2020 reads in 745 unique sequences.
Sample 10 - 1699 reads in 592 unique sequences.
Sample 11 - 2413 reads in 728 unique sequences.
Sample 12 - 2140 reads in 708 unique sequences.
Sample 13 - 2291 reads in 704 unique sequences.
Sample 14 - 2651 reads in 826 unique sequences.
Sample 15 - 3407 reads in 1016 unique sequences.
Sample 16 - 2634 reads in 879 unique sequences.
Sample 17 - 3159 reads in 996 unique sequences.
Sample 18 - 3530 reads in 1065 unique sequences.
Sample 19 - 5693 reads in 1353 unique sequences.
Sample 20 - 1233 reads in 512 unique sequences.
Sample 21 - 3728 reads in 1063 unique sequences.
Sample 22 - 5412 reads in 1690 unique sequences.
Sample 23 - 3671 reads in 1261 unique sequences.
Sample 24 - 836 reads in 339 unique sequences.
Sample 25 - 2997 reads in 940 unique sequences.
Sample 26 - 4283 reads in 1141 unique sequences.
Sample 27 - 1238 reads in 437 unique sequences.
Sample 28 - 2821 reads in 818 unique sequences.
Sample 29 - 2495 reads in 760 unique sequences.
Sample 30 - 1808 reads in 601 unique sequences.
Sample 31 - 1364 reads in 717 unique sequences.
Sample 32 - 2202 reads in 1026 unique sequences.
Sample 33 - 2007 reads in 1027 unique sequences.
Sample 34 - 2882 reads in 1366 unique sequences.
Sample 35 - 3136 reads in 1576 unique sequences.
Sample 36 - 691 reads in 446 unique sequences.
Sample 37 - 3389 reads in 1531 unique sequences.
Sample 38 - 1113 reads in 671 unique sequences.
Sample 39 - 1786 reads in 817 unique sequences.
Sample 40 - 3751 reads in 1225 unique sequences.
Sample 41 - 4777 reads in 1460 unique sequences.
Sample 42 - 4764 reads in 1335 unique sequences.
Sample 43 - 2354 reads in 791 unique sequences.
Sample 44 - 4506 reads in 1337 unique sequences.
Sample 45 - 1802 reads in 622 unique sequences.
Sample 46 - 3429 reads in 1063 unique sequences.
Sample 47 - 3344 reads in 1251 unique sequences.
Sample 48 - 6310 reads in 1688 unique sequences.
Sample 49 - 3613 reads in 1208 unique sequences.
Sample 50 - 2134 reads in 772 unique sequences.
Sample 51 - 4063 reads in 1280 unique sequences.
Sample 52 - 4510 reads in 1420 unique sequences.
Sample 53 - 3246 reads in 1054 unique sequences.
Sample 54 - 3785 reads in 908 unique sequences.
Sample 55 - 2096 reads in 795 unique sequences.
Sample 56 - 1133 reads in 627 unique sequences.
Sample 57 - 2788 reads in 1353 unique sequences.
Sample 58 - 966 reads in 559 unique sequences.
Sample 59 - 4889 reads in 1335 unique sequences.
Sample 60 - 1668 reads in 608 unique sequences.
Sample 61 - 2047 reads in 778 unique sequences.
Sample 62 - 1070 reads in 609 unique sequences.
Sample 63 - 1368 reads in 837 unique sequences.
Sample 64 - 682 reads in 462 unique sequences.
Sample 65 - 506 reads in 289 unique sequences.
Sample 66 - 473 reads in 304 unique sequences.
Sample 67 - 1279 reads in 704 unique sequences.
Sample 68 - 1481 reads in 858 unique sequences.
Sample 69 - 1025 reads in 584 unique sequences.
Sample 70 - 2005 reads in 975 unique sequences.
Sample 71 - 2722 reads in 1283 unique sequences.
Sample 72 - 2398 reads in 1120 unique sequences.
Sample 73 - 1896 reads in 946 unique sequences.
Sample 74 - 1413 reads in 772 unique sequences.
Sample 75 - 3275 reads in 1305 unique sequences.
Sample 76 - 1439 reads in 816 unique sequences.
Sample 77 - 1024 reads in 565 unique sequences.
Sample 78 - 2705 reads in 1396 unique sequences.
Sample 79 - 949 reads in 597 unique sequences.
Sample 80 - 2125 reads in 1124 unique sequences.
Sample 81 - 15327 reads in 6637 unique sequences.
Sample 82 - 8343 reads in 3301 unique sequences.
Sample 83 - 6575 reads in 2927 unique sequences.
Sample 84 - 7871 reads in 3459 unique sequences.
Sample 85 - 11515 reads in 4788 unique sequences.
   Consist step 2 
   Consist step 3 
   Consist step 4 


Convergence after  4  rounds.
```

```python
>>> dadaRs <- dada(derepRs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
Sample 1 - 1273 reads in 303 unique sequences.
Sample 2 - 947 reads in 224 unique sequences.
Sample 3 - 2039 reads in 492 unique sequences.
Sample 4 - 2051 reads in 483 unique sequences.
Sample 5 - 3347 reads in 701 unique sequences.
Sample 6 - 585 reads in 189 unique sequences.
Sample 7 - 5292 reads in 686 unique sequences.
Sample 8 - 1200 reads in 297 unique sequences.
Sample 9 - 2020 reads in 374 unique sequences.
Sample 10 - 1699 reads in 334 unique sequences.
Sample 11 - 2413 reads in 450 unique sequences.
Sample 12 - 2140 reads in 405 unique sequences.
Sample 13 - 2291 reads in 356 unique sequences.
Sample 14 - 2651 reads in 376 unique sequences.
Sample 15 - 3407 reads in 448 unique sequences.
Sample 16 - 2634 reads in 495 unique sequences.
Sample 17 - 3159 reads in 408 unique sequences.
Sample 18 - 3530 reads in 564 unique sequences.
Sample 19 - 5693 reads in 740 unique sequences.
Sample 20 - 1233 reads in 235 unique sequences.
Sample 21 - 3728 reads in 470 unique sequences.
Sample 22 - 5412 reads in 783 unique sequences.
Sample 23 - 3671 reads in 616 unique sequences.
Sample 24 - 836 reads in 160 unique sequences.
Sample 25 - 2997 reads in 479 unique sequences.
Sample 26 - 4283 reads in 604 unique sequences.
Sample 27 - 1238 reads in 245 unique sequences.
Sample 28 - 2821 reads in 429 unique sequences.
Sample 29 - 2495 reads in 429 unique sequences.
Sample 30 - 1808 reads in 264 unique sequences.
Sample 31 - 1364 reads in 280 unique sequences.
Sample 32 - 2202 reads in 495 unique sequences.
Sample 33 - 2007 reads in 463 unique sequences.
Sample 34 - 2882 reads in 538 unique sequences.
Sample 35 - 3136 reads in 614 unique sequences.
Sample 36 - 691 reads in 179 unique sequences.
Sample 37 - 3389 reads in 617 unique sequences.
Sample 38 - 1113 reads in 261 unique sequences.
Sample 39 - 1786 reads in 373 unique sequences.
Sample 40 - 3751 reads in 589 unique sequences.
Sample 41 - 4777 reads in 680 unique sequences.
Sample 42 - 4764 reads in 715 unique sequences.
Sample 43 - 2354 reads in 434 unique sequences.
Sample 44 - 4506 reads in 606 unique sequences.
Sample 45 - 1802 reads in 290 unique sequences.
Sample 46 - 3429 reads in 539 unique sequences.
Sample 47 - 3344 reads in 647 unique sequences.
Sample 48 - 6310 reads in 809 unique sequences.
Sample 49 - 3613 reads in 605 unique sequences.
Sample 50 - 2134 reads in 381 unique sequences.
Sample 51 - 4063 reads in 648 unique sequences.
Sample 52 - 4510 reads in 598 unique sequences.
Sample 53 - 3246 reads in 496 unique sequences.
Sample 54 - 3785 reads in 448 unique sequences.
Sample 55 - 2096 reads in 395 unique sequences.
Sample 56 - 1133 reads in 250 unique sequences.
Sample 57 - 2788 reads in 553 unique sequences.
Sample 58 - 966 reads in 254 unique sequences.
Sample 59 - 4889 reads in 613 unique sequences.
Sample 60 - 1668 reads in 315 unique sequences.
Sample 61 - 2047 reads in 370 unique sequences.
Sample 62 - 1070 reads in 263 unique sequences.
Sample 63 - 1368 reads in 329 unique sequences.
Sample 64 - 682 reads in 185 unique sequences.
Sample 65 - 506 reads in 124 unique sequences.
Sample 66 - 473 reads in 150 unique sequences.
Sample 67 - 1279 reads in 343 unique sequences.
Sample 68 - 1481 reads in 369 unique sequences.
Sample 69 - 1025 reads in 246 unique sequences.
Sample 70 - 2005 reads in 390 unique sequences.
Sample 71 - 2722 reads in 553 unique sequences.
Sample 72 - 2398 reads in 400 unique sequences.
Sample 73 - 1896 reads in 375 unique sequences.
Sample 74 - 1413 reads in 369 unique sequences.
Sample 75 - 3275 reads in 556 unique sequences.
Sample 76 - 1439 reads in 346 unique sequences.
Sample 77 - 1024 reads in 281 unique sequences.
Sample 78 - 2705 reads in 589 unique sequences.
Sample 79 - 949 reads in 250 unique sequences.
Sample 80 - 2125 reads in 480 unique sequences.
Sample 81 - 15327 reads in 2165 unique sequences.
Sample 82 - 8343 reads in 1408 unique sequences.
Sample 83 - 6575 reads in 1207 unique sequences.
Sample 84 - 7871 reads in 1655 unique sequences.
Sample 85 - 11515 reads in 2122 unique sequences.
   Consist step 2 
   Consist step 3 
   Consist step 4 


Convergence after  4  rounds.
```

```python
>>> bimFs <- sapply(dadaFs, isBimeraDenovo, verbose=TRUE)
>>> bimRs <- sapply(dadaRs, isBimeraDenovo, verbose=TRUE)
>>> print(unname(sapply(bimFs, mean)), digits=2)
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 4 bimeras out of 6 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 2 bimeras out of 5 input sequences.
Identified 0 bimeras out of 5 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 2 bimeras out of 4 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 1 bimeras out of 3 input sequences.
Identified 1 bimeras out of 4 input sequences.
Identified 2 bimeras out of 4 input sequences.
Identified 4 bimeras out of 7 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 1 bimeras out of 3 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 10 bimeras out of 12 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 5 input sequences.
Identified 3 bimeras out of 6 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 1 bimeras out of 4 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 3 bimeras out of 7 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 4 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 1 bimeras out of 5 input sequences.
Identified 0 bimeras out of 6 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 2 bimeras out of 5 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 1 bimeras out of 3 input sequences.
Identified 0 bimeras out of 5 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 3 bimeras out of 6 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 3 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 1 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 0 bimeras out of 2 input sequences.
Identified 8 bimeras out of 11 input sequences.
Identified 0 bimeras out of 5 input sequences.
Identified 0 bimeras out of 9 input sequences.
Identified 2 bimeras out of 5 input sequences.
Identified 0 bimeras out of 3 input sequences.
 [1] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
[16] 0.00 0.00 0.00 0.00 0.00 0.00 0.67 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
[31] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.40 0.00 0.00 0.50 0.00
[46] 0.00 0.00 0.00 0.33 0.25 0.50 0.57 0.00 0.00 0.00 0.00 0.33 0.00 0.00 0.00
[61] 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00
[76] 0.00 0.00 0.00 0.00 0.00 0.83 0.00 0.00 0.50 0.00
```

```python
>>> mergers <- mapply(mergePairs, dadaFs, derepFs, dadaRs, derepRs, SIMPLIFY=FALSE)
```

```python
>>> mergers.nochim <- mapply(function(mm, bF, bR) mm[!bF[mm$forward] & !bR[mm$reverse],], mergers, bimFs, bimRs, SIMPLIFY=FALSE)
```

```python
>>> unqs.mock <- getUniques(mergers.nochim[["C-gunnisoni-BLS01-PDMHC-Amplicon"]])
>>> cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
DADA2 inferred 1 sample sequences present in the Mock community.
```

```python
>>> mockRef <- readFasta("6Alleles_alignment.fas")
>>> match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, as.character(sread(mockRef))))))
>>> cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
Of those, 1 were exact matches to the expected reference sequences.
```

```python
>>> mockRef
class: ShortRead
length: 6 reads; width: 271 cycles
```

```python
>>> seqtab <- makeSequenceTable(mergers.nochim[names(mergers.nochim)])
Duplicate sequences detected.
The sequences being tabled vary in length.
```

```python
>>> table(nchar(colnames(seqtab)))

240 249 256 266 271 287 298 
  3   5   1   1  10   1   1
```

```python
>>> for ( i in 1:22) {
...     seq.names <- c(seq.names, paste0("seq", i))
...     }
```

```python
>>> fasta.file = file('full_amplicon_seqs.fasta', "w")
...
>>> for (i in 1:length(colnames(seqtab))) {
...     id_ = paste0(">", seq.names[i])
...     write(id_, fasta.file, append=TRUE)
...     write(colnames(seqtab)[i], fasta.file, append=TRUE)
>>> }
>>> close(fasta.file)
```

```python
>>> test.tab <- seqtab
>>> colnames(test.tab) <- seq.names
>>> write.table(test.tab, "full_amplicon_table.txt", sep="\t", col.names = NA)
```

```python
>>> test.tab
                                  seq1 seq2 seq3 seq4 seq5 seq6 seq7 seq8 seq9
C-gunnisoni-BLS01-PDMHC-Amplicon  1273    0    0    0    0    0    0    0    0
C-gunnisoni-BLS02-PDMHC-Amplicon     0  947    0    0    0    0    0    0    0
C-gunnisoni-BLS03-PDMHC-Amplicon     0  936 1099    0    0    0    0    0    0
C-gunnisoni-BLS04-PDMHC-Amplicon     0 1044 1006    0    0    0    0    0    0
C-gunnisoni-BLS05-PDMHC-Amplicon     0 1575 1761    0    0    0    0    0    0
C-gunnisoni-BLS06-PDMHC-Amplicon   193  392    0    0    0    0    0    0    0
C-gunnisoni-CRL01-PDMHC-Amplicon     0    0 5138    0    0    0    0    0    0
C-gunnisoni-CRL02-PDMHC-Amplicon     0    0  575  625    0    0    0    0    0
C-gunnisoni-CRL03-PDMHC-Amplicon     0    0 1950    0    0    0    0    0    0
C-gunnisoni-CRL04-PDMHC-Amplicon     0    0  826  873    0    0    0    0    0
C-gunnisoni-CRL05-PDMHC-Amplicon     0    0 1158 1253    1    1    0    0    0
C-gunnisoni-CRL06-PDMHC-Amplicon     0    0 1061 1079    0    0    0    0    0
C-gunnisoni-CRL07-PDMHC-Amplicon     0    0 2220    0    0    0    0    0    0
C-gunnisoni-CRL08-PDMHC-Amplicon     0    0 2651    0    0    0    0    0    0
C-gunnisoni-CRL09-PDMHC-Amplicon     0    0 3311    0    0    0    0    0    0
C-gunnisoni-CRL10-PDMHC-Amplicon     0    0 2562    0    0    0    0    0    0
C-gunnisoni-CRL11-PDMHC-Amplicon     0    0 3063    0    0    0    0    0    0
C-gunnisoni-CRL12-PDMHC-Amplicon     0    0 3423    0    0    0    0    0    0
C-gunnisoni-CRL13-PDMHC-Amplicon     0    0 2759 2757    0    0    0    0    0
C-gunnisoni-CRL14-PDMHC-Amplicon     0    0 1233    0    0    0    0    0    0
C-gunnisoni-CRL15-PDMHC-Amplicon     0    0 3611    0    0    0    0    0    0
C-gunnisoni-CRL16-PDMHC-Amplicon     0    0 1893    0    0    0 3287    0    0
C-gunnisoni-CRL17-PDMHC-Amplicon     0    0 1171    0    0    0 2415    0    0
C-gunnisoni-CRL18-PDMHC-Amplicon     0    0  813    0    0    0    0   23    0
C-gunnisoni-CRL19-PDMHC-Amplicon     0    0 1476 1499    0    0    0    0   16
C-gunnisoni-CRL20-PDMHC-Amplicon     0    0 2176 2048    0    0    0    0    0
C-gunnisoni-CRL21-PDMHC-Amplicon     0    0  601  637    0    0    0    0    0
C-gunnisoni-CRL22-PDMHC-Amplicon     0    0 2732    0    0    0    0    0    0
C-gunnisoni-CRL23-PDMHC-Amplicon     0    0 1185 1310    0    0    0    0    0
C-gunnisoni-CRL24-PDMHC-Amplicon     0    0 1808    0    0    0    0    0    0
C-gunnisoni-ELMA01-PDMHC-Amplicon    0    0 1364    0    0    0    0    0    0
C-gunnisoni-ELMA02-PDMHC-Amplicon    0    0  589    0    0    0    0 1604    0
C-gunnisoni-ELMA04-PDMHC-Amplicon    0    0  833    0    0    0    0    0    0
C-gunnisoni-ELMA05-PDMHC-Amplicon    0    0 1160    0    0    0    0    0    0
C-gunnisoni-ELMA06-PDMHC-Amplicon    0    0  829    0    0    0    0    0    0
C-gunnisoni-ELMA07-PDMHC-Amplicon    0    0  691    0    0    0    0    0    0
C-gunnisoni-ELMA08-PDMHC-Amplicon    0    0 3389    0    0    0    0    0    0
C-gunnisoni-ELMA09-PDMHC-Amplicon    0    0  455    0    0    0    0    0    0
C-gunnisoni-ELMA10-PDMHC-Amplicon    0    0    0    0    0    0    0 1004    0
C-gunnisoni-HOVE01-PDMHC-Amplicon    0    0 1751    0    0    0 1920    0    0
C-gunnisoni-HOVE02-PDMHC-Amplicon    0    5    0    0    0    0 2476 2193    0
C-gunnisoni-HOVE03-PDMHC-Amplicon    0    0    6    0    0    0    6 2292 2379
C-gunnisoni-HOVE04-PDMHC-Amplicon    0    0  798    0    0    0 1553    0    0
C-gunnisoni-HOVE05-PDMHC-Amplicon    0    0 1480    0    0    0 2860    0    0
C-gunnisoni-HOVE06-PDMHC-Amplicon    0    0 1796    0    0    0    6    0    0
C-gunnisoni-HOVE07-PDMHC-Amplicon    0    4    0    0    0    0 1768 1654    0
C-gunnisoni-HOVE08-PDMHC-Amplicon    0    0 1187    0    0    0 2142    0    0
C-gunnisoni-HOVE09-PDMHC-Amplicon    0    0    0    0    0    0    0 3064 3073
C-gunnisoni-HOVE10-PDMHC-Amplicon    0    0 1260    0    0    0 2344    0    0
C-gunnisoni-HOVE11-PDMHC-Amplicon    0    0  934    0    0    0    5 1189    0
C-gunnisoni-HOVE12-PDMHC-Amplicon    0    0 1511    0    0    0    0 2491    0
C-gunnisoni-HOVE13-PDMHC-Amplicon    0    0 1609    0    0    0 2735    0    0
C-gunnisoni-HOVE14-PDMHC-Amplicon    0    0  997    0    0    0 2181    0    0
C-gunnisoni-HOVE15-PDMHC-Amplicon    0    0    0    0    0    0    0 3687    0
C-gunnisoni-HOVE16-PDMHC-Amplicon    0    0  765    0    0    0 1328    0    0
C-gunnisoni-HOVE18-PDMHC-Amplicon    0    0 1133    0    0    0    0    0    0
C-gunnisoni-HOVE19-PDMHC-Amplicon    0    0 1022    0    0    0    0 1755    0
C-gunnisoni-HOVE20-PDMHC-Amplicon    0    0  962    0    0    0    0    0    0
C-gunnisoni-HOVE21-PDMHC-Amplicon    0    0 4742    0    0    0    0    0    0
C-gunnisoni-HUTR01-PDMHC-Amplicon    0    0  510    0    0    0    0 1145    0
C-gunnisoni-HUTR02-PDMHC-Amplicon  912    0 1128    0    0    0    0    0    0
C-gunnisoni-HUTR04-PDMHC-Amplicon  272    0    0    0    0    0  797    0    0
C-gunnisoni-HUTR05-PDMHC-Amplicon    0    0  482    0    0    0  883    0    0
C-gunnisoni-HUTR06-PDMHC-Amplicon    0    0  270    0    0    0  412    0    0
C-gunnisoni-HUTR07-PDMHC-Amplicon    0    0  506    0    0    0    0    0    0
C-gunnisoni-HUTR08-PDMHC-Amplicon  202    0  271    0    0    0    0    0    0
C-gunnisoni-HUTR09-PDMHC-Amplicon    0    0  478    0    0    0  799    0    0
C-gunnisoni-HUTR11-PDMHC-Amplicon    0  726  753    0    0    0    0    0    0
C-gunnisoni-HUTR12-PDMHC-Amplicon    0    0  341    0    0    0  683    0    0
C-gunnisoni-HUTR14-PDMHC-Amplicon    0    0 2005    0    0    0    0    0    0
C-gunnisoni-HUTR15-PDMHC-Amplicon    0    0 2722    0    0    0    0    0    0
C-gunnisoni-HUTR16-PDMHC-Amplicon    0 2398    0    0    0    0    0    0    0
C-gunnisoni-HUTR17-PDMHC-Amplicon    0    0  962  934    0    0    0    0    0
C-gunnisoni-HUTR20-PDMHC-Amplicon    0  743  665    0    0    0    0    0    0
C-gunnisoni-HUTR21-PDMHC-Amplicon    0    0    0    0    0    0    0 3275    0
C-gunnisoni-HUTR22-PDMHC-Amplicon    0  698  738    0    0    0    0    0    0
C-gunnisoni-HUTR24-PDMHC-Amplicon    0 1020    0    0    0    0    0    0    0
C-gunnisoni-RSF01-PDMHC-Amplicon     0 1227 1476    0    0    0    0    0    0
C-gunnisoni-RSF02-PDMHC-Amplicon     0  378    0  567    0    0    0    0    0
C-gunnisoni-RSF03-PDMHC-Amplicon     0    0  906    0    0    0 1216    0    0
Cynomys-gunnisoni-105-55             0    0    0    0    0    0    0    0    0
Cynomys-gunnisoni-61-30              0    0    0    0    0    0    0    0    0
Cynomys-gunnisoni-74-35              0    0    0    0    0    0    0    0    0
Cynomys-gunnisoni-77-1               0    0    0    0    0    0    0    0    0
Cynomys-gunnisoni-80-34              0    0    0    0    0    0    0    0    0
                                  seq10 seq11 seq12 seq13 seq14 seq15 seq16
C-gunnisoni-BLS01-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-BLS02-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-BLS03-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-BLS04-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-BLS05-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-BLS06-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL01-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL02-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL03-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL04-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL05-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL06-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL07-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL08-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL09-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL10-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL11-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL12-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL13-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL14-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL15-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL16-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL17-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL18-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL19-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL20-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL21-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL22-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL23-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-CRL24-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-ELMA01-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-ELMA02-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-ELMA04-PDMHC-Amplicon  1173     0     0     0     0     0     0
C-gunnisoni-ELMA05-PDMHC-Amplicon  1715     0     0     0     0     0     0
C-gunnisoni-ELMA06-PDMHC-Amplicon  2291     0     0     0     0     0     0
C-gunnisoni-ELMA07-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-ELMA08-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-ELMA09-PDMHC-Amplicon   654     0     0     0     0     0     0
C-gunnisoni-ELMA10-PDMHC-Amplicon   778     0     0     0     0     0     0
C-gunnisoni-HOVE01-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE02-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE03-PDMHC-Amplicon     0     3     0     0     0     0     0
C-gunnisoni-HOVE04-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE05-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE06-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE07-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE08-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE09-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE10-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE11-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE12-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE13-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE14-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE15-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE16-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE18-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE19-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE20-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HOVE21-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR01-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR02-PDMHC-Amplicon     0     0     4     0     0     0     0
C-gunnisoni-HUTR04-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR05-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR06-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR07-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR08-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR09-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR11-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR12-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR14-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR15-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR16-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR17-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR20-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR21-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR22-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-HUTR24-PDMHC-Amplicon     0     0     0     0     0     0     0
C-gunnisoni-RSF01-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-RSF02-PDMHC-Amplicon      0     0     0     0     0     0     0
C-gunnisoni-RSF03-PDMHC-Amplicon      0     0     0     0     0     0     0
Cynomys-gunnisoni-105-55              0     0     0  8225  6159     0     0
Cynomys-gunnisoni-61-30               0     0     0     0     0  8311     2
Cynomys-gunnisoni-74-35               0     0     0     0  3281  3217     0
Cynomys-gunnisoni-77-1                0     0     0    21     0  3486     0
Cynomys-gunnisoni-80-34               0     0     0    65     0  5572     0
                                  seq17 seq18 seq19 seq20 seq21 seq22
C-gunnisoni-BLS01-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-BLS02-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-BLS03-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-BLS04-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-BLS05-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-BLS06-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL01-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL02-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL03-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL04-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL05-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL06-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL07-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL08-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL09-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL10-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL11-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL12-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL13-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL14-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL15-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL16-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL17-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL18-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL19-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL20-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL21-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL22-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL23-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-CRL24-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-ELMA01-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA02-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA04-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA05-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA06-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA07-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA08-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA09-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-ELMA10-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE01-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE02-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE03-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE04-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE05-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE06-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE07-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE08-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE09-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE10-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE11-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE12-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE13-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE14-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE15-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE16-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE18-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE19-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE20-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HOVE21-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR01-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR02-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR04-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR05-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR06-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR07-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR08-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR09-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR11-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR12-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR14-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR15-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR16-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR17-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR20-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR21-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR22-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-HUTR24-PDMHC-Amplicon     0     0     0     0     0     0
C-gunnisoni-RSF01-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-RSF02-PDMHC-Amplicon      0     0     0     0     0     0
C-gunnisoni-RSF03-PDMHC-Amplicon      0     0     0     0     0     0
Cynomys-gunnisoni-105-55              0     0     0     0     0     0
Cynomys-gunnisoni-61-30               2     0     0     0     0     0
Cynomys-gunnisoni-74-35               0     7     5     2     0     0
Cynomys-gunnisoni-77-1                0     0     0     0  4180     0
Cynomys-gunnisoni-80-34               0     0     0     0     0  5861
```

```python

```

```python
>>> colnames(seqtab)
 [1] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAGACACAACTATGGGGTGTTTGAG"                           
 [2] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTACAACCTGGAGGAGTACGTGCGCTTTGACAGCGACGTGGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTGGTGAG"                           
 [3] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTGGTGAG"                           
 [4] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTATGGGGTTGGTGAG"                           
 [5] "CCTTCTCCTTCTTTTCCCCCCCTCCCTTTTCCTTTCTCCCTTTTCCCATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTATGGGGTTGGTGAG"                           
 [6] "CCTTCTCCTTCTTTTCCCCCCCTCCCTTTTCCTTTCTCCCTTTTCCCATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTGGTGAG"                           
 [7] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGCGGTTCCTGGACAGATACTTCTACAACCGAGAGGAGCTCCTGCGCTTCGACAGCGACGTGGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCAGCAAGCCGAGAACTGGAACAGCCAGAAGGACCTCCTGGAGGATGAGCGGGCCGCGGTGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG"                           
 [8] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCACGGGCAGCGGCGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG"                           
 [9] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTGGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCACGGGCAGCAGCGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG"                           
[10] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTGGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCATGGGCAGCAGCGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG"                           
[11] "CCGGATCCTTCGTGTCCCCACAGCACGAGCGGCACCTCGCTGGTGTTGCGGTAGTTGGCCATCCGGCTGTAGTAGTGGTACACGGTCTGGAAACTGATCTCGTCCTCCAGGCTGAAGACAAATATGACAGCGTCCACCCACATAGCGAACTGCGGAGAGAGCGGGTGAGAGCTTCGCTGTGGAGCGGAGAGGAAGACAAGTCAAGTGCGTGGCTGACTGACTTAGTCTTGATCTCGTATG"                                                          
[12] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAGACACAACTATGGGGTGTTTGAGAGCTTCGCTGTGGAGCGGAGAGGTGAG"
[13] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAG"                                                 
[14] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTACAACCTGGAGGAGTACGTGCGCTTTGACAGCGACGTGGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCG"                                                 
[15] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCG"                                                 
[16] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTG"                                
[17] "CGACAGGCATCGGGCGACTACTACCAGGGGACGTGTCCACGTGGGACTCCTTGATGAATTGCGGCGGATTGTATGCGGGACGCATGAAAATGGGCAAGTAAAAATGACATAAACGTTGCGTTGACGCGCACGAGCGCGCGCGTCCGATCAGGCCCGAAAAGCCCATGAAATACGGGGATTCCTTTCGATTTGGTTTCCGCGCGAGTTCGCCTGTAAGCTTGCGCGCCGCTTGCATCGCGCCGAAGCTTTCGGCGTG"                                          
[18] "GACCTGCAGCAGGTATTCGACATTATCGTTCGGCATACGAAAGAAGTGATTCCCGGCCTCGATCGGCACGCGTTCCGGACGACCGACGCGCTGAAGGATCTCGGGGCGAACTCGCTCGATCGCTCGGAAATCGTCATCATGACGCTCGAATCGCTTTCGCTGCGCATTCCGCTCGTCGAGTTCGCCGGCGCGCGCAACATCGGCGAGCTGGCGGAAGTGCTGTATGCGAAGCTGAAGACA"                                                          
[19] "CCGCCGAATCCGATGCTCAATTTCAGGGCGCGCCGGATGTCCGACCCGCAGGTCGAGCGCACCCAGTCGAACCGCTCGTCGATCGGCCGGTCGAGATTGCGGCTCGGATGCAGGCGCCGCGCCCTCATTTGCAGCAGCGTCGCGATCAGCTCGACGATGCCCGCCGCGCTCAGGCCGTGCCCGACGAGCGACTTCGTCGTGTTGATGCGCGCGTGCGCGAGGCCGCTGCGGCCGAGCGCGTCGAGCTCGGTCGCGTCGCCGAGCGTCGAGCCGCTGCCGTGCGGATT"           
[20] "CTGACCGAACGATGGCTGGAGATACATGCCGAAATCGAAGCGTTGCGCCGCAAGGCGAAGGCGGCGGCGATCATGCAGATCGAACTGCTGATGACGCAATACGGCATCTCGCCGGGCGAACTGAAGGGCGTCGACGGTGCGCCGGCGCATCGCAAGCCGAAGGCCAGGTACTGGAATCCGGAGACGGGGCAGACGTGGTCGGGGCGCGGCCGGATGCCGAAATGGCTGGTCGGAAGGGAG"                                                          
[21] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCACGGGCAGCGGCGGACACGTTCTGCAG"                                                 
[22] "CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCAAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCG"
```

```python
>>> unqs.mock
CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAGACACAACTATGGGGTGTTTGAG 
                                                                                                                                                                                                                                                                           1273
```

```python
>>> cat("DADA2 inferred", length('CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAGACACAACTATGGGGTGTTTGAG'), "sample sequences present in the Mock community.\n")
DADA2 inferred 1 sample sequences present in the Mock community.
```

```python
>>> match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, as.character(sread(mockRef))))))
>>> cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
Of those, 1 were exact matches to the expected reference sequences.
```

```python
>>> mockRef <- readFasta("6Alleles_alignment.fas")
...
>>> for (seq in colnames(seqtab)) {
...     match.ref <- sum(sapply(seq, function(x) any(grepl(x, as.character(sread(mockRef))))))
...     print(match.ref)
>>> }
[1] 1
[1] 1
[1] 1
[1] 0
[1] 0
[1] 0
[1] 0
[1] 1
[1] 0
[1] 0
[1] 0
[1] 0
[1] 1
[1] 1
[1] 1
[1] 1
[1] 0
[1] 0
[1] 0
[1] 0
[1] 1
[1] 1
```

```python
>>> names(mockRef)
NULL
```

```python
>>> my_list = c('a', 'b', 'c', 'd')
```

```python
>>> my_list[1:2]
[1] "a" "b"
```

```python
>>> my_list[3]
[1] "c"
```

```python
>>> seqtab
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAGACACAACTATGGGGTGTTTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1273
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                              193
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                             912
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                             272
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                             202
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTACAACCTGGAGGAGTACGTGCGCTTTGACAGCGACGTGGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                              947
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                              936
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1044
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1575
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                              392
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               5
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               4
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                             726
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2398
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                             743
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                             698
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1020
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1227
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                              378
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1099
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1006
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1761
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                             5138
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                              575
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1950
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                              826
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1158
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1061
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2220
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2651
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                             3311
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2562
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                             3063
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                             3423
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2759
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1233
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                             3611
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1893
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1171
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                              813
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1476
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2176
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                              601
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2732
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1185
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1808
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1364
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                             589
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                             833
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1160
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                             829
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                             691
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                            3389
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                             455
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1751
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               6
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                             798
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1480
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1796
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1187
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1260
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                             934
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1511
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1609
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                             997
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                             765
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1133
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1022
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                             962
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                            4742
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                             510
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1128
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                             482
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                             270
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                             506
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                             271
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                             478
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                             753
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                             341
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2005
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2722
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                             962
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                             665
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                             738
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1476
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                              906
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTATGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                              625
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                              873
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1253
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1079
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2757
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1499
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2048
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                              637
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1310
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                             934
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                              567
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCTTCTCCTTCTTTTCCCCCCCTCCCTTTTCCTTTCTCCCTTTTCCCATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTATGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                1
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCTTCTCCTTCTTTTCCCCCCCTCCCTTTTCCTTTCTCCCTTTTCCCATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                1
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGCGGTTCCTGGACAGATACTTCTACAACCGAGAGGAGCTCCTGCGCTTCGACAGCGACGTGGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCAGCAAGCCGAGAACTGGAACAGCCAGAAGGACCTCCTGGAGGATGAGCGGGCCGCGGTGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                             3287
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                             2415
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1920
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2476
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               6
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1553
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2860
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               6
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1768
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2142
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2344
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               5
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2735
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2181
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1328
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                             797
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                             883
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                             412
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                             799
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                             683
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                             1216
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCACGGGCAGCGGCGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               23
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1604
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1004
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2193
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2292
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1654
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                            3064
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1189
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2491
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                            3687
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1755
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1145
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                            3275
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTGGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCACGGGCAGCAGCGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               16
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2379
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                            3073
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTGGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCATGGGCAGCAGCGGACACGTTCTGCAGACACAACTACGGGGTTGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1173
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                            1715
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                            2291
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                             654
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                             778
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                         0
                                  CCGGATCCTTCGTGTCCCCACAGCACGAGCGGCACCTCGCTGGTGTTGCGGTAGTTGGCCATCCGGCTGTAGTAGTGGTACACGGTCTGGAAACTGATCTCGTCCTCCAGGCTGAAGACAAATATGACAGCGTCCACCCACATAGCGAACTGCGGAGAGAGCGGGTGAGAGCTTCGCTGTGGAGCGGAGAGGAAGACAAGTCAAGTGCGTGGCTGACTGACTTAGTCTTGATCTCGTATG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                3
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                           0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                          0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAGACACAACTATGGGGTGTTTGAGAGCTTCGCTGTGGAGCGGAGAGGTGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          4
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                           0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                                                    0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                                                    0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                                                     0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                                                    0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCTTGAAGCAAGAGTCATTTGAGTGTCACTTTTCCAACGGGACGGAGCGGATACGGTTCCTACACAGATACATCTATAACCGGGAGGAGGTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGATGCTGAGTACTGGAACAGGCAGAAGGACATCCTGGAGCAGAGGCGGGCCGAGGTGGACACAGTGTGCAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                               8225
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                   21
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                  65
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCTACAACCTGGAGGAGTACGTGCGCTTTGACAGCGACGTGGGGGAGTACCGCGCGGTGACCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                               6159
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                3281
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                    0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                   0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                  0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                8311
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                3217
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                 3486
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                5572
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCGACACAACTACGGGGTTG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                           0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                    2
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                    0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                     0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                    0
                                  CGACAGGCATCGGGCGACTACTACCAGGGGACGTGTCCACGTGGGACTCCTTGATGAATTGCGGCGGATTGTATGCGGGACGCATGAAAATGGGCAAGTAAAAATGACATAAACGTTGCGTTGACGCGCACGAGCGCGCGCGTCCGATCAGGCCCGAAAAGCCCATGAAATACGGGGATTCCTTTCGATTTGGTTTCCGCGCGAGTTCGCCTGTAAGCTTGCGCGCCGCTTGCATCGCGCCGAAGCTTTCGGCGTG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                 0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                          2
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                           0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                          0
                                  GACCTGCAGCAGGTATTCGACATTATCGTTCGGCATACGAAAGAAGTGATTCCCGGCCTCGATCGGCACGCGTTCCGGACGACCGACGCGCTGAAGGATCTCGGGGCGAACTCGCTCGATCGCTCGGAAATCGTCATCATGACGCTCGAATCGCTTTCGCTGCGCATTCCGCTCGTCGAGTTCGCCGGCGCGCGCAACATCGGCGAGCTGGCGGAAGTGCTGTATGCGAAGCTGAAGACA
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                          7
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                           0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                          0
                                  CCGCCGAATCCGATGCTCAATTTCAGGGCGCGCCGGATGTCCGACCCGCAGGTCGAGCGCACCCAGTCGAACCGCTCGTCGATCGGCCGGTCGAGATTGCGGCTCGGATGCAGGCGCCGCGCCCTCATTTGCAGCAGCGTCGCGATCAGCTCGACGATGCCCGCCGCGCTCAGGCCGTGCCCGACGAGCGACTTCGTCGTGTTGATGCGCGCGTGCGCGAGGCCGCTGCGGCCGAGCGCGTCGAGCTCGGTCGCGTCGCCGAGCGTCGAGCCGCTGCCGTGCGGATT
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                                                               0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                                                                0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                                                        0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                                                         5
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                                                         0
                                  CTGACCGAACGATGGCTGGAGATACATGCCGAAATCGAAGCGTTGCGCCGCAAGGCGAAGGCGGCGGCGATCATGCAGATCGAACTGCTGATGACGCAATACGGCATCTCGCCGGGCGAACTGAAGGGCGTCGACGGTGCGCCGGCGCATCGCAAGCCGAAGGCCAGGTACTGGAATCCGGAGACGGGGCAGACGTGGTCGGGGCGCGGCCGGATGCCGAAATGGCTGGTCGGAAGGGAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                 0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                         0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                          2
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                           0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                          0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGGTTTCAGTGAGTGTCATTTCTCCAATGGGACGCAGCGGGTGCGGTTCCTGGAGAGACACTTCTACAACCAGGAGGAGTACCTGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCGAGCTGGGGCGGCGGTCAGCCAAGTACTTCAACAGCCAGAAGGACGCCCTGGAGCACGCACGGGCAGCGGCGGACACGTTCTGCAG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                  0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                 4180
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                   0
                                  CCGGATCCTTCGTGTCCCCACAGCACGTTTCCTGGAGCAAGTTTCACATGAGTGTCATTTCTCCAACGGGACGGAGCGGGTGCGGTTCCTGGACAGATACTTCCACAACCGAGAGGAGTTCGCGCGCTTCGACAGCGACGTCGGGGAGTACCGCGCGGTGAGCAAGCTGGGGCGGCCGGACGCCGAGTACTGGAACAGCCAGAAGGACCTCCTGGAGCAGAAGCGGGCCGCGGTGGACAACTACTGCCG
C-gunnisoni-BLS01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-BLS06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL04-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL05-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL06-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL07-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL08-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL09-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL10-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL11-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL12-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL13-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL14-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL15-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL16-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL17-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL18-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL19-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL20-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL21-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL22-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL23-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-CRL24-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-ELMA01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-ELMA10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE03-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE10-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE13-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE18-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE19-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HOVE21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR01-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR02-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR04-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR05-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR06-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR07-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR08-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR09-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR11-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR12-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR14-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR15-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR16-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR17-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR20-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR21-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR22-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-HUTR24-PDMHC-Amplicon                                                                                                                                                                                                                                                         0
C-gunnisoni-RSF01-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF02-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
C-gunnisoni-RSF03-PDMHC-Amplicon                                                                                                                                                                                                                                                          0
Cynomys-gunnisoni-105-55                                                                                                                                                                                                                                                                  0
Cynomys-gunnisoni-61-30                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-74-35                                                                                                                                                                                                                                                                   0
Cynomys-gunnisoni-77-1                                                                                                                                                                                                                                                                    0
Cynomys-gunnisoni-80-34                                                                                                                                                                                                                                                                5861
```

```python

```
