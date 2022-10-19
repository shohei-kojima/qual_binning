# About
This is a simple code that converts a BAM/CRAM file to a lossy CRAM file with quality scoring by 8- or 4-binning.  

# Usage
- You need to [install Rust](https://www.rust-lang.org/tools/install) to use this code.  
- This code is only applicable for Phred+33 (i.e. Sanger and Illumina 1.8+ Phred scaled quality scoring).  

```
# compile
git clone https://github.com/shohei-kojima/qual_binning
cargo build --release

# help message and version
target/release/qual_binning -h
target/release/qual_binning --version

# typical use case, by default, 4-binning
target/release/qual_binning \
-i input.cram \
-o lossy.cram \
-r GRCh38DH.fa \
-p 4

# Illumina 8-binning requires "-e" option
target/release/qual_binning \
-i input.cram \
-o lossy.cram \
-r GRCh38DH.fa \
-e \
-p 4

# make symlink
ln -s target/release/qual_binning /usr/local/bin/qual_binning
```

# Binning

## 4-binning (default)
See details in below:  
https://topmed.nhlbi.nih.gov/sites/default/files/CORE_YR_1_WGS_CRAM_data_file_standardsv4.pdf
```
0, 0, 2, 3, 4, 5, 6, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20,
20, 20, 20, 20, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30,
30, 30, 30, 30, 30, 30
```

## 8-binning (`-e` option)
This is the one so called Illumina 8-bin. For detail, see below:  
https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf
```
0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 15, 15, 15, 15, 15, 15, 15, 15, 15,
15, 22, 22, 22, 22, 22, 27, 27, 27, 27, 27, 33, 33, 33, 33, 33, 37,
37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
40, 40, 40, 40, 40, 40
```

# Author
Shohei Kojima @ RIKEN


