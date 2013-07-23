# Amazon EC2 AMI setup notes

We start with the Ubuntu 13.04 (Raring Ringtail) EBS-backed AMI for the us-east-1 region, selected from the [Ubuntu Amazon EC2 AMI Locator ](http://cloud-images.ubuntu.com/locator/ec2) page.

## Boot AMI

AMI was started as type **m3.xlarge** with a 120GB EBS boot volume. This instance type has 4 CPU cores and ~16GB of memory. To see all instance types, see [EC2Instances.info](http://www.ec2instances.info/).

## Set up data directory

```bash
sudo mkdir /data
sudo chown -R ubuntu:ubuntu /data
```

## Install R

Create /etc/apt/sources.list.d/cran.list and add the line:

```
deb http://cran.revolutionanalytics.com/bin/linux/ubuntu raring/
```

Execute the following to add the appropriate signing key:

```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
```

Update packages and install R:

```bash
sudo apt-get update
sudo apt-get install r-base r-base-dev
sudo apt-get install openjdk-7-jdk
sudo R CMD javareconf 
sudo apt-get install libxml2-dev libcurl3-dev # for Rsamtools
sudo apt-get install htop samtools pigz fastx-toolkit sra-toolkit parallel git python-numpy python-dev
```

Create ~/.Rprofile:

```S
options(repos=c(CRAN="http://cran.rstudio.com"))

q <- function(save = "no", status = 0, runLast = TRUE) {
  .Internal(quit(save, status, runLast))
}

update_bioc = function() {
  library(BiocInstaller)
  update.packages(repos=biocinstallRepos(), ask=FALSE)
}

update_all = function() {
  message("Updating R packages...")
  update.packages(ask=F)
  message("Updating Bioconductor packages...")
  update_bioc()
}
````

Install R and Bioconductor packages by sourcing the following script from a root R session:

```S
update.packages(ask=F)
install.packages(c("gtools", "reshape", "ggplot2", "xtable", 
                   "xlsx", "fastcluster", "data.table", "dataframe",
                   "matrixStats", "stringr", "optparse", "xts",
                   "knitr", "knitcitations", "gridExtra"))

source("http://bioconductor.org/biocLite.R")
biocLite()

bioc_packages <- c("GenomicRanges", "BSgenome.Dmelanogaster.UCSC.dm3", "rtracklayer", 
                   "simpleaffy", "Rsamtools", "chipseq", "org.Dm.eg.db", "GOstats",
                   "seqLogo")
biocLite(bioc_packages)

install.packages("Vennerable", repos="http://R-Forge.R-project.org")
````

## Download analysis code

```bash
cd /home/ubuntu
git clone https://github.com/zeitlingerlab/chen_elife_2013.git analysis_code
```

## Install Bowtie v0.12.8

```bash
mkdir /data/download
cd /data/download
wget 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie/0.12.8/bowtie-0.12.8-linux-x86_64.zip'
unzip bowtie*
mkdir /data/apps
mv bowtie-0.12.8 /data/apps/bowtie

```

## Install Bowtie 2.0.2

```bash
cd /data/download
wget 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.0.2/bowtie2-2.0.2-linux-x86_64.zip'
unzip bowtie2*
mv bowtie2-2.0.2 /data/apps/bowtie2
```

## Add both bowtie versions to path

Modify ~/.profile:

```bash
export PATH=/data/apps/bowtie:/data/apps/bowtie2:$PATH
````

## Install Tophat v2.0.6 and Cufflinks v2.0.2

```bash
cd /data/download
wget 'http://tophat.cbcb.umd.edu/downloads/tophat-2.0.6.Linux_x86_64.tar.gz'
tar xf tophat*
mv tophat-2.0.6.Linux_x86_64 /data/apps/tophat

wget 'http://cufflinks.cbcb.umd.edu/downloads/cufflinks-2.0.2.Linux_x86_64.tar.gz'
tar xf cufflinks*
mv cufflinks-2.0.2.Linux_x86_64 /data/apps/cufflinks
````

## Install MACS

```bash
cd /data/download
git clone https://github.com/taoliu/MACS.git
cd MACS
git checkout v2.0.10_6_6_2012
sudo python setup.py install
```

## Install MEME 4.8.1

```bash

sudo apt-get install libhtml-template-perl libxml-simple-perl libsoap-lite-perl imagemagick

cd /data/download
wget 'ftp://ftp.ebi.edu.au/pub/software/MEME/4.8.1/meme_4.8.1.tar.gz'
tar xf meme_4.8.1.tar.gz
cd meme_4.8.1
./configure --prefix=/data/apps/meme --enable-build-libxml2 --with-url=http://meme.nbcr.net/meme 
make -j 4
make install
```

## Build UCSC dm3 reference genome for Bowtie

```bash
cd ~/download
mkdir dm3
cd dm3
wget 'http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz'
tar xf chromFa.tar.gz
cat *.fa > dm3.fa
/data/apps/bowtie/bowtie-build -o 1 dm3.fa dm3
mkdir /data/dm3_genome
mv dm3.* /data/dm3_genome
```

## Build FlyBase r5.47 reference genome for Bowtie2

```bash
cd ~/download
mkdir flybase
cd flybase
wget 'ftp://ftp.flybase.net/releases/FB2012_05/dmel_r5.47/fasta/dmel-all-chromosome-r5.47.fasta.gz'
gunzip dmel-all-chromosome-r5.47.fasta.gz 
/data/apps/bowtie2/bowtie2-build -o 1 dmel-all-chromosome-r5.47.fasta fb_547_genome
mkdir /data/tophat_cufflinks/genome
mv fb_547_genome.* /data/tophat_cufflinks/genome
```

## Download Eisen RNA-seq data

Create fetch\_rnaseq.sh in /data/eisen_rnaseq:

```bash
#!/bin/bash

for f in 10 11 12 13 14A 14B 14C 14D
do
	female_file="F${f}_r1_reads.fq"
	male_file="M${f}_r1_reads.fq"
	wget http://www.eisenlab.org/dosage/data/Reads/RNAseq/$female_file
	pigz -p 4 $female_file
	wget http://www.eisenlab.org/dosage/data/Reads/RNAseq/$male_file
	pigz -p 4 $male_file
done
````

## Download modENCODE whole-embryo 04to06h and 06to08h RNA-seq data

```bash
mkdir /data/modencode_rnaseq
cd /data/modencode_rnaseq
# SRX008027: 04-06h embryos
wget --mirror ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX008/SRX008027
# SRX008025: 06-08h embryos
wget --mirror ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX008/SRX008025
mkdir we04to06h
mkdir we06to08h
mv `find ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX008/SRX008027 -name "*.sra"` we04to06h
mv `find ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX008/SRX008025 -name "*.sra"` we06to08h
cd we04to06h
parallel -uj 4 fastq-dump --split-files --gzip {} -- *.sra
cd ../we06to08h
parallel -uj 4 fastq-dump --split-files --gzip {} -- *.sra

# remove SRA files to save space
rm -f /data/modencode_rnaseq/*/*.sra
```

## ChIP-seq alignment

Align all single-end samples:

```bash
cd /data/fastq
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(fastx_trimmer -Q 33 -l 40 -i <(gunzip -c preMBT_k4_1.fastq.gz)) | samtools view -F 4 -Sbo preMBT_k4_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_k4_2.fastq.gz) | samtools view -F 4 -Sbo preMBT_k4_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_pol_1.fastq.gz) | samtools view -F 4 -Sbo preMBT_pol_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_pol_2.fastq.gz) | samtools view -F 4 -Sbo preMBT_pol_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(fastx_trimmer -Q 33 -l 40 -i <(gunzip -c preMBT_pol_3.fastq.gz)) | samtools view -F 4 -Sbo preMBT_pol_3.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_pol_4.fastq.gz) | samtools view -F 4 -Sbo preMBT_pol_4.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_tbp_1.fastq.gz) | samtools view -F 4 -Sbo preMBT_tbp_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_tbp_2.fastq.gz) | samtools view -F 4 -Sbo preMBT_tbp_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c preMBT_wce_1.fastq.gz) | samtools view -F 4 -Sbo preMBT_wce_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(fastx_trimmer -Q 33 -l 40 -i <(gunzip -c MBT_h3ac_1.fastq.gz)) | samtools view -F 4 -Sbo MBT_h3ac_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_h3ac_2.fastq.gz) | samtools view -F 4 -Sbo MBT_h3ac_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(fastx_trimmer -Q 33 -l 40 -i <(gunzip -c MBT_k27_1.fastq.gz)) | samtools view -F 4 -Sbo MBT_k27_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_k4_1.fastq.gz) | samtools view -F 4 -Sbo MBT_k4_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_k4_2.fastq.gz) | samtools view -F 4 -Sbo MBT_k4_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_pol_1.fastq.gz) | samtools view -F 4 -Sbo MBT_pol_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_pol_2.fastq.gz) | samtools view -F 4 -Sbo MBT_pol_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(fastx_trimmer -Q 33 -l 40 -i <(gunzip -c MBT_pol_3.fastq.gz)) | samtools view -F 4 -Sbo MBT_pol_3.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_tbp_1.fastq.gz) | samtools view -F 4 -Sbo MBT_tbp_1.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_tbp_2.fastq.gz) | samtools view -F 4 -Sbo MBT_tbp_2.bam -
/data/apps/bowtie/bowtie  -S -p 8 -m 1 -l 40 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c MBT_wce_1.fastq.gz) | samtools view -F 4 -Sbo MBT_wce_1.bam -
/data/apps/bowtie/bowtie --solexa1.3-quals -S -p 8 -m 1 -l 18 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c postMBT_k27_1.fastq.gz) | samtools view -F 4 -Sbo postMBT_k27_1.bam -
/data/apps/bowtie/bowtie --solexa1.3-quals -S -p 8 -m 1 -l 18 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c postMBT_pol_1.fastq.gz) | samtools view -F 4 -Sbo postMBT_pol_1.bam -
/data/apps/bowtie/bowtie --solexa1.3-quals -S -p 8 -m 1 -l 18 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c postMBT_tbp_1.fastq.gz) | samtools view -F 4 -Sbo postMBT_tbp_1.bam -
/data/apps/bowtie/bowtie --solexa1.3-quals -S -p 8 -m 1 -l 18 -k 1 -n 2 --best --strata /data/dm3_genome/dm3 <(gunzip -c postMBT_wce_1.fastq.gz) | samtools view -F 4 -Sbo postMBT_wce_1.bam -
````

Align all paired-end samples:

```bash
# K27 postMBT (fragment size 146)
/data/apps/bowtie/bowtie -S -p 8 -m 1 -l 50 -k 1 -n 2 --best --strata --chunkmbs 1024 -X 300 -1 <(gunzip -c 06to08h_k27_1_read1.fastq.gz) -2 <(gunzip -c 06to08h_k27_1_read2.fastq.gz) /data/dm3_genome/dm3 | samtools view -F 4 -Sbo 06to08h_k27_1.bam -

# MNase MBT (fragment size 146)
/data/apps/bowtie/bowtie -S -p 8 -m 1 -l 50 -k 1 -n 2 --best --strata --chunkmbs 1024 -X 300 -1 <(gunzip -c MBT_mnase_1_read1.fastq.gz) -2 <(gunzip -c MBT_mnase_1_read2.fastq.gz) /data/dm3_genome/dm3 | samtools view -F 4 -Sbo MBT_mnase_1.bam -
```

## Align Eisen RNA-seq samples

Create /data/tophat\_cufflinks/run\_tophat\_for\_directory.sh:

```bash
#!/bin/bash

cd $1
fastq_files=`ls -1 *.fq.gz | paste -s -d,`
/data/apps/tophat/tophat -G /data/tophat_cufflinks/fb-r5.47.gtf -p 8 -I 5000 -z pigz --segment-length=20 -o tophat \
                    --transcriptome-index=/data/tophat_cufflinks/tophat_index/index \
                    /data/tophat_cufflinks/genome/fb_547_genome \
                    $fastq_files
cd /data/tophat_cufflinks
```

Run the file for all Eisen RNA-seq samples:

```bash
cd /data/tophat_cufflinks
for sample in 10 11 12 13 14A 14B 14C 14D
do
  mkdir $sample
  cd $sample
  ln -s /data/eisen_rnaseq/?${sample}_r1_reads.fq.gz .
  cd ..
done

./run_tophat_for_directory.sh 10
./run_tophat_for_directory.sh 11
./run_tophat_for_directory.sh 12
./run_tophat_for_directory.sh 13
./run_tophat_for_directory.sh 14A
./run_tophat_for_directory.sh 14B
./run_tophat_for_directory.sh 14C
./run_tophat_for_directory.sh 14D
````

## Align modENCODE RNA-seq samples (paired-end)

Create /data/tophat\_cufflinks/run\_pe_tophat\_for\_directory.sh:

```bash
#!/bin/bash

cd $1
fastq_files_r1=`ls -1 *_1.fastq.gz | paste -s -d,`
fastq_files_r2=`ls -1 *_2.fastq.gz | paste -s -d,`
/data/apps/tophat/tophat -G /data/tophat_cufflinks/fb-r5.47.gtf -p 8 -I 5000 -z pigz --segment-length=38 -o tophat \
                    --transcriptome-index=/data/tophat_cufflinks/tophat_index/index \
                    /data/tophat_cufflinks/genome/fb_547_genome \
                    $fastq_files_r1 $fastq_files_r2
cd /data/tophat_cufflinks
```

Run the file for both modENCODE samples:

```bash
cd /data/tophat_cufflinks
for sample in we04to06h we06to08h
do
  mkdir $sample
  cd $sample
  ln -s /data/modencode_rnaseq/$sample/*.fastq.gz .
  cd ..
done

./run\_pe\_tophat\_for\_directory.sh we04to06h
./run\_pe\_tophat\_for\_directory.sh we06to08h
```

## Process RNA-seq alignments with Cufflinks

Run cuffdiff on the Eisen RNA-seq samples:

```bash
cd /data/tophat_cufflinks
/data/apps/cufflinks/cuffdiff -p 8 -u -b genome/fb_547_genome.fa -L cc10,cc11,cc12,cc13,cc14a,cc14b,cc14c,cc14d \
    -o cuffdiff2 \
    fb-r5.47.gtf \
    10/tophat/accepted_hits.bam \
    11/tophat/accepted_hits.bam \
    12/tophat/accepted_hits.bam \
    13/tophat/accepted_hits.bam \
    14A/tophat/accepted_hits.bam \
    14B/tophat/accepted_hits.bam \
    14C/tophat/accepted_hits.bam \
    14D/tophat/accepted_hits.bam 
gzip cuffdiff2/*
```

Run cufflinks individually on the two whole-embryo RNA-seq samples:

```bash
cd /data/tophat_cufflinks/we04to06h
/data/apps/cufflinks/cufflinks -p 8 --compatible-hits-norm -I 5000 \
	-G /data/tophat_cufflinks/fb-r5.47.gtf \
	-u \
	-b /data/tophat_cufflinks/genome/fb_547_genome.fa \
	-o cufflinks2 \
	tophat/accepted_hits.bam
gzip cufflinks2/*

cd /data/tophat_cufflinks/we06to08h
/data/apps/cufflinks/cufflinks -p 8 --compatible-hits-norm -I 5000 \
	-G /data/tophat_cufflinks/fb-r5.47.gtf \
	-u \
	-b /data/tophat_cufflinks/genome/fb_547_genome.fa \
	-o cufflinks2 \
	tophat/accepted_hits.bam
gzip cufflinks2/*
```

## R object generation for ChIP-seq samples

Generate R objects for all aligned samples:

```bash
cd /data/fastq
Rscript ~/scripts/bamtoolsr.r -f preMBT_k4_1.bam   -e 123 -n preMBT_k4_1
Rscript ~/scripts/bamtoolsr.r -f preMBT_k4_2.bam   -e 115 -n preMBT_k4_2
Rscript ~/scripts/bamtoolsr.r -f preMBT_pol_1.bam  -e 106 -n preMBT_pol_1
Rscript ~/scripts/bamtoolsr.r -f preMBT_pol_2.bam  -e 127 -n preMBT_pol_2
Rscript ~/scripts/bamtoolsr.r -f preMBT_pol_3.bam  -e 118 -n preMBT_pol_3
Rscript ~/scripts/bamtoolsr.r -f preMBT_pol_4.bam  -e 150 -n preMBT_pol_4
Rscript ~/scripts/bamtoolsr.r -f preMBT_tbp_1.bam  -e 91  -n preMBT_tbp_1
Rscript ~/scripts/bamtoolsr.r -f preMBT_tbp_2.bam  -e 110 -n preMBT_tbp_2
Rscript ~/scripts/bamtoolsr.r -f preMBT_wce_1.bam  -e 146 -n preMBT_wce_1
Rscript ~/scripts/bamtoolsr.r -f MBT_h3ac_1.bam    -e 90  -n MBT_h3ac_1
Rscript ~/scripts/bamtoolsr.r -f MBT_h3ac_2.bam    -e 131 -n MBT_h3ac_2
Rscript ~/scripts/bamtoolsr.r -f MBT_k27_1.bam     -e 120 -n MBT_k27_1
Rscript ~/scripts/bamtoolsr.r -f MBT_k4_1.bam      -e 84  -n MBT_k4_1
Rscript ~/scripts/bamtoolsr.r -f MBT_k4_2.bam      -e 113 -n MBT_k4_2
Rscript ~/scripts/bamtoolsr.r -f MBT_pol_1.bam     -e 161 -n MBT_pol_1
Rscript ~/scripts/bamtoolsr.r -f MBT_pol_2.bam     -e 113 -n MBT_pol_2
Rscript ~/scripts/bamtoolsr.r -f MBT_pol_3.bam     -e 95  -n MBT_pol_3
Rscript ~/scripts/bamtoolsr.r -f MBT_tbp_1.bam     -e 173 -n MBT_tbp_1
Rscript ~/scripts/bamtoolsr.r -f MBT_tbp_2.bam     -e 86  -n MBT_tbp_2
Rscript ~/scripts/bamtoolsr.r -f MBT_wce_1.bam     -e 133 -n MBT_wce_1
Rscript ~/scripts/bamtoolsr.r -f postMBT_k27_1.bam -e 147 -n postMBT_k27_1
Rscript ~/scripts/bamtoolsr.r -f postMBT_pol_1.bam -e 78  -n postMBT_pol_1
Rscript ~/scripts/bamtoolsr.r -f postMBT_tbp_1.bam -e 81  -n postMBT_tbp_1
Rscript ~/scripts/bamtoolsr.r -f postMBT_wce_1.bam -e 113 -n postMBT_wce_1

# paired-end
Rscript ~/scripts/bamtoolsr.r -f 06to08h_k27_1.bam -p     -n 06to08h_k27_1
Rscript ~/scripts/bamtoolsr.r -f MBT_mnase_1.bam   -p     -n MBT_mnase_1

# move files
mkdir /data/bam
mkdir /data/rdata
mkdir /data/bigwigs

mv /data/fastq/*.bam /data/bam
mv /data/fastq/*.bw /data/bigwigs
mv /data/fastq/*.RData /data/rdata

```

## Generate enrichment tracks

```bash
mkdir /data/rdata/enrichment
cd /data/rdata/enrichment
parallel -uj 1 -a ~/analysis_code/enrichment_tracks/build_enrichment_tracks.txt 
~/analysis_code/enrichment_tracks/build_minimum_enrichment_coverage.sh
```

## Run MACS

```bash
mkdir /data/macs
cd /data/macs

# pre-MBT Pol II replicates
parallel -uj 2 macs2 callpeak -t /data/bam/preMBT_pol_{}.bam -c /data/bam/preMBT_wce_1.bam -g dm -n preMBT_pol_{} -- 1 2 3 4

# pre-MBT TBP replicates
parallel -uj 2 macs2 callpeak -t /data/bam/preMBT_tbp_{}.bam -c /data/bam/preMBT_wce_1.bam -g dm -m 2 50 -n preMBT_tbp_{} -- 1 2
```

