# Procrustes-based ANcestry Determination and Assignment (PANDA)

This pipeline classifies targested next generation sequencing samples into gnomAD v3 population groups from GRCh37 BAM or CRAM files.

It does the following:

1. Creates pileup files at selected reference sites using `samtools mpileup`.
2. Converts pileup files into LASER `.seq` format.
3. Runs LASER to project samples into reference PC space.
4. Combines LASER `.coord` files.
5. Runs an R random forest model to assign population labels.

## Dependencies

The following command-line tools are required:

- `samtools`
- GNU `parallel`
- `Rscript`
- `python2`
- `tar`
- `gzip`
- `curl` or `wget`

The R package `randomForest` is also required.

Install it in R with:

```r
install.packages("randomForest", dependencies = TRUE)
```

Example conda installation for command-line dependencies:

```bash
conda install -c bioconda samtools parallel
conda install -c conda-forge r-base
conda install -c conda-forge python=2.7
```

## Resources

Required resources are downloaded automatically by the wrapper script when needed.

These include:

- the custom reference files and classification model;
- LASER software;
- the GRCh37 reference FASTA.

To download and check the resources without running the full pipeline:

```bash
./run_pipeline.sh --setup-only
```

To force resources to be downloaded again:

```bash
./run_pipeline.sh --setup-only --force-download
```

The reference FASTA index is also created automatically if it is missing.

## Input

The input is a text file containing one BAM or CRAM file path per line.

Example `bams.txt`:

```text
/path/to/sample1.bam
/path/to/sample2.bam
/path/to/sample3.cram
```

The BAM or CRAM files should already be indexed.

Input files must be aligned to GRCh37. If your BAM or CRAM files are aligned to another genome build, realign the reads to GRCh37 before running this pipeline. Ideally, realign from the original FASTQ files. If only BAM/CRAM files are available, first convert them back to FASTQ and then align to GRCh37:

```bash
samtools sort -n -o sample.name_sorted.bam sample.input.bam

samtools fastq \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -0 /dev/null \
  -s sample_singletons.fastq.gz \
  -n \
  sample.name_sorted.bam

bwa mem resources/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa \
  sample_R1.fastq.gz sample_R2.fastq.gz | \
  samtools sort -o sample.GRCh37.bam

samtools index sample.GRCh37.bam
```

## Usage

Basic usage:

```bash
./run_pipeline.sh --input bams.txt --output cohort1
```

Run with multiple parallel jobs:

```bash
./run_pipeline.sh --input bams.txt --output cohort1 --threads 16
```

Keep temporary files:

```bash
./run_pipeline.sh --input bams.txt --output cohort1 --threads 16 --keep-temp
```

## Options

```text
-i, --input FILE        Text file containing one BAM/CRAM path per line
-o, --output PREFIX     Output prefix
-t, --threads INT       Number of parallel jobs to use
--threshold FLOAT       Classification probability threshold
--setup-only            Download/check resources and exit
--force-download        Re-download resources even if files already exist
--keep-temp             Keep temporary files after the run
--resource-dir DIR      Alternative resources directory
--bin-dir DIR           Alternative binary/software directory
```

## Output

The pipeline writes:

```text
PREFIX.coord
PREFIX.population.classification.txt
unprocessed_ids.txt
PREFIX.make_pileup.log
```

`PREFIX.coord` contains the combined LASER coordinates.

`PREFIX.population.classification.txt` contains the predicted population labels and class probabilities.

`unprocessed_ids.txt` contains samples that could not be classified because of missing PC data.

`PREFIX.make_pileup.log` contains the GNU parallel job log from the pileup step.