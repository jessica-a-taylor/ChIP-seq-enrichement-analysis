# nextflow-peakcaller
Peak calling pipeline following the method used in the ReMap database.

## Installation

This pipeline is implimented as a nextflow pipeline and setup to automatically install
all needed dependencies with conda. Therefore, the recommened way to install it is by installing
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) in a Unix
environment, then using conda to install Nextflow into an envrionment using the command:

```sh
$ conda create -c bioconda -n nextflow nextflow
```

There are also a number of non-Python dependencies, which can be installed with:

```sh
$ sudo apt install trim-galore bowtie2 samtools bedtools python3-biopython macs sra-toolkit rclone
```

The Nextflow environment can then be activated by running:

```sh
$ conda activate nextflow
```

You may also need to add execution permissions to the executable files included in the repository:

```sh
$ chmod +x bin/getFastaLength.py
$ chmod +x ./run_analysis.sh
```

## Usage with `run_analysis.sh`

This version of the repository includes a shell script which, in addition to running the Nextflow pipeline, will download the relevant data given the correct experiment names and pass them into the pipeline for you. It is also capable of uploading results to a cloud storage provider using the `rclone` Linux package. To execute the script, simply launch it with:

```sh
$ ./run_analysis.sh
```

It will prompt you for the information it needs, download and extract datasets using `prefetch` and `fastq-dump`, execute the pipeline, tidy up any large files left by the execution and upload the result to the storage provider using `rclone`. Multiple datasets can be specified by separating the experiment names with a space, and the experiment names should look like "SRR22894115". You will still need to provide a reference genome as explained in the next section.

You will also likely need to modify the variables at the top of the script to match your setup.

### Cloud upload

In order to configure cloud upload, you will need to add a remote to `rclone` if you don't have one already. This can be done with:

```sh
$ rclone config
```

Alternatively, the cloud upload step can be skipped by leaving the relevant variables (`PATH_TO_LOCAL_OUTPUT_DIR_FOR_UPLOAD` and `PATH_TO_CLOUD_UPLOAD_DIR`) blank in `run_analysis.sh`.

## Usage without `run_analysis.sh`

To call peaks you need a reference genome in FASTA format (`genome.fasta` in the example) and one or more directories containing FASTQ files with ChIP-seq reads 
against an interesting target (`/path/to/H3K4me3` and `/path/to/H3K9me3` in the example). Then you run the pipeline by passing the 
genome to `--genome-fasta` and the directories containing the ChIP-seq FASTQ files as a comma separated list to `--chip-seq-fastq`. For example:

```
$ nextflow run callPeaks.nf --genome-fasta genome.fasta --chip-seq-fastq /path/to/H3K4me3,/path/to/H3K9me3
```

If you want to run on paired end data, you need both end's FASTQ files in the same directory with the same name, except one ends with `1.fastq` and the other ends with `2.fastq`. This can then be run with the `--paired` option. (e.g. as follows)

```
$ nextflow run callPeaks.nf --genome-fasta genome.fasta --chip-seq-fastq /path/to/H3K4me3,/path/to/H3K9me3 --paired
```

Finally, to deal with cases where a ChIP-seq input file is provided to reduce false positives, the directory containing this can be specified in the pipeline with `--control-fastq` e.g.:

```
$ nextflow run callPeaks.nf --genome-fasta genome.fasta --chip-seq-fastq /path/to/H3K4me3,/path/to/H3K9me3 --control-fastq /path/to/Input
```

These will both call peaks separately for all sets of FASTQ files within the directories `H3K4me3` and `H3K9me3`. There are a couple of other parameters that can be set if needed:

- `--threads` sets the number of CPU cores to use, defaults to 8

- `--output-dir` sets the name of the output directory produced, defaults to `output/`

- `--use-rmdup` changes the tool to remove PCR duplicates to `samtools rmdup`. This was used in the orignal REMAP pipeline, but has now been deprecated, so isn't recommended for use, unless exactly replicating the pipeline's output.

- `--thread-samtools` adds a threading argument to all samtools commands, which may make them run faster, but requires a lot of RAM

## Output Format

For each input file provided, two output directories are produced containing output in the same layout as the [MACS2](https://github.com/macs3-project/MACS)
peakcaller. 

- `{filename}peakCallingResultsNarrow` contains narrow peaks called by MACS
- `{filename}peakCallingResultsBroad` contains broad peaks called by MACS
