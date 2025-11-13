# nf-assemflow: Genome & Metagenome Assembly Pipeline (Nextflow)

**nf-assemflow** is a modular and reproducible **Nextflow pipeline** for assembling **genomes** or **metagenomes** from Illumina and/or long-read sequencing data.

It processes quality-controlled, host-removed paired-end or single-end reads (short and/or long), which can be generated with any NGS reads quality control and host removal pipelines (e.g., from **[nf-qcflow](https://github.com/xiaoli-dong/nf-qcflow)**), and integrates multiple assemblers — **SPAdes**, **SKESA**, **Unicycler**, **Shovill**, and **Flye** — with optional **reorientation** and **polishing** using **DNAapler**, **Medaka**, **Polypolish**, and **PyPolca**.

---

## Pipeline summary

- Supports **genome** and **metagenome** assembly workflows  
- Handles **Illumina-only**, **long-read-only**, or **hybrid (Illumina + long-read)** data  
- Assemblers:
  - [SPAdes](https://github.com/ablab/spades)
  - [SKESA](https://github.com/ncbi/SKESA)
  - [Unicycler](https://github.com/rrwick/Unicycler)
  - [Shovill](https://github.com/tseemann/shovill)
  - [Flye](https://github.com/fenderglass/Flye)
- **Post-assembly refinement:**
  - [DNAapler](https://github.com/widdowquinn/DNAAppler) — reorient circular genomes  
  - [Medaka](https://github.com/nanoporetech/medaka) — long-read polishing  
  - [Polypolish](https://github.com/rrwick/Polypolish) — short-read polishing  
  - [PyPolca](https://github.com/alekseyzimin/py-polca) — Illumina correction  
- Containerized (Docker / Singularity)
- Scalable for local, HPC, or cloud environments  
---

## Quick start

### Prepare required samplesheet input
The nf-assemflow pipeline requires user to provide a csv format samplesheet, which contains the sequenence information for each sample, as input. See below for what the samplesheet looks like:

```samplesheet.csv```

```
sample,fastq_1,fastq_2,long_fastq,basecaller_mode
sample1,shortreads_1.fastq.gz,shortreads_2.fastq.gz,longreads.fastq.gz,r1041_e82_400bps_hac_v4.2.0
sample2,shortreads.fastq,NA,longreads.fastq.gz,r1041_e82_400bps_sup_v4.2.0
sample3,NA,NA,longreads.fastq.gz,NA
sample4,shortreads_1.fastq.gz,shortreads_2.fastq.gz,NA
```
The csv format samplesheet has five required columns:
* The first row of the csv file is the header describing the columns
* Each row represents a unique sample to be processed, the first colum is the unique sample id
* When the information for a particular column is missing, please fill the column with "NA"
* The "fastq_1" and "fastq_2" columns are reserved for supplying the short sequence files
* "basecaller_mode" is for user to provide medaka inference model
* 
> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Run the pipeline:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run xiaoli-dong/nf-assemflow \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

xiaolid-ong/assemflow was originally written by Xiaoli Dong.
