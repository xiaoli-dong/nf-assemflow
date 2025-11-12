# nf-assemflow: Genome & Metagenome Assembly Pipeline (Nextflow)

**nf-assemflow** is a modular and reproducible **Nextflow pipeline** for assembling **genomes** or **metagenomes** from Illumina and/or long-read sequencing data.

It processes quality-controlled, host-removed paired-end or single-end reads (short and/or long), which can be generated using **[nf-qcflow](https://github.com/xiaoli-dong/nf-qcflow)**, and integrates multiple assemblers — **SPAdes**, **SKESA**, **Unicycler**, **Shovill**, and **Flye** — with optional **reorientation** and **polishing** using **DNAapler**, **Medaka**, **Polypolish**, and **PyPolca**.

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
- Integrated **MultiQC** and **QUAST** reporting  

---

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

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
