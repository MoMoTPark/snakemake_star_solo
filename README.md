## Snakemake STARSolo 10X 3' and 5' single-cell RNA-seq data processing pipeline

This pipeline is designed to automate 10X derived single-cell data by processing by utilising STARSolo. Note the quantitative outputs are generated for the entire droplet set and therefore, filtering outlier droplets from droplets with real cells should be performed before downstream analysis. A working example code for performing cell calling downstream in R is provided in [Usage - Downstream cell calling](#downstream-cell-calling).

**Input:**

- Read 1 and 2 fastq files

**Output:**

- Mapped sequencing reads
- Seurat compatible gene and transcript counts matrix

**System requirements:**

- Snakemake v7.0 and above

### Usage

If `snakemake` is not installed on your local system, simplest way to install `snakemake` is by creating a new conda environment with an isolated Snakemake installation:

```
mamba create -n snakemake -c conda-forge -c bioconda snakemake;
mamba activate snakemake
```

Then in the activated `snakemake` environment execute the following command from [pipeline_root](#pipeline-directory-structure):

```
snakemake -j 16 --use-conda --rerun-incomplete
```

The pipeline's layout has been designed around Snakemake recommended best practices [Snakemake.readthedocs.io](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). To that extend, any supplied input file and resource to the pipeline can be provided as described below:

- `config/config.yaml`: All inputs are accessible in the pipeline from this config file. This includes reference sequences, annotations and other supplementary input files. Note that 10X barcode whitelist can be obtained from [CellRanger](https://github.com/10XGenomics/cellranger.git) by clonining the repository and going to the following path: `cellranger-x.x.x/lib/python/cellranger/barcodes`.

**Example config file:**
```
samples: config/samples.tsv
units: config/units.tsv
star_index: <star_index_dir_path>
ref_gtf: <reference_annotation_gtf_file_path>
bc_whitelist: <10x_barcode_whitelist_file_path>
```

#### Pipeline directory structure

User specific inputs should be provided via `config/` files. Input reference files (e.g., species specific annotation and genome fasta files) must be specified by user in `config/config.yaml`. Sample metadata including a `sample_id` must be supplied in `config/samples.tsv`. Sample input files (i.e., genome mapped sorted and indexed BAM) must be supplied in `config/units.tsv`. **Note that `sample_id` is used as primary key and therefore, must be the same value for a given sample in both `config/samples.tsv` and `config/units.tsv` files.**

Pipeline output files are stored in `results/` and each sample output is prefixed with the supplied `sample_id` value.

```
pipeline_root
├── workflow
│   ├── rules 
|   │   ├── commons.smk
|   │   └── starsolo.smk
│   ├── envs
│   ├── scripts
|   └── Snakefile
├── config
│   ├── config.yaml (modify config and provide required reference input files)
│   └── samples.tsv (add your sample meta data)
|   └── units.tsv (add your sample input file path)
├── results (pipeline output)
└── resources
```

#### Downstream cell calling

```
# Cell calling script with ont_pb_single_cell
library(DropletUtils)
library(Seurat)
library(Matrix)
library(tidyr)
library(dplyr)

# Count data import block
# character vector of input count directories to pass onto Read10X()
dir <- c(sample_1 = "",
         sample_2 = "")

# seurat object list to store all sample objects
so <- list()
# import expression data per sample
# NOTE inspect barcodes.tsv.gz prior to import
for (sample in names(dir)) {
  count_matrix <- Read10X(dir[[sample]], strip.suffix = TRUE)
  so[[sample]] <- CreateSeuratObject(count_matrix,
                                          assay = "RNA",
                                          project = sample,
                                          min.cells = 0,
                                          min.features = 0)
  remove(count_matrix)
}

# Cell calling block
counts = list()
edrops = list()
for (name in names(so)) {
  counts[[name]] <- LayerData(so[[name]], assay='RNA', layer='counts')
  edrops[[name]] <- emptyDrops(counts[[name]],
                               lower=500)
}
solo_edrops = list()
# Provide n.expected.cells parameter based on experiment
for (name in names(so)) {
  solo_edrops[[name]] <- emptyDropsCellRanger(counts[[name]], 
                                              n.expected.cells=5000, 
                                              max.percentile=0.95,
                                              umi.min = 250)
}

# Create Seurat object with count matrix of real cells
for (name in names(solo_edrops)){
  cc <- solo_edrops[[name]]
  cells <- rownames(cc)[which(cc[['FDR']] <= 0.05)]
  counts_mtx <- counts[[name]][,cells]
  so[[name]] <- CreateSeuratObject(counts_mtx,
                                          assay = "RNA",
                                          project = name,
                                          min.cells = 0,
                                          min.features = 0)
}
```