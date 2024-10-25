### Snakemake STARSolo 10X 3' and 5' single-cell RNA-seq data processing pipeline



#### Directory structure

Initial directory structure is the minimum requirement for a bug-free pipeline execution; however, more complex pipelines can be implemented based on expanding current directory structure.

```
# Example directory structure
├── .gitignore
├── README.md
├── LICENSE.md
├── workflow
│   ├── rules
|   │   ├── commons.smk
|   │   └── test.smk
│   ├── envs
|   │   ├── bamstats.sif
|   │   ├── samtools_env.yaml
|   │   └── global_env.yaml
│   ├── scripts
|   │   └── txt_to_fasta.py
|   └── Snakefile
├── config
│   ├── config.yaml
│   └── samples.tsv
|   └── units.tsv
├── results
└── resources
```

#### How to run

The run command depends on the use case of the pipeline. For example if a Singularity container with a bash command is used the full directory path for the output should be passed onto Snakemake for mapping local directory to container file system. The following command would run a test pipeline and generates a few test text file outputs:  
`snakemake -j 2 --rerun-incomplete --use-conda --use-singularity --singularity-args "--bind <result's dir full path>"`

**Note:** Full path of input files should be supplied in `config/units.tsv` prior to running the test process.