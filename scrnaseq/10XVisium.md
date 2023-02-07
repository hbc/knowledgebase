## Analysis of 10X Visium data

> Download and install spaceranger: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/installation

Helpful resource: https://lmweber.org/OSTA-book/

### Analysis software/packages 
* [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)
* [Spatial transcriptomics with Seurat](https://yu-tong-wang.github.io/talk/sc_st_data_analysis_R.html)
* [Spatial single-cell quantification with alevin-fry](https://combine-lab.github.io/alevin-fry-tutorials/2021/af-spatial/)

## Start with the raw data: preprocessing


### 1. BCL to FASTQ

`spaceranger mkfastq` can be used here. Input is the flow cell directory.

Note that **if your SampleSheet is formatted for BCL Convert**, which is Illumina's new demultiplexing software that is soon going to replace bcl2fastq, you will get an error.
 
You will need to change the formatting slightly.
 
https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/mkfastq#simple_csv

Additionally, if you find you have `AdpaterRead1` or `AdapterRead2` under "Setting" in this file (example below), you will want to remove that. For any 10x library regardless of how the demultiplexing is being done, we do not recommend trimming adapters in the Illumina -- **this will cause problems with reads in downstream analyses**.

```
[BCLConvert_Settings]		
SoftwareVersion	3.8.4	
NoLaneSplitting	TRUE	
AdapterRead2	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT	
FastqCompressionFormat	gzip	
BarcodeMismatchesIndex1	0	
BarcodeMismatchesIndex2	0	![image](https://user-images.githubusercontent.com/6545708/217354159-731e67d3-d9a9-4ee3-8dbf-8aa12fe4e6d7.png)
```

**To create FASTQ files:**

```bash
spaceranger mkfastq --run data/11-7-2022-DeVries-10x-GEX-Visium/Files 
         --samplesheet samplesheets/11-7-2022-DeVries-10x-GEX-Visium_Samplesheet.csv 
         --output-dir fastq/11-7-2022-DeVries-10x-GEX-Visium

```


