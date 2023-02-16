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

If you are creating the simple csv samplesheet with specific oligo sequences for each sample index you may need to make edits. The I2 auto-orientation detector cannot be activated when supplying a simple csv. For example, teh NextSeq instrument needs the index 2 in reverse complement. So if you had:

```
Lane,Sample,Index,Index2
*,CD_3-GEX_03,GCGGGTAAGT,TAGCACTAAG
```

You can either: 

1. Specify the reverse complement Index2 oligo 

```
*,CD_3-GEX_03,GCGGGTAAGT,CTTAGTGCTA
```

2. Use the 10x index names, e.g. SI-TT-G6

```
*,CD_3-GEX_03,SI-TT-G6
```

You can find more information [linked here on the 10X website](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/bcl2fastq-direct#sample-sheet)

Additionally, if you find you have `AdpaterRead1` or `AdapterRead2` under "Setting" in this file (example below), you will want to remove that. For any 10x library regardless of how the demultiplexing is being done, we do not recommend trimming adapters in the Illumina -- **this will cause problems with reads in downstream analyses**.


**To create FASTQ files:**

```bash
spaceranger mkfastq --run data/11-7-2022-DeVries-10x-GEX-Visium/Files 
         --simple-csv t samplesheets/11-7-2022-DeVries-10x-GEX-Visium_Samplesheet.csv 
         --output-dir fastq/11-7-2022-DeVries-10x-GEX-Visium

```

### 2. Image files

Each slide has 4 capture areas and therefore for a single slide you should have 4 image files.

More on image types [from 10X docs here](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/image-recommendations)

* Check what type of image you have (you will need to specify in `spaceranger` with the correct flag)
* Open up the image to make sure you have the fiducial border. It's probably done for you. If there are issues with the fiducial alignment (i.e. too tall, too wide) given to you, you may need to manually align using the Loupe browser
* 

Walker dataset:
* We have four slides total
* Slides 1 &2 are combined and there are 8 image files; 4 for each slide. Presumably “Field 1 -4” are for Slide 1 and “Field 5-8” are for slide 2
    * A similar situation for Slides3&4
* Issue: We have 5 folders for the FASTQ files - 

Typically the spaceranger is run for each individual capture area; **Problem is, how do we map the correct FASTQ files to a single capture area??**


