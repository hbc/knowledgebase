## Analysis of 10X Visium data

> Download and install spaceranger: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/installation

Helpful resource: https://lmweber.org/OSTA-book/

### Analysis software/packages 
* [scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html)
* [Spatial transcriptomics with Seurat](https://yu-tong-wang.github.io/talk/sc_st_data_analysis_R.html)
* [Spatial single-cell quantification with alevin-fry](https://combine-lab.github.io/alevin-fry-tutorials/2021/af-spatial/)


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


### 3. Counting expression data

The next step is to quantify expression for each capture area. To do this we will use [`spaceranger count`](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/tutorials/count-ff-tutorial). This command will need to be run for each capture area. Below is the command for a single capture area (in this case Slide 1, capture area A0. You may find your files have not been named with A-D, so map them accordingly.

A few things to note if you have samples that were run on multiple flow cells:

* include the `--sample` argument to specify the samplename which corresponds to the capture area
* for the `--fastqs` you can add multiple paths to the different flow cell folders and separate them by a comma

```bash

spaceranger count --id="CD_Visium_01" \
                   --sample=CD_Visium_01 \
                   --description="Slide1_CaptureArea1" \
                   --transcriptome=refdata-gex-mm10-2020-A \
                   --fastqs=mkfastq/11-10-2022-DeVries-10x-GEX-Visium/AAAW33YHV/,mkfastq/11-4-2022-Devries-10x-GEX-Visium/AAAW3N3HV/,mkfastq/11-7-2022-DeVries-10x-GEX-Visium/AAAW352HV/,mkfastq/11-8-2022-DeVries-10x-GEX-Visium/AAAW3FCHV/,mkfastq/11-9-2022-DeVries-10x-GEX-Visium/AAAW3F3HV/ \
                   --image=images/100622_Walker_Slides1_and_2/V11S14-092_20221006_01_Field1.tif\
                   --slide=V11S14-092 \
                   --area=A1 \
                   --localcores=6 \
                   --localmem=20

```
