# Demultiplexing SC data using SNP information

## Overview

## Methods:

* scSplit:


  **References**:
  
  * [Paper](https://doi.org/10.1186/s13059-019-1852-7) 
  
  * [Repo](https://github.com/jon-xu/scSplit)

* Demuxlet/Freemuxlet/Popscle:

  Demuxlet is the first iteration of the software. Popscle is a suite that includes an improved version of demuxlet and also freemuxlet. It is recommended 
  by the authors to use popscle.
    
  **Running it:**
  
  Installing demuxlet is not straightforward with very particular instructions. A similar situation might happen with popscle which it's not published yet.
  I recommend using Docker. The repo contains the [Dockerfile](https://github.com/statgen/popscle/blob/master/Dockerfile). You can use it to create your own 
  docker image. One available is [here](https://hub.docker.com/repository/docker/vbarrerab/popscle). This image can also be used to create a singularity container 
  on O2.
  
  _Running on O2: singularity_



singularity exec -B <local_folder_bam_files>:/bam_files,<local_folder_vcf_files>:/vcf_files,<local_folder_dsc_pileup_results>:/results
/n/app/singularity/containers/<user>/<popscle_singularity_container> popscle dsc-pileup --sam /bam_files/<bam_file> --vcf /vcf_files/<vcf_file> 
 --out /results/<pileup_file_output>

**Recommendations:**

  It is highly reccomended to reduce the number of reads and SNPs before running 
  
  

**References**:

  _Demuxlet_
  
  * [Paper - Demuxlet](https://www.nature.com/articles/nbt.4042) 
  
  * [Repo](https://github.com/statgen/demuxlet)

  _Popscle (Demuxlet/Freemuxlet)_
  
  * [Repo](https://github.com/statgen/popscle)

  _popscle helper tools_
  
  * [Repo](https://github.com/aertslab/popscle_helper_tools)
