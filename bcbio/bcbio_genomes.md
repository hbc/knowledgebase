**Homo Sapiens + Covid19**

*GRCh38_SARSCov2 * - built from ensembl
  - Assembly GRCh38, release 99.
  - Files:
    - Genomic sequence: ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    - Annotation: ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
    - Covid sequence and annotation added: https://www.ncbi.nlm.nih.gov/nuccore/MN988713.1?report=GenBank 

**Mus musculus**

*GRCm38_98* - built from ensembl
  - Assembly GRCm38, release 98.
  - Files:
    - Genomic sequence: ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz

    - Annotation: ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz

**Caenorhabditis elegans**

*WBcel235_WS272* - built from wormbase
  - Assembly WBcel235, relase WS272. Project PRJNA13578 (N2 strain)
  - Files:
      - Genomic sequence: ftp://ftp.wormbase.org/pub/wormbase/releases/WS272/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS272.genomic.fa.gz

      - Annotation: ftp://ftp.wormbase.org/pub/wormbase/releases/WS272/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS272.canonical_geneset.gtf.gz

**Drosophila melanogaster**

*DGP6* - built from Flybase
  - has a different format for annotations for non-coding genes in the gtf
  - only protein coding genes will make it into Salmon and downstream
  
*DGP6.92* - built from Ensembl info
  - will have all non-coding RNAs in Salmon and downstream results
  - shows lower gene detection rates than Flybase
 
 **Updating supported transcriptomes**
1. clone cloudbiolinux
2. update transcriptome
```bash
bcbio_python cloudbiolinux/utils/prepare_tx_gff.py --cores 8 --gtf Macaca_mulatta.Mmul_8.0.1.95.chr.gtf.gz --fasta /n/app/bcbio/biodata/genomes/Mmulatta/mmul8noscaffold/seq/mmul8noscaffold.fa Mmulatta mmul8noscaffold
```
3. upload the xz file to the bucket
```bash
aws s3 cp hg19-rnaseq-2019-02-28_75.tar.xz s3://biodata/annotation/ --grants read=uri=http://acs.amazonaws.com/groups/global/AllUsers full=emailaddress=chapmanb@50mail.com
```
4. edit cloudbiolinux ggd transcripts.yaml recipe to point to the new file uploaded on the bucket
5. edit the cloudbiolinux ggd gtf.yaml to show where you got the GTF from and what you did to it
6. test before pushing
```bash
   mkdir tmpbcbio-install
   ln -s `pwd`/cloudbiolinux tmpbcbio-install/cloudbiolinux
   log into bcbio user: sudo -su bcbio /bin/bash
   bcbio_nextgen.py upgrade --data
```
7. push changes back to cloudbiolinux

**Factual list of genomes in O2:/n/shared_db/bcbio/biodata/genomes as of 2020-03-13**
```
.
├── Ad37
│   ├── GW7619026
│   └── GW76-19026
├── Adenovirus
│   └── Ad37
├── Amexicanus
│   └── Amexicanus2
├── Amis
│   ├── ASM28112v4
│   └── ASM28112v4.a
├── Anidulans
│   └── FGSC_A4
├── Atta_cephalotes
│   └── Attacep1.0
├── bcbiotx
├── Btaurus
│   └── UMD3.1
├── Celegans
│   ├── WBcel235
│   ├── WBcel235_90
│   ├── WBcel235_raw
│   └── WBcel235_WS272
├── Dmelanogaster
│   ├── BDGP6
│   ├── BDGP6.15
│   ├── BDGP6.19
│   ├── BDGP6.92
│   ├── flybase
│   └── flybase_dmel_r6.28
├── Drerio
│   ├── Zv10
│   ├── Zv11
│   └── Zv9
├── Ecoli
│   ├── EDL933
│   ├── k12
│   ├── MB0009
│   ├── MB2409
│   ├── MB2455
│   ├── MG1655
│   ├── MG1655_v2
│   ├── MG1655_virus
│   ├── MG1655_wrong_name
│   └── NC_000913.3
├── Gallus_gallus
│   └── galgal5
├── gdc-virus
│   └── gdc-virus-hsv
├── haD37
│   └── DQ900900.1
├── Hsapiens
│   ├── GRCh37
│   ├── hg19
│   ├── hg19-ercc
│   ├── hg19-mt
│   ├── hg19-subset
│   ├── hg19-test
│   └── hg38
├── humanAd37
│   └── Ad37.hg19
├── kraken
│   ├── bcbio
│   ├── micro
│   ├── minikraken_20141208
│   ├── minimal
│   └── old_20141302
├── Lafricana
│   └── loxAfr3
├── Macaca
│   ├── Mfascicularis
│   ├── Mmul8
│   └── mmul8noscaffold
├── Mmulatta
│   ├── mmul8
│   └── mmul8noscaffold
├── Mmusculus
│   ├── cloudbiolinux
│   ├── GRCm38_90
│   ├── GRCm38_98
│   ├── greenberg-mm9
│   ├── mm10
│   └── mm9
├── Oaires
│   └── Oar_v31
├── phiX174
│   └── phix
├── Pintermedia
│   └── ASM195395v1
├── Rnorvegicus
│   └── rn6
├── Scerevisiae
│   └── sacCer3
├── Spombe
│   ├── ASM284v2.25
│   └── ASM284v2.30
├── spombe
│   └── ASM294v2
├── Sscrofa
    ├── ss11.1
    └── Sscrofa10.2
```

**How to install a custom genome in O2**
- `sudo -su bcbio /bin/bash`
- `cd /n/app/bcbio`
- https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#reference-genome-files
- https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#adding-custom-genomes

**Workflow4: Whole genome trio (50x) - hg38**

Inputs (FASTQ files) and results (BAM files, etc) of the [whole genome BWA alignment and GATK variant calling workflow](https://bcbio-nextgen.readthedocs.io/en/latest/contents/germline_variants.html#workflow4-whole-genome-trio-50x-hg38) are stored in `/n/data1/cores/bcbio/shared/NA12878-trio-eval`

**Use an updated hg38 transcriptome**
```bash
wget ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
gtf=Homo_sapiens.GRCh38.101.chr.gtf.gz
remap_url=http://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2UCSC.txt
wget --no-check-certificate -qO- $remap_url | awk '{if($1!=$2) print "s/^"$1"/"$2"/g"}' > remap.sed
gzip -cd ${gtf} | sed -f remap.sed | grep -v "*_*_alt" > hg38-remapped.gtf
```
Then pass `hg38-remapped.gtf` as the `transcriptome_gtf` option.
