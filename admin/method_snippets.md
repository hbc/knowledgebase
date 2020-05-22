# Method snippets

## Expression Microarrays

### Illumina arrays

Arrays were processed using the 'oligo' [Carvalho B. S., and Irizarry, R. A. (2010). A Framework for Oligonucleotide Microarray Preprocessing. Bioinformatics, 26(19):2363-7.] BioConductor package, quality-controlled with arrayQualityMetrics [Kauffmann, A., Gentleman, R.,, Huber, W. (2009) arrayQualityMetrics--a bioconductor package for quality assessment of microarray data. Bioinformatics, 25(3):415-6.]
 and RMA [Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15]
 normalized. Differentially expressed genes (as defined by a minimun 2 fold change in expression and  a  Benjamini-Hochberg FDR adjusted pvalue of less than 0.1) were identified using limma  [Smyth, GK (2005). Limma: linear models for microarray data. In:  'Bioinformatics and Computational Biology Solutions using R and  Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber  (eds), Springer, New York, pages 397-420.].”

### Affymetrix arrays

Arrays were processed using the 'oligo' [Carvalho B. S., and Irizarry, R. A. (2010). A Framework for Oligonucleotide Microarray Preprocessing. Bioinformatics, 26(19):2363-7.] BioConductor package, quality-controlled with arrayQualityMetrics [Kauffmann, A., Gentleman, R.,, Huber, W. (2009) arrayQualityMetrics--a bioconductor package for quality assessment of microarray data. Bioinformatics, 25(3):415-6.]
 and RMA [Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15]
 normalized. Differentially expressed genes (as defined by a minimun 2 fold change in expression and  a  Benjamini-Hochberg FDR adjusted pvalue of less than 0.1) were identified using limma  [Smyth, GK (2005). Limma: linear models for microarray data. In:  'Bioinformatics and Computational Biology Solutions using R and  Bioconductor'. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber  (eds), Springer, New York, pages 397-420.].”


## Variant calling
Overall, the parameters of our workflows are based on GATK best practices (https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows), contributions from bcbio community (https://github.com/bcbio/bcbio-nextgen) and our own validations (https://github.com/bcbio/bcbio_validations/).

### Read alignment
We align reads with `bwa mem` (1), (https://github.com/lh3/bwa), using samtools (2), (https://github.com/samtools/) and sambamba (3), (https://github.com/biod/sambamba) to sort bam files and mark duplicate reads.

### Quality control
We run many tools to gather QC metrics: 
- peddy (https://github.com/brentp/peddy)
- verifybamid (https://genome.sph.umich.edu/wiki/VerifyBamID)
- DKFZ bias filter(https://github.com/DKFZ-ODCF/DKFZBiasFilter)
- fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- qualimap (http://qualimap.bioinfo.cipf.es/)
- samtools (https://github.com/samtools/)
- bcftools (http://www.htslib.org/doc/bcftools.html)

We aggregate all metrics in a single QC report with multiqc (4-8), (https://multiqc.info/).

### Coverage and callable regions
We calculate coverage using mosdepth (https://github.com/brentp/mosdepth) and calculate callable regions based on real coverage and bed files provided.

### SNP and indels in germline (WES, WGS, gene panels)
We support variant calling in germline with:
- gatk4x (https://github.com/broadinstitute/gatk/)
- gatk3.8x (https://console.cloud.google.com/storage/browser/gatk-software/package-archive)
- freebayes (https://github.com/ekg/freebayes)
- samtools (https://github.com/samtools/)
- octopus (https://github.com/luntergroup/octopus).

We support:
- single sample variant calling
- batch variant calling
- population variant calling (joint genotyping).

### Structural and copy number variants in germline (WGS data)
We call structural variants with 
- manta (11), (https://github.com/Illumina/manta)
- lumpy (12), (https://github.com/arq5x/lumpy-sv)
- delly (13), (https://github.com/dellytools/delly)
- wham (14), (https://github.com/zeeev/wham)

We annotate structural variant calls with coverage information using duphold (https://github.com/brentp/duphold)

### Somatic small variants
We call somatic variants in tumor only or tumor/normal mode with:
- mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2)
- vardict (https://github.com/AstraZeneca-NGS/VarDict)
- strelka2 (https://github.com/Illumina/strelka)
- varscan2 (Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, & Wilson RK (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Research PMID: 22300766)

We support ensemble approach to combine somatic calls from several callers (https://github.com/bcbio/bcbio.variation.recall)

### Somatic copy number variants
We use 
- gatk-cnv (https://gatkforums.broadinstitute.org/gatk/discussion/9143/how-to-call-somatic-copy-number-variants-using-gatk4-cnv)
- pureCN (https://bioconductor.org/packages/release/bioc/html/PureCN.html)
- seq2c (https://github.com/AstraZeneca-NGS/Seq2C)
- titanCNA (https://github.com/gavinha/TitanCNA)

### Variant annotation (15,16)
We annotate variants with
- VEP (https://useast.ensembl.org/info/docs/tools/vep/index.html)
- snpEff (A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672].
- using vcfanno (https://github.com/brentp/vcfanno) and many annotation sources:
  - gnomad (https://gnomad.broadinstitute.org/)
  - topmed (https://bravo.sph.umich.edu/freeze5/hg38/)
  - cosmic (https://cancer.sanger.ac.uk/cosmic)
  - dbsnp  (https://www.ncbi.nlm.nih.gov/snp/)
  - dbnsfp (https://sites.google.com/site/jpopgen/dbNSFP)
  - clinvar (https://www.ncbi.nlm.nih.gov/clinvar/)

We create gemini database (17) as output (https://gemini.readthedocs.io/en/latest/). We support any internal vcf or bed based annotation (internal frequency database) via vcfanno.

1. Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. 2013 arXiv:1303.3997.
2. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9 [19505943].
3. A. Tarasov, A. J. Vilella, E. Cuppen, I. J. Nijman, and P. Prins. Sambamba: fast processing of NGS alignment formats. Bioinformatics, 2015.
4. DKFZ bias filter (https://github.com/DKFZ-ODCF/DKFZBiasFilter).
5. fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
6. qualimap (http://qualimap.bioinfo.cipf.es/).
7. bcftools (http://www.htslib.org/doc/bcftools.html)
8. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. Epub 2016 Jun 16. PMID: 27312411; PMCID: PMC5039924.
9. Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018;34(5):867‐868. doi:10.1093/bioinformatics/btx699.
10. GATK4 (https://github.com/broadinstitute/gatk/).
11. Chen, X. et al. (2016) Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32, 1220-1222. doi:10.1093/bioinformatics/btv710.
12. Layer RM, Chiang C, Quinlan AR, Hall IM. LUMPY: a probabilistic framework for structural variant discovery. Genome Biol. 2014;15(6):R84. Published 2014 Jun 26. doi:10.1186/gb-2014-15-6-r84.
13. Rausch T, Zichner T, Schlattl A, Stütz AM, Benes V, Korbel JO. DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics. 2012;28(18):i333‐i339. doi:10.1093/bioinformatics/bts378.
14. Kronenberg ZN, Osborne EJ, Cone KR, et al. Wham: Identifying Structural Variants of Biological Consequence. PLoS Comput Biol. 2015;11(12):e1004572. Published 2015 Dec 1. doi:10.1371/journal.pcbi.1004572.
15. Pedersen BS, Layer RM, Quinlan AR. Vcfanno: fast, flexible annotation of genetic variants. Genome Biol. 2016;17(1):118. Published 2016 Jun 1. doi:10.1186/s13059-016-0973-5.
16. McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biol. 2016;17(1):122. Published 2016 Jun 6. doi:10.1186/s13059-016-0974-4.
17. Paila U, Chapman BA, Kirchner R, Quinlan AR. GEMINI: integrative exploration of genetic variation and genome annotations. PLoS Comput Biol. 2013;9(7):e1003153. doi:10.1371/journal.pcbi.1003153.

### RNA-Seq

HiSeq Illumina sequencing will be performed on our behalf by ###. All samples
will be indexed so that pools can be run across ### lanes of an 8-laned Illumina
flow cell, providing an estimated 20-30 million single-end reads per sample.

All samples will be processed using an RNA-seq pipeline implemented in the
bcbio-nextgen project (https://bcbio-nextgen.readthedocs.org/en/latest/). Raw
reads will be examined for quality issues using FastQC
(http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to ensure library
generation and sequencing are suitable for further analysis.

If necessary, adapter sequences, other contaminant sequences such as polyA tails
and low quality sequences with PHRED quality scores less than five will be
trimmed from reads using atropos 
[https://github.com/jdidion/atropos; 10.5281/zenodo.596588].
Trimmed reads will be aligned to UCSC build ### of the ## # genome (###),
augmented with transcript information from Ensembl release ### using STAR
[Dobin, Alexander, Carrie A Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R Gingeras. 2013. “STAR: Ultrafast Universal RNA-Seq Aligner..” Bioinformatics (Oxford, England) 29 (1). Oxford University Press: 15–21. doi:10.1093/bioinformatics/bts635.].

Alignments will be checked for evenness of coverage, rRNA content, genomic
context of alignments (for example, alignments in known transcripts and
introns), complexity and other quality checks using a combination of FastQC,
Qualimap
[García-Alcalde, F., Okonechnikov, K., Carbonell, J., Cruz, L. M., Götz, S., Tarazona, S., et al. (2012). Qualimap: evaluating next-generation sequencing alignment data. Bioinformatics (Oxford, England), 28(20), 2678–2679. http://doi.org/10.1093/bioinformatics/bts503],
MultiQC (https://github.com/ewels/MultiQC) and custom tools.

Counts of reads aligning to known genes are generated by featureCounts
[Liao, Yang, Gordon K Smyth, and Wei Shi. 2014. “featureCounts: an Efficient General Purpose Program for Assigning Sequence Reads to Genomic Features..” Bioinformatics (Oxford, England) 30 (7). Oxford University Press: 923–30. doi:10.1093/bioinformatics/btt656.].
In parallel, Transcripts Per Million (TPM) measurements per isoform will be
generated by quasialignment using Salmon [Patro, R., Duggal, G., & Kingsford, C. (2015). Accurate, fast, and model-aware transcript expression quantification with Salmon. bioRxiv, 021592. http://doi.org/10.1101/021592].

If necessary for a project, novel transcripts will be identified via reference-guided assembly with
Stringtie [Pertea, Mihaela, Geo M. Pertea, Corina M. Antonescu, Tsung-Cheng
Chang, Joshua T. Mendell, and Steven L. Salzberg. 2015. “StringTie Enables
Improved Reconstruction of a Transcriptome from RNA-Seq Reads.” Nature
Biotechnology, February. Nature Publishing Group.
http://www.nature.com/doifinder/10.1038/nbt.3122]., with novel transcripts filtered for coding potential agreement [Wang, Liguo, Hyun Jung Park, Surendra Dasari, Shengqin Wang, Jean-Pierre Kocher, and Wei Li. 2013. “CPAT: Coding-Potential Assessment Tool Using an Alignment-Free Logistic Regression Model..” Nucleic Acids Research 41 (6). Oxford University Press: e74–e74. doi: 10.1093/nar/gkt006.] with known genes to reduce false positive assemblies.

Differential expression at the gene level will be called with DESeq2 [Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2..” Genome Biology 15 (12): 550. doi:10.1186/PREACCEPT-8897612761307401.], preferring to use counts per gene estimated from the Salmon quasialignments by tximport [Soneson, C., Love, M. I., & Robinson, M. D. (2016). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 4. http://doi.org/10.12688/f1000research.7563.1]. Quantitating at the isoform level has been shown to produce more accurate results at the gene level.

Isoform differential expression will be called using the TPM quasialignments from Salmon converted to HDF5 format for use with Sleuth (https://github.com/pachterlab/sleuth).
If necessary for a project, splicing event level calls will be generated using a combination of DEXSeq [Anders, S., Reyes, A., & Huber, W. (2012). Detecting differential usage of exons from RNA-seq data. Genome Research, 22(10), 2008–2017. http://doi.org/10.1101/gr.133744.111] and MISO [Analysis and design of RNA sequencing experiments for identifying isoform regulation. (2010). Analysis and design of RNA sequencing experiments for identifying isoform regulation].

Lists of differentially expressed genes will be examined for gene ontology (GO)
and KEGG term enrichment with clusterProfiler [Yu G, Wang L, Han Y and He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), pp. 284-287.]. Functional redundancy in enriched
GO terms will be reduced with GOSemSim [Yu G, Li F, Qin Y, Bo X, Wu Y and Wang S (2010). “GOSemSim: an R package for measuring semantic similarity among GO terms and gene products.” Bioinformatics, 26(7), pp. 976-978.] In addition, a cut-off-free gene set
enrichment analysis (GSEA) will be performed using clusterProfiler and weighted fold change calculations from DESeq2.


### mRNA

RNA quality will be performed using the Agilent Bioanalyzer (Agilent, Santa Clara, CA). Transcriptional profiling using the ### array will be performed at ###. All data analysis will be performed using the BioConductor framework [1]. Array quality will be assessed using the arrayQualityMetrics package [2], batch adjusted [3], normalized and summarized using Robust Multiarray Averaging (RMA) [4] and tested for genes with significant differential expression using limma [5].

[1]: Reimers, Mark, and Vincent J Carey. “Bioconductor: an Open Source Framework for Bioinformatics and Computational Biology.” Methods in Enzymology 411: 119–134. doi:10.1016/S0076-6879(06)11008-3.
[2]: Kauffmann, Audrey, Robert Gentleman, and Wolfgang Huber. “arrayQualityMetrics--a Bioconductor Package for Quality Assessment of Microarray Data.” Bioinformatics (Oxford, England) 25, no. 3: 415–416. doi:10.1093/bioinformatics/btn647.
[3]: Leek, Jeffrey T, Robert B Scharpf, Héctor Corrada Bravo, David Simcha, Benjamin Langmead, W Evan Johnson, Donald Geman, Keith Baggerly, and Rafael A Irizarry. “Tackling the Widespread and Critical Impact of Batch Effects in High-Throughput Data.” Nature Reviews Genetics 11, no. 10: 733–739. doi:10.1038/nrg2825.
[4]: Wilson, Claire L, and Crispin J Miller. “Simpleaffy: a BioConductor Package for Affymetrix Quality Control and Data Analysis..” Bioinformatics (Oxford, England) 21, no. 18: 3683–3685. doi:10.1093/bioinformatics/bti605.
[5]: Smyth, G. K. (2005). Limma: linear models for microarray data. In: Bioinformatics and Computational Biology Solutions using R and Bioconductor, R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W. Huber (eds.), Springer, New York, pages 397-420.


### small RNA-seq

All samples will be processed using an small RNA-seq pipeline implemented in the
bcbio-nextgen project (https://bcbio-nextgen.readthedocs.org/en/latest/). Raw
reads will be examined for quality issues using FastQC
(http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to ensure library
generation and sequencing are suitable for further analysis.

Adapter sequences, other contaminant sequences such as polyA tails and low quality sequences with PHRED quality scores less than five will be trimmed from reads using atropos (http://dx.doi.org/10.14806/ej.17.1.200).

Trimmed reads will be aligned to miRBase v21(Enright et al., n.d.) to the specific species with seqbuster (Lorena Pantano, Estivill, and Martí 2010). As well, they will be aligned to ## # genome (version ##) using STAR (Dobin et al. 2013). The aligned genomes will be used with seqcluster (Lorena Pantano, Estivill, and Martí 2011) to characterize the whole small RNA transcriptome and classify reads into rRNA, miRNA, repeats, genes, tRNAs and others from USCC annotation (Mangan et al. 2014). Finally, aligned reads were used with miRDeep2 (Mackowiak 2011), an algorithm that assesses the fit of sequenced RNAs to a biological model of miRNA generation and correct folding.
Alignments will be checked for evenness of coverage, rRNA content, genomic context of alignments (for example, alignments in known transcripts and introns), complexity and other quality checks using a combination of FastQC,  MultiQC  (Ewels et al. 2016) and custom code inside bcbio-nextgen pipeline.

Data will be loaded into R with bcbioSmallRna R package (https://github.com/lpantano/bcbioSmallRna) and isomiRs BioC package (Ramos et al. 2017; L. Pantano, Escaramis, and Argyropoulos 2016).

Differential expression at the gene level will be called with DESeq2 [Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2..” Genome Biology 15 (12): 550. doi:10.1186/PREACCEPT-8897612761307401.], after loading the data with bcbioSmallRna R package [https://github.com/lpantano/bcbioSmallRna] and isomiRs BioC package [Lorena Pantano, . (2016, January 31). isomiRs v0.99.12 (Version v0.99.12). Zenodo. http://doi.org/10.5281/zenodo.45382].


### miRNA/mRNA integration

To comprehensively integrate the interaction of differentially expressed genes and miRNAs we will use miRTrail [ Laczny, Cedric, Petra Leidinger, Jan Haas, Nicole Ludwig, Christina Backes, Andreas Gerasch, Michael Kaufmann, et al. “miRTrail--a Comprehensive Webserver for Analyzing Gene and miRNA Patterns to Enhance the Understanding of Regulatory Mechanisms in Diseases..” BMC Bioinformatics 13 (2012): 36. doi:10.1186/1471-2105-13-36.] which integrates information on 20.000 genes and almost 1.000 miRNAs. A secondary screen will be performed using MMIA [Nam, Seungyoon, Meng Li, Kwangmin Choi, Curtis Balch, Sun Kim, and Kenneth P Nephew. “MicroRNA and mRNA Integrated Analysis (MMIA): a Web Tool for Examining Biological Functions of microRNA Expression.” Nucleic Acids Research 37, no. Web Server issue (July 1, 2009): W356–62. doi:10.1093/nar/gkp294.], an algorithm incorporating three miRNA prediction algorithms (TargetScan, PITA and PicTar) into the mRNA expression analysis.


### ChIP-Seq

ChIP-seq data quality will be evaluated using FASTQC [6], and if required filtering and trimming of reads will be performed with Atropos [7]. High quality reads will be mapped to the current mouse genome build using Bowtie2 [8]. Multi-mapping reads will be excluded. We will use MACS2 [9] to call peaks on unique reads only and assess peak quality using ChIPQC [10]. Depending on signal strength and the width of the peaks, EPIC2 [11] may also be used to verify broad peak calls. Likely peak artifacts will be filtered out using the ENCODE blacklist [12]. Data will be visualized using IGV [13]. DeepTools [14] will be used to assess coverage and the reproducibility of peaks across replicates. Consensus peak sets will be generated by expanding each summit 250 bases and choosing among overlapping peaks the peak with the highest score [15]. ChIPseeker [16] will be used to annotate ChIP-seq peaks to their closest transcriptional start site and perform functional enrichment analysis. Binding sites will be evaluated using the MEME [17] suite of tools for motif discovery and motif enrichment. 

[6]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[7]: Didion JP, Martin M, Collins FS. Atropos: specific, sensitive, and speedy trimming of sequencing reads. PeerJ. 2017 Aug 30;5:e3720. doi: 10.7717/peerj.3720. PMID: 28875074; PMCID: PMC5581536.

[8]: Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

[9]: Zhang, Yong, Tao Liu, Clifford A Meyer, Jérôme Eeckhoute, David S Johnson, Bradley E Bernstein, Chad Nussbaum, et al. “Model-Based Analysis of ChIP-Seq (MACS).” Genome Biology 9, no. 9: R137. doi:10.1186/gb-2008-9-9-r137.

[10]: Carroll TS, Liang Z, Salama R, Stark R, de Santiago I. Impact of artifact removal on ChIP quality metrics in ChIP-seq and ChIP-exo data. Front Genet. 2014;5(APR):1–11.

[11]: Stovner EB, Sætrom P. epic2 efficiently finds diffuse domains in ChIP-seq data. Bioinformatics. 2019 Nov 1;35(21):4392-4393. doi:
10.1093/bioinformatics/btz232. PubMed PMID: 30923821.

[12]: Amemiya, H.M., Kundaje, A. & Boyle, A.P. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. Sci Rep 9, 9354 (2019). https://doi.org/10.1038/s41598-019-45839-z.

[13]: Robinson JT, Thorvaldsdóttir H, Winckler W, Guttman M, Lander ES, Getz G, Mesirov JP. Integrative genomics viewer. Nat Biotechnol. 2011 Jan;29(1):24-6. doi: 10.1038/nbt.1754. PMID: 21221095; PMCID: PMC3346182.

[14]: Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike Dündar, and Thomas Manke. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.

[15]: Neph S, Kuehn MS, Reynolds AP, et al. BEDOPS: high-performance genomic feature operations. Bioinformatics. 2012;28(14):1919‐1920. doi:10.1093/bioinformatics/bts277.

[16]: Yu G, Wang L, He Q (2015). “ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.” Bioinformatics, 31(14), 2382-2383. doi: 10.1093/bioinformatics/btv145.
[17]: Bailey, Timothy L, Mikael Bodén, Fabian A Buske, Martin Frith, Charles E Grant, Luca Clementi, Jingyuan Ren, Wilfred W Li, and William S Noble. “MEME SUITE: Tools for Motif Discovery and Searching..” Nucleic Acids Research 37, no. Web Server issue: W202–8. doi:10.1093/nar/gkp335.


### TF

Likely transcription factors (TF) associated with differentially expressed genes will be identified using oPOSSUM-3[1], a web-accessible software system for identification of over-represented transcription factor binding sites (TFBS) and TFBS families in DNA sequences of co-expressed genes.

[1]: Kwon, Andrew T, David J Arenillas, Rebecca Worsley Hunt, and Wyeth W Wasserman. “oPOSSUM-3: Advanced Analysis of Regulatory Motif Over-Representation Across Genes or ChIP-Seq Datasets..” G3 (Bethesda, Md.) 2, no. 9 (September 2012): 987–1002. doi:10.1534/g3.112.003202.


### RRBS

HiSeq Illumina sequencing will be performed on our behalf by ###. All samples will be indexed so that ### RRBS pools can be run on a single lane of an 8-laned Illumina flow cell, providing an estimated ### million reads per sample.

Illumina sequence quality will be surveyed using FastQC [Andrews, S.: FastQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/] to ensure library generation and sequencing are suitable for further analysis. Particular attention will be paid to per base sequence quality, overrepresented sequences (an indicator of adapter contamination) and the per base cytosine content (an indicator of conversion efficiency). Using a bioinformatics pipeline constructed in bpipe [Sadedin, S. P., Pope, B., & Oshlack, A. (2012). Bpipe: a tool for running and managing bioinformatics pipelines. Bioinformatics (Oxford, England), 28(11), 1525–1526. doi:10.1093/bioinformatics/bts167], sequences will then be trimmed, aligned and their cytosine methylation status determined. Adapters will be removed and low quality bases (<25–30 Phred quality scores) adaptively trimmed from reads with Trim Galore [Krueger F: Trim Galore!. http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/]. In addition, bases containing a cytosine artificially introduced in the end-repair step during the MspI-RRBS library preparation will be removed. Sequence quality will be re-assessed with FastQC and trimming parameters adjusted if necessary. Trimmed reads will be aligned to the appropriate reference genome with Bismark [Krueger, F., & Andrews, S. R. (2011). Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics (Oxford, England), 27(11), 1571–1572. doi:10.1093/bioinformatics/btr167], a ‘three-letter’ bisulfite aligner based on the gapped aligner Bowtie2 [Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. doi:10.1038/nmeth.1923]. Aligned reads will be automatically trimmed with BSeQC to remove nucleotides subject to further technical biases that can result in inaccurate methylation estimation, such as higher bisulfite conversion failure at the 5’ end of reads [Lin, X., Lin, X., Sun, D., Sun, D., Rodriguez, B., Rodriguez, B., et al. (2013). BSeQC: quality control of bisulfite sequencing experiments. Bioinformatics, 29(24), 3227–3229. doi:10.1093/bioinformatics/btt548]. Cytosine methylation will be determined with Bis-SNP [Liu, Y., Siegmund, K. D., Laird, P. W., & Berman, B. P. (2012). Bis-SNP: Combined DNA methylation and SNP calling for Bisulfite-seq data. Genome Biology, 13(7), R61. doi:10.1186/gb–2012–13–7-r61
], a GATK [McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., et al. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. doi:10.1101/gr.107524.110] based framework capable of simultaneous genotyping and DNA methylation estimation; methylation at non-CpG cytosines will be used to assess conversion rates. DNA methylation will be validated by bisulfite PCR or pyrosequencing on a subset of CpGs.

### inDrop RNA-seq
Sequencing via the inDrop method (XXX link for inDrop) will be performed by the
Single Cell Core at HMS. Prior to library preparation, cells will be rigorously
quality controlled for viability. A total of XXX cells per sample will be
captured, encapsulated in hydrogel and libraries will be prepared using a
modified version of the protocol outlined previously [Klein AM, Mazutis L,
Akartuna I, Tallapragada N, Veres A, Li V, Peshkin L, Weitz DA, Kirschner MW
(2015). Droplet Barcoding for Single-Cell Transcriptomics Applied to Embryonic
Stem Cells.]. Each cell will be marked by a unique cellular barcode and all
transcript fragments for each cell will have a universal molecular identifier
(UMI) attached. Sequencing of the inDrop libraries will be performed at
(XXX sequencing core), with a target of (XXX) reads per sample.

Post sequencing, reads will be assigned to each cell via identifying the
cellular barcodes, and the UMI will be extracted from each read
[Power Analysis of Single Cell RNA‐Sequencing Experiments View ORCID ProfileValentine Svensson, Kedar N Natarajan, Lam-Ha Ly, Ricardo J Miragaia, Charlotte Labalette, Iain C Macaulay, Ana Cvejic, Sarah A Teichmann
doi: http://dx.doi.org/10.1101/073692]. The
distribution of reads per cell will be used to identify a cutoff for total reads
sequenced that is the hallmark of high quality RNA from a cell. Reads from cells
passing the quality filter will be aligned to (XXX genome/build) using Rapmap [RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes. (2016). RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes., 32(12), i192–i200. http://doi.org/10.1093/bioinformatics/btw277]
and counts of reads per transcript per unique UMI will be generated for each
cell. Reads will also be aligned to (XXX genome/build) using kallisto  in
single-cell mode to generate transcript compatibility counts (TCC) [Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts. (2016). Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts., 17(1), 112. http://doi.org/10.1186/s13059-016-0970-8]for
each cell.

Heterogeneity analysis of the UMI disambiguated counts per gene per cell will be
performed using the Seurat R package [Spatial reconstruction of single-cell gene expression data. (2015). Spatial reconstruction of single-cell gene expression data., 33(5), 495–502. http://doi.org/10.1038/nbt.3192]. This includes normalization and transformation of the raw gene counts per cell to account for differences in sequencing depth, identification of high variance genes, regression of sources of unwanted variation (e.g. cell cycle phrase), identification of primary sources of heterogeneity using PCA analysis, and clustering of cells based on significant PCs (metagenes). 

Non-linear dimensional reduction using tSNE and UMAP [Etienne Becht, Charles-Antoine Dutertre, Immanuel W.H. Kwok, Lai Guan Ng, Florent Ginhoux, Evan W Newell (2018). Evaluation of UMAP as an alternative to t-SNE for single-cell data. bioRxiv 298430; doi: https://doi.org/10.1101/298430.] will be performed to visualize and explore data. Cluster quality control will be performed to assess possible cluster artifacts (variance correlated with UMI counts, mitochondrial ratio, cell cycle batch effects and any other principle component biases). After assigning clusters using sets of marker genes, differential expression analysis will be performed using zinbwave [Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert J (2018). “A general and flexible method for signal extraction from single-cell RNA-seq data.” Nature Communications, 9, 284. https://doi.org/10.1038/s41467-017-02554-5] and edgeR [Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140].

### 10X single cell analysis

We used Cell Ranger 3.0.2 (10X Genomics) to process the raw sequencing data. This pipeline converted Illumina basecall files to fastq format, aligned sequencing reads to the mm10 transcriptome using the STAR aligner (Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PubMed PMID: 23104886; PubMed Central PMCID: PMC3530905.), and quantified the expression of transcripts in each cell using Chromium barcodes. 
 
We carried out analyses of processed scRNASeq data in R version [#.#.#] using the Seurat suite version [#.#.#] (Butler A, Hoffman P, Smibert P, Papalexi E, Satija R. Integrating single-cell  transcriptomic data across different conditions, technologies, and species. Nat Biotechnol. 2018 Jun;36(5):411-420. doi: 10.1038/nbt.4096. Epub 2018 Apr 2. PubMed PMID: 29608179.; Macosko EZ, Basu A, Satija R, Nemesh J, Shekhar K, Goldman M, Tirosh I, Bialas AR, Kamitaki N, Martersteck EM, Trombetta JJ, Weitz DA, Sanes JR, Shalek AK, Regev A, McCarroll SA. Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets. Cell. 2015 May 21;161(5):1202-1214. doi: 10.1016/j.cell.2015.05.002. PubMed PMID: 26000488; PubMed Central PMCID: PMC4481139) and tidyverse packages (Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. R package version 1.2.1. https://CRAN.R-project.org/package=tidyverse). We obtained data from [#] cells from [#] samples that passed quality control steps implemented in Cell Ranger. 
  
As a further quality-control measure, we filtered out cells meeting any of the following criteria: <### or >#### unique genes expressed, <#### or >#### UMIs, or >#% of reads mapping to mitochondria. These steps removed an additional #### cells, resulting in a final dataset of ### cells. We quantified gene expression across ###  genes, with ### cells expression > # genes. 
 
To explore transcriptional heterogeneity and to undertake initial cell clustering, we reduced dimensionality using PCA. For the overall dataset, we selected ## PCs that explained more variability than expected by chance using a permutation-based test as implemented in Seurat (Macosko et al., 2015). We used PC loadings as input for a graph-based approach to cluster cells by cell type (Villani et al., 2017) and as input for t-distributed stochastic neighbor embedding (t-SNE; L van der Maaten, G Hinton (2008). “Visualizating data using t-SNE.” The Journal of Machine Learning Research 9 (2579-2605), 85) for reduction to two dimensions for visualization purposes.  To maintain a standard procedure for clustering, we used a value of #.# for the resolution. This resulted in ## clusters.
 
Specific gene markers for each cluster were identified with the FindAllMarkers() function using the Wilcoxon rank sum test. Cluster quality was assessed for possible cluster artifacts (variance correlated with UMI counts, mitochondrial ratio, batch effects, and any other principle component biases). Clusters were assigned based on specific gene markers as visualized in Spring (Weinreb C, Wolock S, Klein AM. SPRING: a kinetic interface for visualizing high dimensional single-cell expression data. Bioinformatics. 2018 Apr 1;34(7):1246-1248. doi: 10.1093/bioinformatics/btx792. Epub 2017 Dec 7. PubMed PMID: 29228172; PubMed Central PMCID: PMC6030950.). Differential expression at the gene level between clusters was performed with DESeq2 1.20 (Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.] (https://bioconductor.org/packages/DESeq2/). 

### DGE
Reads were processed to counts through the bcbio-nextgen single cell/DGE RNA-seq analysis pipeline (1). A brief description follows: The well barcode and UMIs were identified for all reads and all reads not within one edit distance of a known well barcode were discarded. Each surviving read was quasialigned to the transcriptome (GRCh38)using RapMap (2) Reads per well were counted using umis (3), discarding duplicated UMIs, weighting multimapped reads by the number of transcripts they aligned to and collapsing counts to genes by adding all counts for each transcript of a gene. The R package DESeq2 (4) will be used for differential expression analysis.

1. https://bcbio-nextgen.readthedocs.io/en/latest/contents/3prime_dge.html.
2. Srivastava, Avi, Hirak Sarkar, Nitish Gupta, and Rob Patro. 2016. “RapMap: A Rapid, Sensitive and Accurate Tool for Mapping RNA-Seq Reads to Transcriptomes.” Bioinformatics 32 (12): i192–200.
3. Svensson, Valentine, Kedar Nath Natarajan, Lam-Ha Ly, Ricardo J. Miragaia, Charlotte Labalette, Iain C. Macaulay, Ana Cvejic, and Sarah A. Teichmann. 2017. “Power Analysis of Single-Cell RNA-Sequencing Experiments.” Nature Methods 14 (4): 381–87.
4. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi:10.1186/s13059-014-0550-8.

### ATAC-seq
Reads were aligned with bwa[1] version 0.7.17-r1188 against mm10. Alignments to mitochondria and non-assembled chromosomes, duplicate alignments and multimapping alignments were removed. Alignment files were split into nucleosome (NF) free, mono/di/tri nucleosome free fractions and peaks were called individually with MACS2[2] version 2.2.6 with the parameters `--nomodel -f BAMPE -g 2730871774 --bdg --nolambda`. Consensus peaks were called on NF regions by expanding each summit 250 bases and choosing among overlapping peaks the peak with the highest score (see https://bedops.readthedocs.io/en/latest/content/usage-examples/master-list.html). Counts per peak per sample were called using featureCounts with the consensus peak file using the NF fraction BAM file. Differential peaks were called with DESeq2[3], fitting a model of the form `~0 + sg` where sg is a combined sex-genotype factor. Differential expression of the various comparisons were made with contrasts of the sex-genotype factor, using a BH adjusted p-value cutoff of 0.1. Peaks were annotated for genomic context using the ChIPseeker[4] package.

[1]: Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].  whole BWA package)
[2]: Zhang et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol (2008) vol. 9 (9) pp. R137
[3]: Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.
[4]: Yu G, Wang L, He Q (2015). “ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.” Bioinformatics, 31(14), 2382-2383. doi: 10.1093/bioinformatics/btv145


## Mass spectrometry based

### Metabolomics
Our analysis strategy for the metabolomic data generally follows the methods implemented in MetaboAnalyst [1], but has the advantage of versatility for customized experimental design. Metabolites with constant values across samples and with more than 80% values missing will be removed before further analysis. Remaining missing values will be imputed by k-means imputation using the impute package in R [2]. A median normalization will be performed for each sample and the resulting data log-transformed. The data will then be examined via Principal Component Analysis (PCA) and hierarchically clustering to identify outliers and potential confounding factors. Finally, differentially abundant metabolites will be identified using Limma [3], which allows for complex customized comparisons.
 
1. Chong, J., Soufan, O., Li, C., Caraus, I., Li, S., Bourque, G., Wishart, D.S. and Xia, J. (2018) MetaboAnalyst 4.0: towards more transparent and integrative metabolomics analysis. Nucl. Acids Res. 46, W486-494. 
2. Trevor Hastie, Robert Tibshirani, Balasubramanian Narasimhan and Gilbert Chu (2019). impute: impute: Imputation for microarray  data. R package version 1.60.0.
3. Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.

