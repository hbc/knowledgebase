# Note on Metadata for chipseq
## Linking the inputs together with one line.
## antibody column matters!
```
samplename,description,batch,phenotype,replicate,treatment,antibody
Lib4.R1.bc.2.WTMTF2.fq,WTMTF2_1,pair1,chip,1,WT,MTF2
Lib9.R1R2.bc.19.WTMTF2.fq,WTMTF2_2,pair2,chip,2,WT,MTF2
Lib3.bc.1.WTH3K27ME3.fq,WTH3k27ME3_1,pair3,chip,1,WT,H3k27ME3
Lib9.R1R2.bc.1.WTH3K27ME3.fq,WTH3k27ME3_2,pair4,chip,2,WT,H3k27ME3
Lib2.bc.1.MKOFLAG.fq,MTF2KO_1,pair5,chip,1,MTF2KO,FLAG
Lib9.R1R2.bc.30.MKOFLAG.fq,MTF2KO_2,pair6,chip,2,MTF2KO,FLAG
Lib2.bc.2.MKOWTFLAG.fq,MTF2KO_WTRES_1,pair7,chip,1,MTF2KO_WTRES,WT-FLAG
Lib9.R1R2.bc.31.MKOWTFLAG.fq,MTF2KO_WTRES_2,pair8,chip,2,MTF2KO_WTRES,WT-FLAG
Lib2.bc.15.MKOMUTFLAG.fq,MTF2KO_MUTRES_1,pair9,chip,1,MTF2KO_MUTRES,MUT-FLAG
Lib9.R1R2.bc.32.MKOMUTFLAG.fq,MTF2KO_MUTRES_2,pair10,chip,2,MTF2KO_MUTRES,MUT-FLAG
Lib3.bc.7.EKOH3K27ME3.fq,EEDKO_1,pair11,chip,1,EEDKO,H3k27ME3
Lib3.bc.8.EKOH3K27ME3.fq,EEDKO_2,pair12,chip,2,EEDKO,H3k27ME3
Lib10.R1R2.bc.3.EKOWTRES.fq,EKOWT_1,pair13,chip,1,EKO_WT,H3k27ME3
Lib10.R1R2.bc.5.EKOWTRES.fq,EKOWT_2,pair14,chip,2,EKO_WT,H3k27ME3
Lib8.R1R2.bc.16.EKOMUTRES.fq,EKOMUT_1,pair15,chip,1,EKO_MUT,H3k27ME3
Lib10.R1R2.bc.4.EKOMUTRES.fq,EKOMUT_2,pair16,chip,2,EKO_MUT,H3k27ME3
Lib2.bc.14.INPUT.fq,input_global,pair1;pair2;pair3;pair4;pair5;pair6;pair7;pair8;pair9;pair10;pair11;pair12;pair13;pair14;pair15;pair16,input,1,WT,Input
```

## I am getting these warnings and some samples ran with broadpeak. (samples that have H3K27ME3)
> Going through the log, I found this.... as I didn't get any peaks for the MTF2 samples.
```
[2021-04-10T05:28Z] h3k27me3 specified, using broad peak settings.
[2021-04-10T05:28Z] h3k27me3 specified, using broad peak settings.
[2021-04-10T05:28Z] h3k27me3 specified, using broad peak settings.

[2021-04-10T05:28Z] mut-flag specified, but not listed as a supported antibody. Valid antibodies are {'h3k36me3', 'narrow', 'h3k4me1', 
'h2afz', 'h3ac', 'h4k20me1', 'h3k4me3', 'h3k4me2', 'h3k9ac', 'h3k79me2', 'h3k9me2', 'h3f3a', 'h3k79me3', 'h3k27me3', 'broad', 'h3k9me3', 'h3k9me1', 'h3k27ac'}. 
If you know your antibody should be called with narrow or broad peaks, supply 'narrow' or 'broad' as the antibody.

[2021-04-10T05:28Z] flag specified, but not listed as a supported antibody. Valid antibodies are {'h3k36me3', 'narrow', 'h3k4me1', 
'h2afz', 'h3ac', 'h4k20me1', 'h3k4me3', 'h3k4me2', 'h3k9ac', 'h3k79me2', 'h3k9me2', 'h3f3a', 'h3k79me3', 'h3k27me3', 'broad', 'h3k9me3', 'h3k9me1', 'h3k27ac'}. 
If you know your antibody should be called with narrow or broad peaks, supply 'narrow' or 'broad' as the antibody.

[2021-04-10T05:28Z] wt-flag specified, but not listed as a supported antibody. Valid antibodies are {'h3k36me3', 'narrow', 'h3k4me1', 
'h2afz', 'h3ac', 'h4k20me1', 'h3k4me3', 'h3k4me2', 'h3k9ac', 'h3k79me2', 'h3k9me2', 'h3f3a', 'h3k79me3', 'h3k27me3', 'broad', 'h3k9me3', 'h3k9me1', 'h3k27ac'}. 
If you know your antibody should be called with narrow or broad peaks, supply 'narrow' or 'broad' as the antibody.

[2021-04-10T05:28Z] h3k27me3 specified, using broad peak settings.
[2021-04-10T05:28Z] h3k27me3 specified, using broad peak settings.
[2021-04-10T05:28Z] h3k27me3 specified, using broad peak settings.
```
