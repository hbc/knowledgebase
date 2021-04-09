# Main differences in CUT&RUN to atacseq
## About trimming
### Cutruntools - excerpt
```
Short fragments are frequently encountered in CUT&RUN experiments due to the fine cutting by pA-MN enzyme.
As a result,  it is common to expect both mates of DNA fragment to overlap. 
When the fragment is shorter than the length of a read, then we can expect that adapter run-through will occur. 
It is thus critical to remove adapter sequences at the end of the reads. 
To deal with the issues caused by the alignment of short fragments, we made two important modifications to the typical adapter trimming and alignment protocol:

1. An initial trimming was first performed with Trimmomatic [9], with settings optimized to detect adapter contamination in short-read sequences. 
Trimmomatic is a template-based trimmer. However, reads containing 6 bp, or less, of adapters are not trimmed. 
Therefore, a separate tool Kseq was developed to trim up to 6-bp adapters from the 3′ end of each read that was not effectively processed by Trimmomatic. 
Note that this trimming does not affect the cut site calculation, which counts only the 5′ end of sequences. 
After trimming, a minimum read length of 25 bp was imposed, as reads smaller than this were hard to align accurately.

2. Dovetail alignment policy. 
Bowtie2 [10] aligns each mate of a pair separately and then discards any pairs that have been aligned inconsistently. 
Dovetail refers to the situation when mates extend past each other. 
In the default setting, these alignments are discarded. Dovetail is unusual but encountered in CUT&RUN experiments. 
The --dove-tail setting [10] was enabled to flag this situation as normal or “concordant” instead of elimination of such reads.
```

### MACS2 vs. SEACR
```
@houghtos
huge variability between the two programs and very sensitive to hte parameter changes.

In some cases, SEACR with IgG control was the best fit while others required tuning a numeric threshold or MACS2. 

Other variables, such as subsetting by insert size, had major effects on peaks. 

For example, in some cases it was optimal to call peaks using SEACR with IgG control on the sorted bam < 120bp.

Pipeline most closely followed cut and run tools:

    1. Trimmomatic for adapter contamination
    2. Kseq further trimming for what was missed by trimmomatic
    3. Bowtie2 align --dovetail
    4. Sort, mark, and remove duplicates (I see a following post about duplicates later)
    5. Subset bam file by < and > 120bp

Peaks were then called using MACS2 and SEACR on all resulting bam files and examined for quality control.
```

### Question regarding the insert size and the IgG samples.
```
@drpatelh
For step 5, you could just create another rendition of the file below to filter based on those insert sizes and pass that to the pipeline using --bamtools_filter_pe_config <new_config>. 

So you only want to keep read pairs that have an insert size of less than 120bp?

Also, just out of interest do you actually get anything useful from your IgG samples. 

In theory, you shouldnt be pulling down much because its an non-specific antibody and shouldnt really pull much down. 

As a result, the library complexity is quite low and most of the reads are filtered out as duplicates and so renders the data unusable as a control...

have you got a MultiQC report or something I can look at where you have analysed these? 

Be interested to see what yours look like.
```
> The answer
```
@houghtos
Honestly, unsure how changing 1-4 will result. 
Cut & run is very sensitive to parameter changes since it's sparse structure. 
I've seen pipelines using BWA etc. so it will likely get you reasonable results.

    1. Trimming: The 2-step trimming protocol was implemented due to the predominance of short 25-50bp fragments which skews normal trimming and subsequent aligning methods. 
    I haven't found much regarding trimming.

    2. Aligning: Dovetail option allows mates to be considered concordant if one extends past the other, which are normally discarded. 
    Dovetail is unusual but encountered in cut&run experiments. 
    The original author did not use dovetail but rather "--local --very-sensitive- local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700".

    3. For read pairs < 120bp - this is specific to for examining TFs in cut&run. 
    This can be made as a default, separate option, or a pipeline that includes multiple outputs. 
    As indicated in the original paper:

Separation of sequenced fragments into ≤120 bp and ≥150 bp size classes provides mapping of the local vicinity of a DNA-binding protein, but this can vary depending on the steric access to the DNA by the tethered MNase. 

Single-end sequencing is not recommended for CUT&RUN, as it sacrifices resolution and discrimination between transcription factors and neighboring nucleosomes.

I've seen most pipelines include some variation of subsetting results or including them as part of the output for further downstream analysis.

    SEACR: There have been questions surrounding SEACR (for example look at reviewer comments to original author here). 
    Here, I considered 'optimal peaks' between 20-70k. Everything outside of this range was either <5k or > 170k. 
    Majority of pipelines I see use MACS2 default. 
    SEACR is meant to address the sparse properties of cut&run (hence the sensitivity issues). 
    With MACS2 I often got very few peaks (100-1,000) with the caveat I used q-value < 0.01. 
    SEACR was incredibly imprecise with peak numbers, but I was more likely able to find optimal peaks. 
    For example, (1) sample with 100k cells had optimal peak calling as <120bp subsetted with SEACR numerical threshold 0.01 (no IgG control). 
    (2) Two samples 100k & 200k cells with H3K27ac had few peaks with MACS2 but too many with SEACR numerical threshold 0.01. Optimal peaks required IgG control.

Some thoughts:

    > Peak calling can be imprecise. 
    
    It's probably best to allow the user flexibility to re-call peaks if there are too many or too few per library. 
    SEACR is available as a web tool for this purpose here: https://seacr.fredhutch.org/. 
    MACS2 or HOMER may be better options for default peak calling with adjusted parameters.
    IgG control is used to determine the numeric threshold. 
    Signal blocks from the input files that overlap with IgG blocks are also filtered out as means to reduce false positives. 
    However, I imagine the latter step can be handled by removing duplicates and multimapping reads beforehand.

HERE: ** He says we should be removing duplicates **
```

### Part where @houghtos talks about why duplicates should not be removed
```
Wanted to add C&R may actually produce identical molecules due to base pair resolution of footprints obtained by MNase (differing from ChIP-seq). 
It's observed high yield experiments produce "PCR duplicates" that are actually biologically meaningful molecules (100k+ cells & histone modification antibodies). 
Something to keep in mind for post alignment QC in removing duplicates.
```

### Primary difference between CUT&RUN, CUT&TAG and atacseq
```
Yes the primary differences are the peak caller used (SEACR is the standard for CUT&TAG) and the spike in normalisation. 

We are building the pipeline with options for automated spike in normalisation. 

There are some other smaller differences but they are the main ones.

While SEACR is standard for CUT&TAG, it isn't for cut&run. We just have to worry about spike in normalization mainly and some of the other parameters.
```

### Discussion again on Trimming and Dovetails
> To sum up, trimming could take care the 'dovetails', but since Chris
```
original question:
@nservant
Hi all. Reading a bit about cut&run data analysis, I did not find any major differences with a standard ChIP-seq pipeline. 
The trimming may be a bit different, and some people recommand using the --dovetail  option for the bowtie2 mapping, + SEACR for peak calling. 
Is there any other differences that I missed ? 
Then, I'm a bit confused about dovetail read pairs (... when one mate alignment extends past the beginning of the other ...), 
I do not understand how this is possible? 
If someone has any reference or explanation, it would be great ! thanks !
```
```
original answer:
@christ cheshire (main dev for nf-core/cutandrun)
Hi thanks for your message - I think you are right in that the format of the analysis is quite similar as it is with most enriched fragment type genomics pipelines (ChIP-seq, ATAC-seq, CUT&RUN, CUT&TAG etc.); 
however, the things you mentioned do make quite a large difference to the final output of the pipeline.
There are two primary differences with CUT&RUN/TAG over ChIP-seq: fragment length and background noise.
In C&R, the read lengths are predominantly shorter and so normal trimming and alignment settings must be adjusted to account for this. 
For this reason, trimming is not advised in C&R and during alignment it is also advised that the dovetailed read mates are considered for alignment. 
I have included a screen shot of the bowtie2 manual below but briefly you get dovetailing when an insert is so short that the sequencer reads almost the same sequence during both read 1 and read 2. 
Usually this is a sign of a fragment you dont want or perhaps some kind of error, but in C&R this is expected and so the aligner must be told to accept these mate types.
The second main difference is the very low background noise obtained from C&R when compared to ChIP-seq. Peak callers like MACS were designed to account for large amounts of noise but they can be quite conservative because of this. 
SEACR is designed specifically for low background-noise experiments and can take advantage of this to produce more accurate peaks for C&R data.
In summary all of this really means is that less information is filtered out during processing so that you can take full advantage of the quality of C&R data over ChIP-seq.
Another minor difference is that spike-in normalisation is almost always used in C&R due to the residual ecoli genetic material that is always present in a sample from the production process of the PA-MNase; 
therefore, most pipelines include auto-spikein calibration against ecoli. 
In ChIP-seq, you must add in a specific spike-in which does not always happen.
The pipeline should be complete very soon and we hope at some point before the summer to include the MACS peak caller in the pipeline anyway, so that people can compare.
Let me know you if have further questions, I would be happy to answer them!
```
> dovetail in bowtie2
<img src="https://i.imgur.com/Oh7I152.png" width="800">

```
@Joon Yoon
Hi Chris, from the past discussions, 
I was getting the impression that we should be trimming with trimmomatic and then kseq to account for the <6bp adapters.
(https://github.com/nf-core/chipseq/issues/127)
But in this answer, you are suggesting that we shouldn't trim at all. 
So, would 'no trimming' be the most up to date suggestion?
Sorry for the naive question as I am a novice at C&R, just reading through the limited information on the web.
Thanks!
```
```@ Chris cheshire
Hey Joon, thanks for your message and my slow reply - the message is that we have to be careful with trimming. 
You still need to look for adaptor sequences and then perform limited trimming for short reads if you see these sequences in your data. 
Normal trimming strategies however may truncate your reads so I would avoid them
```

#### still has issues
```
@nservant
Many thanks @Chris Cheshire for the detailed answer. It's very clear. 
Just one last question regarding dovetail reads. Actually, I do ont understand how this can happen ! 
not only for C&R/C&T but for sequencing in general. I understand that because of fragment size, you can have a large overlap between paired reads. 
But as your insert is flanked by adapter sequences, as soon as you have one mate which overtakes the beginning of the other mate, 
it means that you start sequencing the adapter sequence, which should be trimmed. 
So except if we trim the 5' end of the reads (which is not the case here if I'm correct), how this can happen ? Sorry if my question is naive, but thanks again for your time!
```
```
@chris cheshire
hi @nservant sorry for the delayed response. 
You are right that the insert is flanked by adaptor sequences; however, these are not read as part of the read. 
In read1, sequencing primer that binds during sequencing actually binds to the R1 site, which is the inner most portion of the adaptor. 
The sequencing read then immediately moves to the insert. 
You get adaptor sequences in your reads at the END of the read if your insert is so short that your read 1 process starts reading the R2 (read2) 
adaptor sequence on the other side of the insert. 
The same applies for read 2 but in the opposite direction. Does this answer your question?
```
```
@nservant
Thanks @Chris Cheshire I fully agree. 
That's exactly why I do not understand why you can have dovetailed reads as defined in the bowtie2 manual ... 
in the sense of R1 (or R2) which goes beyond their mates ! I'm sorry if I'm not clear. I try to make a picture.
```
> the_nservant_plot
<img src="https://i.imgur.com/OBFGSI8.png" width="600">

```
@nservant
At the top, this is the definition of dovetailed reads according to bowtie2. 
At the bottom, what I would expect after trimming ... as you said the end of the reads should be adapter sequences ... 
so after trimming, we cannot have reads that extend beyond the start position of their mates. Does it make sense?
```
```
@chris cheshire
yeah thats a really good point! 
I guess the bowtie2 diagram does not talk about adapter trimming and so its only possible if you dont trim the sequences completely. 
In CUT&RUN/TAG I guess you could have dovetails depending on how aggressive your trimming is. 
As we dont trim aggressively as a rule in CUT&RUN, I suppose there is more chance of a dovetail scenario as the whole adapter sequence has not been trimmed off
which is why the setting is included specifically. 
good spot though I hadnt really thought about it properly.
```
```
@nservant
Thanks again for the nice discussion !
```
```
@Joon Yoon
I am no expert at this, and I am still interested in this discussion. 
Would there be problems with the adapters that are not removed because they are <6bp? If trimming of the reads are the issue with the dovetails, 
I was wondering if the discussion about kseq of cutruntools in this github (https://github.com/nf-core/chipseq/issues/127) is somehow relevant to this issue.
```
```
@chris cheshire
Yes perhaps, I am knee deep trying to get the pipeline compete right now, 
but this kind of thing will be critical to test when we move to the testing phase
```


## Useful Links
### First Mail Meeta received from Qian
```
Hi Meeta,
Thanks for your questions. Only a small part uses ChIP-seq pipeline. I would say that probably peak calling is the only part shared with ChIP-seq (MACS2).
There are several major considerations unique to pA-MNase CUT&RUN since it is an enzyme digestion based approach, different from sonication based ChIP-seq.
1) short fragment size. All CR are paired end and give exact fragment length. However, most fragments of interest are short (50-100bp), and when using 150bp paired end sequencing, will leave a long adapter overhang in reads. Trimming is absolutely crucial.
2) alignment is similar to most high-throughput sequencing, but we use --dove-tail setting in bowtie2 to optimize alignment.
3) identification of binding sites is based on cut frequency matrix of pA-MNase, not based on peaks alone. Different from ChIP-seq where peaks are the endpoint of an analysis, we calculate cut-frequency matrix in peak regions to further distinguish if a peak is a direct binding site or not. This allows us to give a number to each peak as to how likely it is a binding site. There is no cut-frequency matrix in ChIP-seq. This cut matrix allows users to generate beautiful motif footprinting figures.
4) the pipeline for histone CUT&RUN is different from TF CUT&RUN. We filter fragments based on fragment size depending if it is histone or TF antibody. Different fragment sizes give us great a deal of information of what regions a TF prefers to bind.
5) we further provide motif enrichment analysis to validate CUT&RUN peaks. 
So there are principled differences between CUT&RUN and ChIP-seq and a couple of nice additions in our pipeline CUT&RUNTools (by the way it is recently published in Genome Biology - users have a paper to cite now).
If there is an interest in integrating it into the Core, I could work with you to have this installed for users of your cluster environment. 
Best wishes,
Qian
```

### Lab who published CUT&RUN comes up with a pipeline - Henikoff lab FHCRC
```
CUT&RUN paper: https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4
pipeline: https://pypi.org/project/henipipe/
```
### Other papers and links
```
activemotif: https://www.activemotif.com/blog-cut-tag
cutandtools: https://bitbucket.org/qzhudfci/cutruntools/src/master/
SEACR peak caller: https://github.com/FredHutch/SEACR
hainer pipeline: https://github.com/sarahhainer/uliCUT-RUN
```

## At HBC
```
Running CUT&RUNtools on O2
<120bp will not be used. Using the all fragments peak calls.
```
### JSON file required for the CUT&RUNtools to run at the moment.
```
{
    "Rscriptbin": "/n/app/R/3.3.3/bin",
    "pythonbin": "/n/app/python/2.7.12/bin/",
    "perlbin": "/n/app/perl/5.24.0/bin",
    "javabin": "/n/app/java/jdk-1.8u112/bin",
    "trimmomaticbin": "/n/app/trimmomatic/0.36/bin",
    "trimmomaticjarfile": "trimmomatic-0.36.jar",
    "bowtie2bin": "/n/app/bowtie2/2.2.9/bin",
    "samtoolsbin": "/n/app/samtools/1.3.1/bin",
    "adapterpath": "/home/qz64/cutrun_pipeline/adapters", 
    "picardbin": "/n/app/picard/2.8.0/bin",
    "picardjarfile": "picard-2.8.0.jar",
    "macs2bin": "/n/app/macs2/2.1.1.20160309/bin",
    "macs2pythonlib": "/n/app/macs2/2.1.1.20160309/lib/python2.7/site-packages",
    "kseqbin": "/home/qz64/cutrun_pipeline", 
    "memebin": "/home/qz64/meme/bin", 
    "bedopsbin": "/n/app/bedops/2.4.30", 
    "bedtoolsbin": "/n/app/bedtools/2.27.1/bin",
    "makecutmatrixbin": "/home/jy256/.local/bin",
    "bt2idx": "/n/groups/shared_databases/bowtie2_indexes",
    "genome_sequence": "/home/qz64/chrom.hg19/hg19.fa",
    "extratoolsbin": "/home/qz64/cutrun_pipeline", 
    "extrasettings": "/home/qz64/cutrun_pipeline", 
    "input/output": {
        "fastq_directory": "/home/jy256/scratch/hbc_cutandrun_shi_yang_violetta_diff_neurons_iPSCs_hbc04185/data",
        "workdir": "/home/jy256/scratch/hbc_cutandrun_shi_yang_violetta_diff_neurons_iPSCs_hbc04185/work",
        "fastq_sequence_length": 42,
        "organism_build": "hg19"
    },
    "motif_finding": {
        "num_bp_from_summit": 150,
        "num_peaks": 5000,
        "total_peaks": 15000,
        "motif_scanning_pval": 0.0005,
        "num_motifs": 20
    },
    "cluster": {
        "email": "jyoon@hsph.harvard.edu",
        "step_alignment": {
            "queue": "short",
            "memory": 32000,
            "time_limit": "0-12:00"
        },
        "step_process_bam": {
            "queue": "short",
            "memory": 32000,
            "time_limit": "0-12:00"
        },
        "step_motif_find": {
            "queue": "short",
            "memory": 32000,
            "time_limit": "0-12:00"
        },
        "step_footprinting": {
            "queue": "short",
            "memory": 32000,
            "time_limit": "0-12:00"
        }
    }
}

```

