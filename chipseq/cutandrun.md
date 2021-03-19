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

## Useful Links
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
