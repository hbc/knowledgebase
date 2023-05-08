# Leafviz - Visualize Leafcutter Results


**All scripts found in `/HBC Team Folder (1)/Resources/LeafCutter_2023`**

### Step 1 - Get annotation files

Annotation files already prepared for hg38 are in the folder listed above in a folder named new_hg38. To make annotation files for another organism run 

```bash
./gtf2leafcutter.pl -o /full/path/to/directory/annotation_directory_name/annotation_file_prefix \
/full/path/to/gencode/annotation.gtf
```

### Step 2 - Make RData files from your results

Use the `prepare_results.R` script in the linked folder. The one on the leafcutter github has not been updated.
Your leafcutter should have output `leafcutter_ds_cluster_significance_XXXXX.txt`, `leafcutter_ds_effect_sizes_XXXXX.txt`, and `XXXX_perind_numers.counts.gz`.
You will need all three of these files and the groups file you made. An example groups file for my PD comparison is below:

```bash
NCIH1568_RB1_KO_DMSO_Replicate_3-ready.bam	DMSO
NCIH1568_RB1_KO_DMSO_Replicate_2-ready.bam	DMSO
NCIH1568_RB1_KO_DMSO_Replicate_1-ready.bam	DMSO
NCIH1568_RB1_KO_PD_Replicate_1-ready.bam	PD
NCIH1568_RB1_KO_PD_Replicate_3-ready.bam	PD
NCIH1568_RB1_KO_PD_Replicate_2-ready.bam	PD
```

Below is an example to make the RData from the PD comparison. This code outputs `PD_new.RData`. 
Note that the annotation has the file path and the annotation file prefix.

```bash
./prepare_results.R -m PD_groups.txt NCIH_PD_perind_numers.counts.gz \
leafcutter_ds_cluster_significance_PD.txt leafcutter_ds_effect_sizes_PD.txt \
new_hg38/new_hg38 -o PD_new.RData
```

### Step 3 - Visualize

Once you have RDatas made for all of your contrasts it is time to actually run leafviz. 
The critical scripts to run leafviz are `run_leafviz.R `, `ui.r`, and `server.R`.  
Before you can run it you **must** change the path on line 41 of `run_leafviz.R`. 
This path needs to reflect the location of the `ui.R` and `server.R` files. 
To run leafviz with my PD data set I give

```bash
./run_leafviz.R PD_new.RData
```

This will open a new tab in your browswer with your results!

For more on what is being shown check the leafviz [documentation](http://davidaknowles.github.io/leafcutter/articles/Visualization.html)
