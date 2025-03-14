3_6_25 Billingsley

### Running Baysor for spatial cell segmentation


For review

https://www.10xgenomics.com/analysis-guides/using-baysor-to-perform-xenium-cell-segmentation

https://kharchenkolab.github.io/Baysor/dev/

https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/flat-file-exports/flat-files-compare.html/

_above has some description of AtoMx transcript file output_<br><br>



https://datadryad.org/stash/dataset/doi:10.5061/dryad.37pvmcvsg#readme/

*above has some description of Baysor output data in the segmentation.csv file*<br><br>
	
https://github.com/Bassi-git/Baysor_edit<br><br>


https://vimeo.com/558564804/

*above is a very helpful and short vid*<br><br><br><br>


### Running Baysor on CosMx output.


Basic usage.



Run Baysor here on o2:

/n/app/bcbio/baysor<br><br>



First run **Baysor preview**. It requires a transcript coordinate file.  With CosMx output the file will look something like this,  BWH_20240509_WC0933_tx_file.csv.gz (unzipped and renamed here to tx_file.csv)


Baysor will estimate and construct spatially nearest neighbor clusters of transcripts into Neighborhood Composition Vectors (NCVs) which can be viewed as "pseudocells" and also estimate which transcripts are noise. It will perform unsupervised clustering of the NCVs, assigning them to clusters based on expression profile similarity. You can actually use NCVs similarly to single cell data and do things like create umaps and perform marker identification. The NCVs don't have centroids so you can't plot them spatially the same way you can cell centroids. The number of expected NCV clusters can be selected or it will use a default of 4. Without even having cells called, this can give you an idea of the different cell types present. 

The output html file will show the individual transcript locations on a spatial plot. The NCVs are colored such that NCVs with similar expression profiles are colored similarly, This can give you an idea of different cell types present and their locations as well as providing a helpful reference when later selecting segmentation parameters to help select correctly sized segments. 


In the Baysor preview command, -x, -y, and -z point to their respective column names in the transcript file. -m gives the minimum number of transcripts to be included in an NCV "pseudocell". 


```../bin/baysor/bin/baysor  preview -x x_global_px -y y_global_px -z z -g target tx_file.csv -m 10 ```<br><br><br><br>
	

	
Next, **Segmentation** can be done using a couple of different approaches. It can use prior information from previous segmentation analyses, or can run without priors. It can use an image or other information as a prior. Here I use CosMx cell identifiers as a prior (each transcript is assigned to a called cell.) You can select how much confidence to give the priors, from 0 none to 1 full. (In an attempt to give full confidence, I tried running with with --prior-segmentation-confidence 1 and it gave very odd results, so you might need to try something like 8 for high confidence.)


-p will plot the segments on an HTML file, you'll need this.

This example gives confidence of 8 to the CosMx cell calls , the cell ids are found in the "cell" column in the tx_file.csv and indicated by :cell in the Baysor run command. I also ask for the number of unsupervised clusters found to be 13 here.


```../bin/baysor/bin/baysor run -x x_global_px -y y_global_px -z z -g target tx_file.csv -p --prior-segmentation-confidence 8  :cell -m 10 --n-clusters=13```<br><br><br><br>
 
A second way to run Baysor is without a prior. This requires supplying an -s "scale" argument which corresponds to expected cell diameter. With my CosMx coordinate system (pixels) s = 5 was much much too small, s = 150 too large.  You can view the segmentation output on the output html file, and compare this to the preview html file to assess quality of segmentation.


```../bin/baysor/bin/baysor  run -x x_global_px -y y_global_px -z z -g target tx_file.csv -p -s 50 -m 10 --n-clusters=13```


Output will look like this

-rw-rw-r-- 1 jmb17 bcbio  16467552 Feb 26 13:20 segmentation_borders.html<br>
-rw-rw-r-- 1 jmb17 bcbio   2986301 Feb 26 13:18 segmentation_cell_stats.csv<br>
-rw-rw-r-- 1 jmb17 bcbio   8105239 Feb 26 13:18 segmentation_counts.loom<br>
-rw-rw-r-- 1 jmb17 bcbio 374830291 Feb 26 13:18 segmentation.csv<br>
-rw-rw-r-- 1 jmb17 bcbio   1084869 Feb 26 13:18 segmentation_diagnostics.html<br>
-rw-rw-r-- 1 jmb17 bcbio       985 Feb 26 13:20 segmentation_log.log<br>
-rw-rw-r-- 1 jmb17 bcbio       651 Feb 26 10:59 segmentation_params.dump.toml<br>
-rw-rw-r-- 1 jmb17 bcbio  11374607 Feb 26 13:18 segmentation_polygons_2d.json<br>
-rw-rw-r-- 1 jmb17 bcbio  51628559 Feb 26 13:18 segmentation_polygons_3d.json<br>

segmentation_borders.html is the segmentation image.

segmentation_cell_stats.csv is the called cell metadata which can be loaded into a Seurat Object. It includes x, y cell centroid coordinates and cell area, which is very useful for filtering out unreasonably small or large segments.

segmentation.csv are the transcript quality metrics which can be used for transcript filtering in Seurat. It also contains the x and y transcript coordinates and cell and NVC transcript assignments, which can be used to create transcript by cell (or NCV) count matrices, which can be used to create seurat objects.

segmentation_polygons_2d.json can be used for plotting the segments.<br><br><br><br>


After running baysor run, check the segmentation html image for segmentation quality, compare to the preview output. 

Then also check the range of the segmentation_cell_stats area column (column 9)

```awk -F',' 'NR>1 {if(min==""){min=max=$9}; if($9<min){min=$9}; if($9>max){max=$9}} END {print "Min:", min, "Max:", max}' segmentation_cell_stats.csv```


My CosMx data use pixels, and a conversion of .12028 px per uM, so .12028^2 px per uM^2.  Adjust the scale -s parameter if you're not returning cell areas you expect in your experiment.



















