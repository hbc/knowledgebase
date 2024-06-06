#/bin/bash
# Written by Will Gammerdinger at HSPH on September 15th, 2022.

# Assign input and output files to variables
metadata_file=/n/data1/cores/bcbio/PIs/andrew_lassar/hbc_lassar_mouse_TF_profiling_FoxC_atac_cutnrun_rnaseq_hbc04930/cutnrun/meta/FOXC.csv
multiqc_file=/n/data1/cores/bcbio/PIs/andrew_lassar/hbc_lassar_mouse_TF_profiling_FoxC_atac_cutnrun_rnaseq_hbc04930/cutnrun/final/2024-02-01_FOXC/multiqc/multiqc_data/multiqc_general_stats.txt
sample_directory=/n/data1/cores/bcbio/PIs/andrew_lassar/hbc_lassar_mouse_TF_profiling_FoxC_atac_cutnrun_rnaseq_hbc04930/cutnrun/final/
sample_prefix=X
input_antibody_label=input
output_file=/n/data1/cores/bcbio/PIs/andrew_lassar/hbc_lassar_mouse_TF_profiling_FoxC_atac_cutnrun_rnaseq_hbc04930/cutnrun/final/summarized_report.txt

# Print Header line
echo -e "sample\tantibody\ttreatment\treads\tmapped_reads\tmapping_percent\tpeaks\tRiP_percent\tPBC1\tPBC2\tBottlenecking\tNRF\tComplexity\tGC_percent" > $output_file;

# Determine the column number for the columns on interest
antibody_column=`awk -F ',' 'NR==1{for (i=1; i<=NF; i++) { if ($i == "antibody") { print i } }}' $metadata_file`
phenotype_column=`awk -F ',' 'NR==1{for (i=1; i<=NF; i++) { if ($i == "phenotype") { print i } }}' $metadata_file`
treatment_column=`awk -F ',' 'NR==1{for (i=1; i<=NF; i++) { if ($i == "treatment") { print i } }}' $metadata_file`
reads_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-Total_reads") { print i } }}' $multiqc_file`
mapped_reads_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-Mapped_reads") { print i } }}' $multiqc_file`
mapping_percent_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "Samtools_mqc-generalstats-samtools-reads_mapped_percent") { print i } }}' $multiqc_file`
RiP_percent_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-RiP_pct") { print i } }}' $multiqc_file`
PBC1_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-PBC1") { print i } }}' $multiqc_file`
PBC2_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-PBC2") { print i } }}' $multiqc_file`
Bottlenecking_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-bottlenecking") { print i } }}' $multiqc_file`
NRF_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-NRF") { print i } }}' $multiqc_file`
Complexity_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "bcbio_mqc-generalstats-bcbio-complexity") { print i } }}' $multiqc_file`
GC_percent_column=`awk 'NR==1{for (i=1; i<=NF; i++) { if ($i == "FastQC_mqc-generalstats-fastqc-percent_gc") { print i } }}' $multiqc_file`

# For each sample gather the various statistics and print to a the summarized report
for i in ${sample_directory}${sample_prefix}*; do 
    sample=`basename $i`
    antibody=`grep $sample $metadata_file | awk -F ',' -v antibody_column=$antibody_column '{print $antibody_column}'`
    sample_type=`grep $sample $metadata_file | awk -F ',' -v phenotype_column=$phenotype_column '{print $phenotype_column}'`
    treatment=`grep $sample $metadata_file | awk -F ',' -v treatment_column=$treatment_column '{print $treatment_column}'`
    reads=`grep $sample $multiqc_file | awk -F'\t' -v reads_column=$reads_column '{print $reads_column}' | sed 's/.0$//g'`
    mapped_reads=`grep $sample $multiqc_file | awk -F'\t' -v mapped_reads_column=$mapped_reads_column '{print $mapped_reads_column}' | sed 's/.0$//g'`
    mapping_percent=`grep $sample $multiqc_file  | awk -F'\t' -v mapping_percent_column=$mapping_percent_column '{print $mapping_percent_column}' | head -c 5`
    if [[ $antibody == $input_antibody_label ]]; then
        peaks="NA"
    else
        peaks=`wc -l ${sample_directory}${sample}/macs2/${sample}_peaks.* | awk 'NR==1{print $1}'`
    fi
    RiP_percent=`grep $sample $multiqc_file | awk -F'\t' -v RiP_percent_column=$RiP_percent_column '{print $RiP_percent_column}'`
    PBC1=`grep $sample $multiqc_file | awk -F'\t' -v PBC1_column=$PBC1_column '{print $PBC1_column}' | head -c 5`
    PBC2=`grep $sample $multiqc_file | awk -F'\t' -v PBC2_column=$PBC2_column '{print $PBC2_column}' | head -c 5`
    Bottlenecking=`grep $sample $multiqc_file | awk -F'\t' -v Bottlenecking_column=$Bottlenecking_column '{print $Bottlenecking_column}'`
    NRF=`grep $sample $multiqc_file | awk -F'\t' -v NRF_column=$NRF_column '{print $NRF_column}' | head -c 5`
    Complexity=`grep $sample $multiqc_file | awk -F'\t' -v Complexity_column=$Complexity_column '{print $Complexity_column}'`
    GC_percent=`grep $sample $multiqc_file | awk -F'\t' -v GC_percent_column=$GC_percent_column '{print $GC_percent_column}'`
    echo -e  "$sample\t$antibody\t$treatment\t$reads\t$mapped_reads\t$mapping_percent\t$peaks\t$RiP_percent\t$PBC1\t$PBC2\t$Bottlenecking\t$NRF\t$Complexity\t$GC_percent" >> $output_file;
done

echo -e "The summarized report has been created and can be found here:\n\t$output_file"
