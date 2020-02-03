**Admin**: [Getting started](admin/getting_started.md) * [Setting up an analysis](admin/setting_up_an_analysis_guidelines.md) * [Reproducible Research](admin/data_management.md) * [Methods for manuscripts](admin/method_snippets.md) * [Downloading data](admin/download_data.md) * [Acknowledging funding](admin/acknowledging_funding) *  [Initial consults](admin/initial_consults.md) * [Chargeback models](admin/chargeback_models.md)

**bcbio**: [docs](https://bcbio-nextgen.readthedocs.io/en/latest/) * [issues](https://github.com/bcbio/bcbio-nextgen/issues) * [running in parallel](https://bcbio-nextgen.readthedocs.io/en/latest/contents/parallel.html) * [genomes](bcbio/bcbio_genomes.md) * [tips](bcbio/bcbio_tips.md) * [RNA-seq workflow by Mary](bcbio/bcbio_workflow_mary.md) * [installation by Michael](https://steinbaugh.com/posts/install-bcbio.html) * [installation by Sergey](https://github.com/naumenko-sa/bioscripts/blob/master/bcbio/bcbio.upgrade.sh) * [RNA-seq workflow by Systems Pharmacology lab](https://labsyspharm.github.io/rnaseq/)

**Bulk RNA-seq**: [Count normalization methods](https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html) * [Dexseq](rnaseq/dexseq.Rmd) * [Failure types](rnaseq/failure_types) * [IRFinder_report](rnaseq/IRFinder_report.md) * [README](rnaseq/README.md) * [RepEnrich2 guide](rnaseq/RepEnrich2_guide.md) * [IRFinder](rnaseq/running_IRFinder.md) * [Strandedness](rnaseq/strandedness.md) * [Tools](rnaseq/tools.md) * [Bibliography](rnaseq/bibliography.md) * [Bibliography.bib](rnaseq/bcbio_rnaseq.bib) * [ASE](rnaseq/ase.md)

**CHIP-seq**: [ENCODE Guidelines](http://genome.cshlp.org/cgi/pmidlookup?view=long&pmid=22955991) * [CHIp-seq intro](https://hbctraining.github.io/Intro-to-ChIPseq/schedule/2-day.html) * [Tools](chipseq/tools.md)

**Misc**: [Git tips](misc/git.md) * [AWS](misc/aws.md) * [Core resources](misc/core_resources.md) * [General NGS](misc/general_ngs.md) * [Orphan improvements](misc/orphan_improvements.md) * [OSX](misc/OSX.md) * [Snakemake](misc/snakemake-example-pipeline) * [miRNA](misc/miRNA.md)

**R**: [.Rprofile](r/.Rprofile) * [RShiny](r/rshiny_server.md) * [htmlwidgets](/r/htmlwidgets) * [Tip and tricks](r/R-tips-and-tricks.md)

**RC**: [Connection to HPC](rc/connection-to-hpc.md) * [Keep command alive](rc/keepalive.md) * [IPython notebook](rc/ipython-notebook-on-O2.md) * [Manage files](rc/manage-files.md) * [O2 tips](rc/O2-tips.md) * [Scheduler](rc/scheduler.md) * [Jupyter notebooks](rc/jupyter_notebooks.md)

**[Single Cell RNA-seq](scrnaseq/README.md):** [Bcbio indrops3](scrnaseq/bcbio_indrops3.md) * [Rstudio docker](scrnaseq/rstudio_sc_docker.md) * [Saturation](scrnaseq/saturation_qc.md) * [Clustering analysis in Seurat](scrnaseq/seurat_clustering_analysis.md) * [Seurat markers](scrnaseq/seurat_markers.md) * [Single Cell conda](scrnaseq/Single-Cell-conda.md) * [Tinyatlas](scrnaseq/tinyatlas.md) * [Tools](scrnaseq/tools.md) * [Tutorials](scrnaseq/tutorials.md) * [Bibliography](scrnaseq/bibliography.md) * [Velocity](scrnaseq/velocity.md)

**Training team**: [github](https://github.com/hbctraining) * [Upcoming workshops](http://bioinformatics.sph.harvard.edu/training) * [Past workshops](http://bioinformatics.sph.harvard.edu/training#past-workshops) * [Workshop materials](https://hbctraining.github.io/main)
- Team members: Radhika Khetani, Mary Piper, Meeta Mistry, Jihe Liu
- Email: [hbctraining (at) hsph.harvard.edu](mailto:hbctraining@hsph.harvard.edu)

**WGS**: [CRISPR offtarget](wgs/crispr-offtarget.md) * [Pacbio genome assembly](wgs/pacbio_genome_assembly.md)

***

# How it is organized
- minimal overhead: feel free to add any piece of knowledge in any way
- easy to navigate in one click
- github is tracking changes for you
- it is ok to have multiple similar recipes
- growing naturally
- no placeholders
- please name all pages in small_letters.md
- order matters: topics come alpabetically, but within topic and within page please put important things first (important articles/tools on top, important links on the left).

# How to contribute
- nagivate to where you want to share your knowledge or fix something
- don't create a new folder unless you have >=2 new pages that don't fit into existing structure
- add actual value (know how), not structure
- hit 'edit file' or 'new file' (add **.md** extension)
- translate knowledge from your mind into markdown
- hit 'preview changes' to make sure your formatting is fine
- hit 'commit changes'
- thank you, you just have saved somebody's day!

# Knowledge we need to extract (aka TODO):
- ATAC-seq
- miRNA
- integrate (link?) hbctraining, tutorials, bcbio docs
- https://github.com/hbc/hbcABC
- https://github.com/vbarrera/bcbioLite
- Lorena's wisdom:
  - information that apply to all templates: http://bioinformatics.sph.harvard.edu/hbcABC/articles/general_start.html
  - List of templates right now: http://bioinformatics.sph.harvard.edu/hbcABC/articles/list_of_templates.html
  - single cell custom templates: http://bioinformatics.sph.harvard.edu/hbcABC/articles/singlecell-tutorial.html

