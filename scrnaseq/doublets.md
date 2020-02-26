# Doublet identification

- gene number filter is not effective in identifying doublets (Scrublet2019 article).
- there is no good unsupervised doublets detection method for now
- DoubletIdentification works for a group of cells we suspect they might be doublets (a cluster or a group of clusters)  - if we see mixed marker signature and we know that these cells are not in transitional state, i.e. expert review of clusters is needed before doublet deconvolution
- dump counts from suspected counts from Seurat
- identify doublets with Scrublet
- get back to Seurat

R based DoubletFinder and DoubletDecon have issues
- https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/64
- https://github.com/EDePasquale/DoubletDecon/issues/21

