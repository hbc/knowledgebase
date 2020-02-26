# Doublet identification

- gene number filter is not effective in identifying doublets (Scrublet2019 article).
- DoubletIdentification works for a group of cells we suspect they might be doublets (a cluster or a group of clusters)  - if we see expression of two cell types, 'mixed' marker signature. 
- dump counts from suspected counts from Seurat
- identify doublets with Scrublet
- get back to Seurat

R based DoubletFinder and DoubletDecon have issues
- https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/64
- https://github.com/EDePasquale/DoubletDecon/issues/21

