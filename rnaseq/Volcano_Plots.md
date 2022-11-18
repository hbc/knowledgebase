## [Enhanced Volcano](https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html) is a great and flexible way to creat volcano plots.

Input is a dataframe of test statistics. It works well with the output of `lfcShrink()`
Below is an example call and output:

```
library(EnhancedVolcano)

EnhancedVolcano(shrunken_res_treatment,
                lab= NA,
    x = 'log2FoldChange',
    y = 'pvalue', title="Volcano Plot for Treatment", subtitle = "")
```

<p align="center">
<img src="/rnaseq/img/volcano.png" width="800">
</p>


Almost every aspect is flexible and changable. 
