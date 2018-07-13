# Import/Export of files
Stop using write.csv, write.table and use the [rio](https://cran.r-project.org/web/packages/rio/index.html) library instead. All rio needs is the file extension to figure out what file type you're dealing with. Easy import and export to Excel files for clients.

# Parsing in R using Tidyverse
This is a link to a nice tutorial from Ista Zahn from IQSS using stringr and tidyverse for parsing files in R. It is from the Computefest 2017 workshop:
http://tutorials-live.iq.harvard.edu:8000/user/zwD2ioESyGbS/notebooks/workshops/R/RProgramming/Rprogramming.ipynb

# Better clean default ggplot
install cowplot (https://cran.r-project.org/web/packages/cowplot/index.html)
```r
library(cowplot)
```

# Nice looking log scales
Example for x-axis
```r
library(scales)
   p + scale_x_log10(
         breaks = scales::trans_breaks("log10", function(x) 10^x),
         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
       annotation_logticks(sides='b')
```

# Read a bunch of files into one dataframe
```r
library(purrr)
library(readr)
library(dplyr)
library(tidyr)
read_files = function(files) {
  data_frame(filename = files) %>%
    mutate(contents = map(filename, ~ read_tsv(.))) %>%
    unnest()
}
```

# remove a layer from a ggplot2 object with ggedit
```
plotGeneSaturation(bcb, interestingGroups=NULL) +
  ggrepel::geom_text_repel(aes(label=description, color=NULL)) 
p %>%
  ggedit::remove_geom('point', 1) +
  geom_point(aes(color=NULL))
```

# [Link to information about count normalization methods](https://github.com/hbc/knowledgebase/wiki/Count-normalization-methods) 
The images currently break, but I will update when the course materials are in a more permanent state.