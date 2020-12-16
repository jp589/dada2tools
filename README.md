# dada2tools
DADA2 Tools

This package was written for researchers using DADA2 for their 16s rRNA sequencing data analysis. It allows easy construction of PCoA plots, heatmaps, barplots, and more. The vignette walks the reader through each of the functions in the order that they should be run with detailed examples and helpful notes.

Installation instructions:

1. Install devtools and load it in R console with:
  `BiocManager::install("devtools"); library(devtools)`
2. Install the `dada2tools` package hosted on github:
  `install_github("jp589/dada2tools", build_vignettes = TRUE)`
3. Load the package with:
  `library(dada2tools)`
4. Start the vignette to discover functions in the package:
  `vignette("dada2tools")`
5. For function descriptions:
  `help(package = 'dada2tools', help_type = 'html')`
