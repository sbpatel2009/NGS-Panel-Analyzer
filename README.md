# NGS-Panel-Analyzer

## Installation

```
source('https://bioconductor.org/biocLite.R')
biocLite('sbpatel2009/NGS-Panel-Analyzer')
```

## Usage


1. Manually download the COSMICMUTANTEXPORT.TSV, VERSION 83, GRCh37 file. Call it `cosmic.tsv` and
place it in the R working directory.
2. Load the package: `library(NGSPanelAnalyzer)`
3. Run the app: `runNGSPanelApp()`
