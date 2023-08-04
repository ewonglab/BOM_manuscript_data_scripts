## [Tutorial](https://ewonglab.github.io/BOM/)

For a comprehensive guide on fully utilizing BOM, please refer to the tutorial ('binary_tutorial.ipynb') provided here.

In this [tutorial](https://ewonglab.github.io/BOM/), we will go through the whole BOM pipeline and use BOM to train a binary model, which is able to predict mouse cardiomyocyte-specific candidate regulatory elements (CRE) from background CRE of other 16 cell types. 

<!---
![BOM_pipeline](///BOM_pipeline.png)
-->

Before you run through this tutorial, make sure you execute the below commands first.

```
git clone https://github.com/ewonglab/BOM.git
cd BOM/Tutorial
unzip Pijuan_etal_table_S6.csv.zip
```

> The snATAC-seq data resources can be found at Tutorial/data/Pijuan_etal_table_S6.csv.zip, from Pijuan-Sala et al. 2020 [1]
> In the data file, the cell type annotation is indicated in the column "celltype_specificity". Some snATAC-seq peaks were annotated to multiple cell types. Please find more details in the tutorial.
