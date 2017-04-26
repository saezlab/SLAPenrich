---
layout: page
title: Usage
---



To start using SLAPenrich, first load the sparse **binary matrix data (EM)**. See the example `LUAD_CaseStudy` included in the package:

```r
data(LUAD_CaseStudy)
class(LUAD_CaseStudy)
LUAD_CaseStudy[1:5, 1:5]
```

Then load the **list of pathway gene sets (PATH_COLLECTION)**. SLAPenrich provides two collections of pathways from two databases: KEGG and PathwayCommons. For example, to load the pathways from KEGG type:

```r
data(SLAPE.MSigDB_KEGG_hugoUpdated)
class(KEGG_PATH)
names(KEGG_PATH)
```

Finally, load information on the **exonic lentgh of the genes (GeneLenghts)**:

```r
data(SLAPE.all_genes_exonic_content_block_lengths_ensemble)
class(GECOBLenghts)
head(GECOBLenghts)
```


To run the analysis, type the `SLAPE.analyse` function:

```r
mySLAPE_analysis = SLAPE.analyse(EM=LUAD_CaseStudy, PATH_COLLECTION=KEGG_PATH, GeneLenghts = GECOBLenghts)
class(mySLAPE_analysis)
names(mySLAPE_analysis)
```


Plot your core components summary with `SLAPE.core_components`. This will generate a series of pdf files in the PATH folder. Each pdf file will represent a core component (i.e. a group of related patways):

```r
SLAPE.core_components(PATH = "./", PFP=mySLAPE_analysis, EM=LUAD_CaseStudy, PATH_COLLECTION=KEGG_PATH)
```



To save your results as a table, run `SLAPE.write.table`:

```r
SLAPE.write.table(filename = "mySLAPE_analysis.csv", PFP=mySLAPE_analysis, EM=LUAD_CaseStudy, PATH_COLLECTION=KEGG_PATH, GeneLenghts = GECOBLenghts)
```

Check here for a full [documentation](https://github.com/saezlab/SLAPenrich/blob/master/SLAPenrich.pdf) of SLAPenrich.




