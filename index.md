---
layout: default
title: Home
---

# SLAPenrich

## Overview

SLAPenrich (**S**ample-population **L**evel **A**nalysis of **P**athway **enrich**ments) is a pathway-level analysis tool set for the characterization of genomic data from heterogeneous populations under positive selection, such as cancer patients.

SLAPenrich implements a statistical framework to identify pathways that tend to be recurrently genomically altered across the samples of a genomic dataset. Differently from traditional over-recurrence analyses, SLAPenrich does not require the genes belonging to a given pathway to be statistically enriched among those altered in the individual samples. Consistently with the mutual exclusivity principle, and differently from other proposed computational tools, our approach assumes that the functionality of a given pathway might be altered in an individual sample if at least one of its genes is genomically altered. The method accounts for the differences in the **mutation rates** between samples and the **exonic lengths** of the genes in the pathways. It statistically tests against the null hypothesis that no associations between a pathway and the disease population under study does exist, assessing analytically the divergence of the total number of samples with alterations in a given pathway from its expectation. Moreover, the used formalism allows SLAPenrich to perform **differential enrichment analysis** of pathway alterations across different clinically relevant sub-populations of samples. SLAPenrich also includes function to **visualise** the identified enriched pathway implementing a heuristic sorting to highlight mutual exclusivity trends among the pattern of alterations of the composing genes.

This package contains the R package and accompanying scripts that implement the method as well as several examples of how to run a SLAPenrich analysis.

## References
Please use this reference to cite SLAPenrich:
>  Iorio F., Garcia-Alonso L., Brammeld J., Martincorena I., Wille D. R.,  McDermott U.,  Saez-Rodriguez J. (2017). [Dissecting the genomic heterogeneity of cancer hallmarks' acquisition with SLAPenrich](http://biorxiv.org/content/early/2017/03/27/077701). _Under Revision_
   

### Applications of SLAPenrich:

 > Brammeld J, Petljak M, Martincorena I, Williams SP, Garcia-Alonso G, Dalmases A, Bellosillo B, Robles-Espinoza CD, Price S, Barthorpe S, Tarpey P, Alifrangis C, Bignell G, Vidal J, Young J, Stebbings L, Beal K, Stratton MR, Saez-Rodriguez J, Garnett M, Montagut C, Iorio F, McDermott U. (2017) [Genome-wide chemical mutagenesis screens allow unbiased saturation of the cancer genome and identification of drug resistance mutations](http://genome.cshlp.org/content/early/2017/03/15/gr.213546.116.abstract?cited-by=yes&legid=genome;gr.213546.116v2) _Genome Res._


