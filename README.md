MBASIC
======

OVERVIEW
--------

MBASIC ( *Matrix Based Analysis for State-space Inference and Clustering* ) package is a statistical framework for integrative analysis of many experiments where observations are collected over a set of units. The MBASIC framework allows for simultaneous projection of the observations onto a discrete state space and clustering the units based on their state-space profiles. 

MBASIC package implements specific functions for the analysis of ChIP-seq experiments. It enables the integration of multiple ChIP-seq data across different experimental conditions (i.e. celltypes and transcription factors). Given a set of loci, MBASIC is able to identify the enrichment signal for each locus under each experimental condition, as well as extract clusters of loci that exhibit highly similar enrichment patterns across all experimental conditions.

MBASIC is being currently used for multiple on-going research projects. It has exhibited the following advantages:

- Adaptable to the varying number of replicates for each experimental condition;

- Flexible in distributional assumptions for the actural data: we have implemented log-normal and negative binomial distributions, but we will extend to more distributions;

- High computational competency: one data set we currently analyze ranges over 10K loci and 200 conditions.


INSTALLATION
------------


MBASIC will be available at Bioconductor. Currently you can download the development version here and install in R by:

    library(devtools)
    install_github("chandlerzuo/mbasic")


REFERENCES
----------

Chandler Zuo and Sunduz Keles (2015). "A hierarchical framework for state-space matrix inference and clustering". *To appear*.
