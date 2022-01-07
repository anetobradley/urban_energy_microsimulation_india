# Population Synthesis <!-- omit in toc -->

This folder contains an R Scripts for using Iterative Proportional Fitting (IPF) to generate a synthetic population for a given city in India.

- [1. Data Sources](#1-data-sources)
- [2. Population Synthesis Overview](#2-population-synthesis-overview)


## 1. Data Sources

> Datasets used are publicly available although you may have to register on the respective platform to download the data.

The contsraint tables are wrangled from the ward level asset ownership tables of the Indian Census which contain information of appliance ownership, access to banking and home construction.

The sample of representative individual level instances (households) can be sourced from the Indian Human Development Survey. Information about this survey as well as the data can be accessed online (https://ihds.umd.edu/)

This was originally carried out using 2011 Census and 2011-12 IHDS data although ta third round of the IHDS was carried out in 2021 as well as the Census and so from 2023 there should be more recent data avilable for generating these synthetic populations.

## 2. Population Synthesis Overview

The aim of the population synthesis is to generate a representative population of households for each ward in the city of interest with socio-economic attributes that may be used to estimate fuel consumption in the absnece of such fine scale data.

The procedure carried out by the R script in this folder performs the following steps:


