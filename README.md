# Urban residential energy microsimulation toolkit for India
Tools for generating urban synthetic populations for India, and running energy transition microsimulations.

REPO IS STILL BEING ASSEMBLED - DOCUMENTATION AND SOURCE CODE BEING UPDATED.

- [1. Data Sources](#1-data-sources)
- [2. Objectives & Functionality](#2-objectives-and-functionality)
- [3. Workflow Overview](#3-workflow-overview)
  - [3.1 Population Synthesis](#31-population-synthesis)
  - [3.2 Bayesian Multi-level Model](#32-bayesian-multi-level-model)
  - [3.3 Spatially Visualise Results](#33-spatially-visualise-results)

## 1. Data Sources
This toolkit is designed to work with the following public data sources which can be found online.

- Census Ward Level Tables
- Indian Human Development Survey (IHDS) 
- NSS Consumer Expenditure Survey (NSS)
- Shapefile of Census City Ward Map

Sample results and outputs have been generated using data fora given city from these sources.

## 2. Objectives and Functionality

This microsimulation aims to estimate likely cooking fuel use and fuel stacking likelihood for city subdivisions (wards), taking into account households practices, local socio-cultural context, and city economic context. To do this a representative synthetic population is generated for each ward in the city and then a bayesian multilevel model is use to estimate fuel consumption for each based on the primary cooking fuel group of the household, and the ward it is in.

## 3. Workflow Overview

![20210414_Microsim_Overview](https://user-images.githubusercontent.com/66263560/115389727-5ae51280-a1d5-11eb-89b8-db5de217be53.png)

### 3.1 Population Synthesis
Population synthesis is carried out using Iterative Proportional Fitting and implemented in R using the `ipfp` package (https://spatial-microsim-book.robinlovelace.net/). Ward-level Census Tables are used as contraint tables and a sample of household level responses from the IHDS are used as microdata to populate the synthetic population from.

### 3.2 Bayesian Multi-level Model
A multilevel model is used to estimate cooking fuel use to allow model coefficients to vary by primary cooking fuel group and city ward, thus capturing effects of household practices, and local socio-cultural and spatial context.

#### Primary Fuel Choice Estimation
Primary fuel choice estimation is carried out using a categorical logit using income, majority religion membership, ration card status, home ownership status, and employment type as predictors.

#### Hierarchical Model
A Bayesian Hierarchical Model is used to estimate cooking fuel use with coefficients varying by cooking fuel group `i` and cityward `j`.

#### Spatial Effects
Two different types of coefficient are used in the model to capture local spatial effects. 

- Random Effects Coefficient: This coefficient captures how the given city varies from the state-wide urban average. This can be due to economic, demographic, social or other factors.
- Intrinsic Conditional Autoregressive (ICAR) Coefficient: This coefficient captures how each individual ward compares to its neighbouring ones, and reflects socio-economic and cultural patterns across city wards.

### 3.3 Spatially Visualise Results
A key benefit of this city scale microsimulation is that results can be visualised at a city ward level. This not only makes results easy to communicate but also can effectively show patterns of inequality.
