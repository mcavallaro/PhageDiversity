This repository contains software used in reference [1]. Please cite [1] if you find this repository useful.

<p align="center">
<img src="https://github.com/mcavallaro/PhageDiversity/blob/clean/cover.png?raw=true" alt="" style="width:70%; height:auto;">
</p>

The software is organised as follows.
  
  * General utility functions for working with species distributions.
    - `utils.R`

  * Scripts with estimators.
    - `FPG_estimator.R`
    - `PYP_estimator.R`
    - `nonparam_estimators.R` 

  * Scripts for validation.
    - `internal_valid_fixedm.R`      
    - `internal_valid_from25set.R`   
      * visualisation and summaries
        - `plots_intval.R`
        - `plots_intval_fixedm.R`
        - `get_table_intval.R`
 
  * Scripts to illustrate diversity and predictions.
    - `diversity.R`
    - `inext_asymp_plot.R`
    - `histograms.R`
    - `predict_ET_monotone.R`
    - `predict_DB24.R`
    - `predict_DB25.R`
    - `add1Ksampling.R`
    - `dump24vs25.R`
    - `predict_cover.R`
    
On top, we provide the ouput data of these scripts (`.RData` objects). The folder `data/` contains
the two database snapshots DB24 and DB25 (can be read with ` import_data.R`) analysed in reference [1].

[1] "How many phage species remain undiscovered? Species sampling approaches to inform phage discovery"
Massimo Cavallaro, Andrew Kinsella, Spyridon Megremis, Andrew Morozov, Andrew D. Millard, Fabian Freund
*bioRxiv* 2026.02.15.704868; doi: \[https://doi.org/10.64898/2026.02.15.704868][https://doi.org/10.64898/2026.02.15.704868]