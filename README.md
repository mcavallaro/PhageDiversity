This repository contains software used in reference [1]. Please cite [1] if you find this repository useful.

[1] How many phage species remain undiscovered? Species sampling approaches to inform phage discovery.
Massimo Cavallaro, Andrew Kinsella, Spyridon Megremis, Andrew Morozov, Andrew D. Millard, Fabian Freund
bioRxiv 2026.02.15.704868; doi: https://doi.org/10.64898/2026.02.15.704868 


<p align="center">
<img src="https://github.com/mcavallaro/PhageDiversity/blob/clean/cover.png?raw=true" alt="" style="width:70%; height:auto;">
</p>

The software is organised as follows (all methods are referenced in [1] above).
  
  * General utility functions for working with species distributions.
    - `utils.R`

  * Scripts to source estimators for the number of unseen species in further samples.
    - `FPG_estimator.R`: parametric Fisher-Poisson-Gamma estimator (FPG)
    - `PYP_estimator.R`: parametric estimator based on the Pitman-Yor distribution family (PYP)
    - `nonparam_estimators.R`: non-parametric estimators Efron-Thisted (ET), Chao-Jost (CJ), Orlitsky-Suresh-Wu (OSW),
    Good-Toulmin (GT)
    - `ugland_logmodel.R`: estimator via species accumulation curve extrapolation following an approach of Ugland et al.
  
  * Scripts for assessing goodness-of-fit for PYP and FPG estimates 
    - fit_FPG.R  
    - fit_PYP.R
  
  * Script for data procurement from subfolder data/ (sourced from other scripts)
    - import_data.R

    
  * Scripts for validation via internal random subsetting.
    - `internal_valid_from25set.R`: for estimators FPG, PYP, GT, ET
        - internal_valid_from25set_chaojost.R: adds CJ
        - internal_valid_from25set_orlitzkyestim.R: adds OSW
        - internal_valid_from25set_ugland.R: adds limited analysis for semi-log species accumulation curve extrapolation (see Supplementary text 2)
    - `internal_valid_fixedm.R`: prediction using FPG and OSW for fixed-size sub-samples for training and prediction sets        
      * visualisation and summaries of internal validation results
        - `plots_intval.R`
        - `plots_intval_fixedm.R`
        - `get_table_intval.R`
 
  * Scripts to illustrate diversity and predictions.
    - `diversity.R`: Analysing species diversity using Hill numbers (Supplementary text 4)
    - `inext_asymp_plot.R`: Analysing species diversity using Hill numbers (Supplementary text 4)
    - `histograms.R`: Descriptive analysis of species abundances in the 8 host genera assessed in [1]
    - `predict_OSWFPGCJ_bt.R`: Predicting unseen species in future samples from DB25 using different estimators (as in [1]), with and without bootstrapping
        - `predict_OSWFPGCJ_loglin_monotone.R`: *NOT IN [1]* Adds further, naive semi-log and log-log extrapolations of the species accumulation curve
    - `add1Ksampling.R`: *NOT IN [1]* Different assessments of predicting unseen species i) a different/random genus w. 1,000 phages already sampled 2) a new, unsampled host
    - `dump24vs25.R`: Predict the additional (unseen) species in DB25 from DB24 using different estimators 
	- `predict_cover.R`: Prediction of the number of new species in future samples using the OSW estimator (produces the cover image of this repo)
    
On top, we provide the ouput data of these scripts (`.RData` objects). The folder `data/` contains
the two database snapshots DB24 and DB25 (can be read with ` import_data.R`) analysed in reference [1].
