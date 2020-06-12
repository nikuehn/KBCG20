# KBCG20

This is a repository that includes files/information for the BGA subduction ground-motion model of Kuehn, Bozorgnia, Campbell, and Gregor (2020) (KBCG20). The repository contains the following directories:

* DATA: contains files with coefficients
  * posterior_coefficients_KBCG_TXX.XX.csv: contains 800 sets of coefficients/parameters for one period.
  * coefficients_KBCG19.csv: contains smoothed mean coefficients/parameters for all periods.
  * data_NGAsub_KBCG19.csv: contans record sequence numbers, earthquake ids, and station sequence numbers for the record used in the regression.
* RESIDUALS: contains plots of event terms against magnitude and depth-to-top-of-rupture; as well as within-event residuals against distance, Vs30, and basin depth.
* Stan: contains Stan files for regression (see http://mc-stan.org for info)
* pictures: contains pictures used in R-examples

The files KBCG20 contain an R-implementation of the functional form for KBCG20.
