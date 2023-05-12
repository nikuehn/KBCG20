# KBCG20

This is a repository that includes files/information for the BGA subduction ground-motion model of Kuehn, Bozorgnia, Campbell, and Gregor (2020) (KBCG20). The repository contains the following directories:

* COEEFS: contains files with coefficients
  * posterior_coefficients_KBCG20_TXX.XXX.csv: contains 800 sets of coefficients/parameters for one period.
  * coefficients_KBCG20.csv: contains smoothed mean coefficients/parameters for all periods.
  * Note: These coefficents are superseded by COEFFS_SEP21, which contain updated coefficients for Alaska based on updated Alaska data.
* DATA: contains files with coefficients
  * data_NGAsub_KBCG20.csv: contans record sequence numbers, earthquake ids, and station sequence numbers for the record used in the regression.
* RESIDUALS: contains plots of event terms against magnitude and depth-to-top-of-rupture; as well as within-event residuals against distance, Vs30, and basin depth.
* Stan: contains Stan files for regression (see http://mc-stan.org for info)
* UNCERTAINTY: Contains values of the standard deviations of epistemic uncerainty for different regions, magnitudes, distances, and tectonic environments, calculated from predictions based on the posterior samples. Can be used for a simplified application of epistemic uncertainty, such as a scaled backbone.
* pictures: contains pictures used in R-examples

The files KBCG20 contain an R-implementation of the functional form for KBCG20. THis function is only used to show he implementation, and compare predictions.

Update September 2021:
The drectory COEFFS_SEP21 contains new sets of cofficiets, based on updated Alaska data.
Only the regional coefficients for Alaska have changed.
These are the coefficients that should be used.
