# Principal Component Analysis (PCA) of Star-Formation Histories (2D-data) of LEGA-C galaxies derived from two Bayesian Tools 

_Keywords_ - PCA, SFH, Galaxies, photometric data, spectrophotometric data, <br>

_References_ - This analysis makes extensive use of [Sparre et al 2015](https://academic.oup.com/mnras/article/447/4/3548/1751599?login=false) method of performing PCA on 2D SFH data. This analysis is also built on my recent publication to understand the similarity in SFHs from two Bayesian SED fitting tools each optimized for LEGA-C dataset. [Kaushal et al 2023](https://doi.org/10.3847/1538-4357/ad0c4e). <br>

## Context 

We can obtain two types of observational data of galaxies - <br>

1. Photometric data (small) 
2. Spectroscopic data (large)
   
Bayesian modeling of these data together will give us information about their Star-formation Histories (SFHs) as well as other internal properties like dust content, metal content, ages of stars, star-formation rates etc. <br>

My recent publication focused on comparing the SFHs of these galaxies at 3 physical timestamps, namely - t10, t50 and t90 - corresponding to the times when a galaxy forms 10%, 50% and 90% of its stellar mass [Kaushal et al 2023](https://doi.org/10.3847/1538-4357/ad0c4e). <br>

In this notebook, we explore a step further and systematically study major star formation modes of galaxies of the LEGA-C Survey at - **not just at few time-stamps, but across their full lifetime/history.** <br>
For that, we adopt PCA technique described in detail in [Sparre et al 2015](https://academic.oup.com/mnras/article/447/4/3548/1751599?login=false) (inspired by Cohn & Voort 2015). <br>
SFH of a galaxy is seen as a vector in an N-dimensional space. <br>
As >90% of the information could be represented in first 3 components, we will study PC0, PC1 and PC2 only. Note that they describe scatter around the mean SFH and diagonalize the scatter matrix. <br>

They characterize the **most important modes of star-formation hisotry in the population** <br>
The PC0-mode accounts for galaxies forming stars early or late, depending on whether the coefficient q0 is positive (earlier than sample mean) or negative (later than sample mean) <br>
PC1 and PC2 determine the more detailed evolution of the SFH and cross the zero-point two and three times, respectively. 

## Questions to answer 

We would like to statistically understand the Prinacipal Components of SFHs of population of galaxies in the LEGA-C Survey. <br>
We will compare Principal Components of QUIESCENT and STAR-FORMING galaxies _separately_ in bins of stellar mass (low to high mass) derived from two Bayesian SED Modeling tools - BAGPIPES and PROSPECTOR <br>
We would like to answer - <br>
1. How similar / different Principal Components (PC0, PC1 and PC2) of recovered SFHs are from (1) Photometric data (2) Spectro-photometric data. 
2. How much more information does spectroscopic data add to SFHs? 
3. Are the inferences similar/different from the two Bayesian modelings? 
4. Is there a correlation between the recovered properties of galaxies (age, metallicity, dust, ssfr) with the SFH? 

We chose **BvrizYJ bands** for photometric modeling and **LEGA-C spectra + BvrizYJ bands** for spectro-photometric modeling. <br>

## Data 

1. LEGA-C Spectroscopic data (https://iopscience.iop.org/article/10.3847/1538-4365/ac1356/meta) 
2. UltraVista Photometric data (https://iopscience.iop.org/article/10.1088/0067-0049/206/1/8)

## Results 

1. For Bayesian code 1 (BAGPIPES) -
   a. SFHs from photometry only data are not able to recover higher moments of SFH as spectro-photometric data does.
   b. Adding spectroscopic information is leading to expected correlation between ages (both mass and light weighted ages) and SFH PC0 coefficient q0, more so in the star-forming galaxies than in the quiescent galaxies.
   c. Contours of photometry only data are primary concentrated around q0 = 0 indicating weak dependance on SFH, except the most massive bin (11.5 < M < 12) where quiescent galaxies show opposite dependence on the SFH.
   d. Ages and Metallicites of most massive quiescent systems from photometric only data are not well recovered and require spectroscopic information. 

2. For Bayesian code 2 (PROSPECTOR) -
   a. We know that PC0 contains the most information of the SFH variance and PC1/PC2 represent information in the higher modes.
   b. Y-axis represents the direction of maximum information variance with respect to the mean SFH of the population.
   c. Adding spectra definitely leads to better recovery of finer variations / higher moments of SFHs in the population, which are missed by photometric data.
   d. Most massive galaxies (11.5 < M < 12) show the largest difference in Principal Components from photometric data and spectrophotometric data.
   e. Early variablility of SFH is captured similarly in star-forming, partly owing to their young stellar populations. 


