# Bayesian statistics 
#### Olivier Gimenez, July 2020

## Slides codes and data

* All material prepared in `R`.
* `R Markdown` used to write reproducible material (`R` code also available).
* Slides available on FigShare [here](https://doi.org/10.6084/m9.figshare.12656894.v2).
* Material available via Github [there](https://github.com/oliviergimenez/Bayesian_Workshop).  

## Objectives

* Try and demystify Bayesian statistics, and what we call MCMC.
* Make the difference between Bayesian and Frequentist analyses.
* Understand the Methods section of ecological papers doing Bayesian stuff.
* Run Bayesian analyses (in Jags), safely hopefully.

## Schedule

1. Bayesian inference: Motivation and simple example.
2. The likelihood.
3. A detour to explore priors.
4. Markov chains Monte Carlo methods (MCMC).
5. Bayesian analyses in R with the Jags software.
6. Contrast ecological hypotheses with model selection.
7. Heterogeneity and multilevel models (aka mixed models).

## To do list

   * Add a plot with several lines from posterior distribution of regression parameters to a plot of mean response function of a covariate; then get the credible interval on the prediction. 
   * Besides `Jags` and `Nimble`, there are other software options to fit models in the Bayesian framework that do not need coding. Check out the [CRAN Task View: Bayesian Inference](https://cran.r-project.org/web/views/Bayesian.html).
   * Mention free Bayesian books: [here](https://www.bookdown.org/home/tags/bayesian/) and Gelman BDA [there](http://www.stat.columbia.edu/~gelman/book/).  





* Add a section on posterior predictive checks, to comply with the 3 steps of a Bayesian analysis as defined by Gelman (set up a probabilistic model, inference and model checking; iterate to improve model).
* Say more on prior predictive checks. 
* Say something about confidence, credible and HPD intervals.
* Add another Metropolis example, with adaptation, with the beta-binomial example, and discuss several levels of acceptance. 
* Add a section on LOO, and discuss complementarity with WAIC.
* Add a section on models with varying slopes. Can we use the [LKJ prior](https://www.sciencedirect.com/science/article/pii/S0047259X09000876) in `Jags` and `Nimble`?




* Write a short introduction to `Nimble` and provide both the `Jags` and `Nimble` codes. Translating `Jags` code in `Nimble` is [easy](https://r-nimble.org/quick-guide-for-converting-from-jags-or-bugs-to-nimble). For now, check out [training materials](https://github.com/nimble-training) and [examples](https://r-nimble.org/examples).
* Do all plots with `ggplot2`; add [short introduction to the `Tidyverse`](https://github.com/oliviergimenez/intro_tidyverse).
* Add a section on population ecology (occupancy models, capture-recapture models). And/or something on hierarchical models, models with hidden variables. Make use of [nimbleEcology](https://cran.r-project.org/web/packages/nimbleEcology/vignettes/Introduction_to_nimbleEcology.html).
* Add a section on [penalized splines](https://www.cambridge.org/core/books/semiparametric-regression/02FC9A9435232CA67532B4D31874412C) (possibly using package `jagam`) and [spatial analyses](https://r-nimble.org/html_manual/cha-spatial.html).



* ~~Add a short section on sequential analysis (today prior is yesterday posterior).~~
* ~~Add an example with Poisson GLM(M) example.~~
