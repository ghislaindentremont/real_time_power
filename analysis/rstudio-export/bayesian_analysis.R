# ---
# title: "Bayesian Analysis"
# output: html_notebook
# author: "Ghislain d'Entremont"
# ---

# Load and Tidy Data

##### Load

# load packages we'll use ----
library(tidyverse)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)

library(devtools)
install_github('mike-lawrence/ezStan')
library(ezStan)


##### Read 

# Read in the data and look at it ----
dat = readr::read_csv('dat.csv')
View(dat)

##### Sort
dat %>%
	dplyr::arrange(
		subj
	) -> dat


##### Within Contrasts

W = get_contrast_matrix(
	data = dat
	, formula = ~  task * trial * laterality 
)
View(W)


##### Between Contrasts

# I don't actually have any between-subject variables for the time being. 

dat %>%
	dplyr::group_by(
		subj
	) %>%
	dplyr::summarize() %>% #View()
	get_contrast_matrix(
		formula = ~ 1
	) -> B
View(B)



## Run COMPLEX model and save

#package in list for Stan
data_for_stan = list(
	#nTrials: num trials
	nTrials = nrow(dat)
	#power: power outcomes
	, power = scale(dat$log_relative_power)[,1] #scaling for easy priors
	#nS: num subjects
	, nS = length(unique(dat$subj))
	#S: trial-by-trial subject labels
	, S = as.numeric(factor(dat$subj))
	#nWmpower: num within predictors on mean power
	, nWmpower = ncol(W)
	#nBmpower: num group predictors on mean power
	, nBmpower = ncol(B)
	#nWspower: num within predictors on sd power
	, nWspower = ncol(W)
	#nBspower: num group predictors on sd power
	, nBspower = ncol(B)
	#Wmpower: within predictors for mean power
	, Wmpower = W
	#Bmpower: between predictors for mean power
	, Bmpower = B
	#Wspower: within predictors for sd power
	, Wspower = W #for just intercepts do: matrix(W[,1],ncol=1)
	#Bspower: between predictors for sd power
	, Bspower = B #for just intercepts do: matrix(B[,1],ncol=1)
)


# Compile & sample the model ----
post_c = rstan::stan(
	file = 'bayesian_analysis.stan'
	, data = data_for_stan
	, seed = 1
	, chains = 16
	, cores = 16 #set this to # of physical cores on your system
	, iter = 2e3
	, init = 0
	, pars = c('normal01','cors_helper') #don't bother saving these variables
	, include = FALSE
)
save(post_c,file='post_c.rdata')

##### complex output

load('post_c.rdata')

stan_summary(
	from_stan = post_c
	, par = 'coefMpower'
	, W = W
	, B = B
)

stan_summary(
	from_stan = post_c
	, par = 'coefSpower'
	, W = W
	, B = B
)

# this is still scaled
stan_summary(
	from_stan = post_c
	, par = 'sdsW'
)

stan_summary(
	from_stan = post_c
	, par = 'corsW'
)


