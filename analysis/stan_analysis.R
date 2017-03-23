# load packages we'll use ----
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)

library(devtools)
install_github('mike-lawrence/ezStan')
library(ezStan)

# Read in the data and look at it ----
dat = readr::read_csv('dat.csv')
View(dat)
