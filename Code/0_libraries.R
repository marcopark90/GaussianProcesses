library(tidyverse)
library(Metrics)
library(magrittr)
library(rstan)
library(bayesplot)
library(scales)
library(gridExtra)
library(ChainLadder)
library(fields)

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)
