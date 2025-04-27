library(tidyverse)
library(solarr)
# Test the residuals under N(0,1) or GM distribution
model <- Bologna

solarModel_test_residuals(model, ci = 0.05, nrep = 100, seed = 1)
