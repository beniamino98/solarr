# Document
#devtools::document()
# Build manual
#devtools::build_manual(path = "/Users/macbook/Documents/University/PhD/Projects/solar-project/R/solarr")

# Sample model
# Control list
control <- control_solarModel(outliers_quantile = 0, seasonal.mean = list(monthly.mean = TRUE),
                              mean.model = list(arOrder = 1, maOrder = 1), garch_variance = TRUE)
# Model specification
#spec <- solarModel_spec("Bologna", from="2005-01-01", control_model = control)
#Bologna <- solarModel$new(spec)
# Model fit
#Bologna$fit()
#save(Bologna, file = "data/Bologna.RData")
