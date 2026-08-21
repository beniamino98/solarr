





















library(solarr)


spec <- solarModel_spec$new()
spec$set_mean.model(arOrder = 1, maOrder = 1)
spec$specification("Bologna")
spec



spec$set_seasonal.mean(include.trend = TRUE, monthly.mean = TRUE)
spec

spec$set_mean.model(2, 1)
spec

spec$set_mixture.model(match.expectation = TRUE, match.variance = TRUE, match.empiric = TRUE)


###
###

spec_train <- Bologna$spec

model <- solarModel$new(spec_train)
model$filter()
model$update_moments()
model$moments$conditional[454,]
Bologna$moments$conditional[454,]

model$update_logLik()
model$loglik
Bologna$loglik

model$spec$mean.model
Bologna$spec$mean.model

model$spec$mixture.model$model[5]
Bologna$spec$mixture.model$model[5]

spec <- solarModel_spec$new()
spec$specification("Bologna", to = "2022-01-01")
spec
spec$variance.model
spec$seasonal.mean

model <- solarModel$new(spec)
model$fit()
model$spec$hessian

model <- solarModel_QMLE(model)
model$.__enclos_env__$private$..hessian
model$.__enclos_env__$private$..jacobian


model$.__enclos_env__$private$..data
spec_train <- model$spec$clone(TRUE)

model2 <- solarModel$new(spec_train)
model2$filter()
model2$update_classification()
model2$update_moments()
model2$update_logLik()

model2

solarOption_model(model, model$moments$conditional[4:565,])
solarOption_model(model2, model2$moments$conditional[4:565,])

