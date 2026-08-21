

# Test sGARCH reparametrization

# Test map (ARCH 1)
test_that('sGARCH_params_to_: check GARCH(1,0) ...', {
  # Original parameters
  params <- c(omega = 0.8, alpha1 = 0.2)
  # Reparametrized
  zeta <-  sGARCH_params_to_zeta(params, 1, 0, mode = "unitOmega")
  # Map back
  params0 <- sGARCH_params_to_phi(zeta, 1, 0, mode = "unitOmega")
  # Jacobian
  J0 <- sGARCH_params_to_zeta_jacobian(zeta, 1, 0, mode = "unitOmega")
  expect_true(dim(J0)[1] == 2 & dim(J0)[2] == 1)
  # Check equality
  expect_true(sum(abs(params - params0)) == 0)
})

# Test map (GARCH 1)
test_that('sGARCH_params_to_: check GARCH(0, 1) ...', {
  # Original parameters
  params <- c(omega = 0.8, beta1 = 0.2)
  # Reparametrized
  zeta <-  sGARCH_params_to_zeta(params, 0, 1, mode = "unitOmega")
  # Map back
  params0 <- sGARCH_params_to_phi(zeta, 0, 1, mode = "unitOmega")
  # Jacobian
  J0 <- sGARCH_params_to_zeta_jacobian(zeta, 0, 1, mode = "unitOmega")
  expect_true(dim(J0)[1] == 2 & dim(J0)[2] == 1)
  # Check equality
  expect_true(sum(abs(params - params0)) <= 10e-10)
})

# Test map (GARCH 1,1)
test_that('sGARCH_params_to_: check GARCH(1,1) ...', {
  # Original parameters
  params <- c(omega = 0.5, alpha1 = 0.2, beta1 = 0.3)
  # Reparametrized
  zeta <-  sGARCH_params_to_zeta(params, 1, 1, mode = "unitOmega")
  # Map back
  params0 <- sGARCH_params_to_phi(zeta, 1, 1, mode = "unitOmega")
  # Jacobian
  J0 <- sGARCH_params_to_zeta_jacobian(zeta, 1, 1, mode = "unitOmega")
  expect_true(dim(J0)[1] == 3 & dim(J0)[2] == 2)
  # Check equality
  expect_true(sum(abs(params - params0)) <= 10e-10)
})

# Test map (GARCH 1,2)
test_that('sGARCH_params_to_: check GARCH(1,2) ...', {
  # Original parameters
  params <- c(omega = 0.5, alpha1 = 0.2, beta1 = 0.25, beta2 = 0.05)
  # Reparametrized
  zeta <-  sGARCH_params_to_zeta(params, 1, 2, mode = "unitOmega")
  # Map back
  params0 <- sGARCH_params_to_phi(zeta, 1, 2, mode = "unitOmega")
  # Jacobian
  J0 <- sGARCH_params_to_zeta_jacobian(zeta, 1, 2, mode = "unitOmega")
  expect_true(dim(J0)[1] == 4 & dim(J0)[2] == 3)
  # Check equality
  expect_true(sum(abs(params - params0)) <= 10e-10)
})

# Test map (GARCH 1,2)
test_that('sGARCH_params_to_: check GARCH(2,1) ...', {
  # Original parameters
  params <- c(omega = 0.5, alpha1 = 0.25, alpha2 = 0.15, beta1 = 0.1)
  # Reparametrized
  zeta <-  sGARCH_params_to_zeta(params, 2, 1, mode = "unitOmega")
  # Map back
  params0 <- sGARCH_params_to_phi(zeta, 2, 1, mode = "unitOmega")
  # Jacobian
  J0 <- sGARCH_params_to_zeta_jacobian(zeta, 1, 2, mode = "unitOmega")
  expect_true(dim(J0)[1] == 4 & dim(J0)[2] == 3)
  # Check equality
  expect_true(sum(abs(params - params0)) <= 10e-10)
})

# Test map (GARCH 2,2)
test_that('sGARCH_params_to_: check GARCH(2,2) ...', {
  # Original parameters
  params <- c(omega = 0.5, alpha1 = 0.25, alpha2 = 0.15, beta1 = 0.07, beta2 = 0.03)
  # Reparametrized
  zeta <-  sGARCH_params_to_zeta(params, 2, 2, mode = "unitOmega")
  # Map back
  params0 <- sGARCH_params_to_phi(zeta, 2, 2, mode = "unitOmega")
  # Jacobian
  J0 <- sGARCH_params_to_zeta_jacobian(zeta, 2, 2, mode = "unitOmega")
  expect_true(dim(J0)[1] == 5 & dim(J0)[2] == 4)
  # Check equality
  expect_true(sum(abs(params - params0)) <= 10e-10)
})


# Test map (GARCH 0) ERROR
params <- c(omega = 1)
zeta <-  sGARCH_params_to_zeta(params, 0, 0, mode = "unitOmega")
zeta
sGARCH_params_to_phi(zeta, 0, 1, mode = "unitOmega")


