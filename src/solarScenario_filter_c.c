#include <R.h>
#include <Rinternals.h>

SEXP solarScenario_filter_c(SEXP z_, SEXP B_, SEXP theta_mix_,
                          SEXP mu1_, SEXP mu2_, SEXP sd1_, SEXP sd2_,
                          SEXP sigma_bar_, SEXP sigma_uncond_,
                          SEXP Yt_tilde_, SEXP eps_, SEXP eps_tilde_,
                          SEXP sigma_, SEXP Yt_tilde_hat_,
                          SEXP omega_, SEXP alpha_, SEXP beta_,
                          SEXP intercept_, SEXP phi_, SEXP ma_theta_,
                          SEXP i_start_) {
  
  // basic sizes
  int n = LENGTH(z_);
  if (LENGTH(B_) != n || LENGTH(theta_mix_) != n ||
      LENGTH(mu1_) != n || LENGTH(mu2_) != n ||
      LENGTH(sd1_) != n || LENGTH(sd2_) != n ||
      LENGTH(sigma_bar_) != n || LENGTH(sigma_uncond_) != n ||
      LENGTH(Yt_tilde_) != n || LENGTH(eps_) != n ||
      LENGTH(eps_tilde_) != n || LENGTH(sigma_) != n ||
      LENGTH(Yt_tilde_hat_) != n) {
    error("All time-series vectors must have the same length.");
  }
  
  int arOrder    = LENGTH(phi_);
  int maOrder    = LENGTH(ma_theta_);
  int archOrder  = LENGTH(alpha_);
  int garchOrder = LENGTH(beta_);
  
  int i_start = INTEGER(i_start_)[0];
  if (i_start < 1 || i_start > n) {
    error("i_start out of range.");
  }
  int istart0 = i_start - 1;  // 0-based index
  
  // ------------------------------------------------------------------
  // Duplicate inputs so we can modify them safely
  // ------------------------------------------------------------------
  SEXP z_out           = PROTECT(duplicate(z_));
  SEXP B_out           = PROTECT(duplicate(B_));
  SEXP theta_mix_out   = PROTECT(duplicate(theta_mix_));
  SEXP mu1_out         = PROTECT(duplicate(mu1_));
  SEXP mu2_out         = PROTECT(duplicate(mu2_));
  SEXP sd1_out         = PROTECT(duplicate(sd1_));
  SEXP sd2_out         = PROTECT(duplicate(sd2_));
  SEXP sigma_bar_out   = PROTECT(duplicate(sigma_bar_));
  SEXP sigma_uncond_out= PROTECT(duplicate(sigma_uncond_));
  SEXP Yt_tilde_out    = PROTECT(duplicate(Yt_tilde_));
  SEXP eps_out         = PROTECT(duplicate(eps_));
  SEXP eps_tilde_out   = PROTECT(duplicate(eps_tilde_));
  SEXP sigma_out       = PROTECT(duplicate(sigma_));
  SEXP Yt_tilde_hat_out= PROTECT(duplicate(Yt_tilde_hat_));
  
  double *z           = REAL(z_out);
  double *B           = REAL(B_out);
  double *theta_mix   = REAL(theta_mix_out);
  double *mu1         = REAL(mu1_out);
  double *mu2         = REAL(mu2_out);
  double *sd1         = REAL(sd1_out);
  double *sd2         = REAL(sd2_out);
  double *sigma_bar   = REAL(sigma_bar_out);
  double *sigma_unc   = REAL(sigma_uncond_out);
  double *Yt_tilde    = REAL(Yt_tilde_out);
  double *eps         = REAL(eps_out);
  double *eps_tilde   = REAL(eps_tilde_out);
  double *sigma       = REAL(sigma_out);
  double *Yt_hat      = REAL(Yt_tilde_hat_out);
  
  double omega      = REAL(omega_)[0];
  double intercept  = REAL(intercept_)[0];
  
  double *alpha = REAL(alpha_);
  double *beta  = REAL(beta_);
  double *phi   = REAL(phi_);
  double *ma_th = REAL(ma_theta_);
  
  // ------------------------------------------------------------------
  // Main time loop
  // ------------------------------------------------------------------
  for (int t = istart0; t < n; t++) {
    
    // 1) GARCH variance recursion
    double sig2 = omega;
    
    for (int k = 0; k < archOrder; k++) {
      int idx = t - 1 - k;
      if (idx < 0) break; // safety guard
      double e = eps_tilde[idx];
      sig2 += alpha[k] * e * e;
    }
    
    for (int k = 0; k < garchOrder; k++) {
      int idx = t - 1 - k;
      if (idx < 0) break;
      double s = sigma[idx];
      sig2 += beta[k] * s * s;
    }
    
    if (sig2 < 0.0) sig2 = 0.0; // numerical guard
    sigma[t] = sqrt(sig2);
    
    // 2) ARMA conditional mean
    double mu = intercept;
    
    for (int k = 0; k < arOrder; k++) {
      int idx = t - 1 - k;
      if (idx < 0) break;
      mu += phi[k] * Yt_tilde[idx];
    }
    
    for (int k = 0; k < maOrder; k++) {
      int idx = t - 1 - k;
      if (idx < 0) break;
      mu += ma_th[k] * eps[idx];
    }
    
    Yt_hat[t] = mu;
    
    // 3) Mixture and residuals
    // shift z by theta_mix
    z[t] += theta_mix[t];
    
    double b  = B[t];
    double zt = z[t];
    
    double u1 = mu1[t] + sd1[t] * zt;
    double u2 = mu2[t] + sd2[t] * zt;
    
    double u_tilde = u1 * b + u2 * (1.0 - b);
    eps_tilde[t]   = sigma[t] * u_tilde;
    eps[t]         = sigma_bar[t] * sigma_unc[t] * eps_tilde[t];
    
    Yt_tilde[t]    = Yt_hat[t] + eps[t];
  }
  
  // ------------------------------------------------------------------
  // Build return list
  // ------------------------------------------------------------------
  SEXP out = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, z_out);
  SET_VECTOR_ELT(out, 1, sigma_out);
  SET_VECTOR_ELT(out, 2, eps_tilde_out);
  SET_VECTOR_ELT(out, 3, eps_out);
  SET_VECTOR_ELT(out, 4, Yt_tilde_out);
  SET_VECTOR_ELT(out, 5, Yt_tilde_hat_out);
  SET_VECTOR_ELT(out, 6, B_out); // B unchanged, but returned for completeness
  
  SEXP names = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names, 0, mkChar("z"));
  SET_STRING_ELT(names, 1, mkChar("sigma"));
  SET_STRING_ELT(names, 2, mkChar("eps_tilde"));
  SET_STRING_ELT(names, 3, mkChar("eps"));
  SET_STRING_ELT(names, 4, mkChar("Yt_tilde"));
  SET_STRING_ELT(names, 5, mkChar("Yt_tilde_hat"));
  SET_STRING_ELT(names, 6, mkChar("B"));
  setAttrib(out, R_NamesSymbol, names);
  
  UNPROTECT(16);
  return out;
}