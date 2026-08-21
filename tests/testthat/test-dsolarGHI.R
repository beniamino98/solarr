library(ggplot2)

alpha <- 0.05
beta <- 0.93
Ct <- 8
df <- 5
grid <- seq(Ct*(1-alpha-beta), Ct*(1-alpha), length.out = 100)[-c(100)]
sd <- c(1, 0.75, 0.5, 0.25, 0.1)


# ********************************************* Normal *********************************************
# Density
pdf1 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dnorm(x, sd = sd[1]))
pdf2 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dnorm(x, sd = sd[2]))
pdf3 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dnorm(x, sd = sd[3]))
pdf4 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dnorm(x, sd = sd[4]))
pdf5 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dnorm(x, sd = sd[5]))
# Compute expected values
integrate(function(x) x*pdf1(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf2(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf3(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf4(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf5(x), lower = min(grid), upper = max(grid))$value

ggplot()+
  geom_line(aes(grid, pdf1(grid), color = "sd1"))+
  geom_line(aes(grid, pdf2(grid), color = "sd2"))+
  geom_line(aes(grid, pdf3(grid), color = "sd3"))+
  geom_line(aes(grid, pdf4(grid), color = "sd4"))+
  geom_line(aes(grid, pdf5(grid), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Solar radiation", y = "Density", title = "Normal")
# Density
cdf1 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pnorm(x, sd = sd[1]))
cdf2 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pnorm(x, sd = sd[2]))
cdf3 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pnorm(x, sd = sd[3]))
cdf4 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pnorm(x, sd = sd[4]))
cdf5 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pnorm(x, sd = sd[5]))
# Plot
ggplot()+
  geom_line(aes(grid, cdf1(grid), color = "sd1"))+
  geom_line(aes(grid, cdf2(grid), color = "sd2"))+
  geom_line(aes(grid, cdf3(grid), color = "sd3"))+
  geom_line(aes(grid, cdf4(grid), color = "sd4"))+
  geom_line(aes(grid, cdf5(grid), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Solar radiation", y = "Distribution", title = "Normal")



# ********************************************* Student t *********************************************
pdf1 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dt(x/sd[1], df = df))
pdf2 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dt(x/sd[2], df = df))
pdf3 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dt(x/sd[3], df = df))
pdf4 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dt(x/sd[4], df = df))
pdf5 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dt(x/sd[5], df = df))
# Compute expected values
integrate(function(x) x*pdf1(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf2(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf3(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf4(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf5(x), lower = min(grid), upper = max(grid))$value
# Plot
ggplot()+
  geom_line(aes(grid, pdf1(grid), color = "sd1"))+
  geom_line(aes(grid, pdf2(grid), color = "sd2"))+
  geom_line(aes(grid, pdf3(grid), color = "sd3"))+
  geom_line(aes(grid, pdf4(grid), color = "sd4"))+
  geom_line(aes(grid, pdf5(grid), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Solar radiation", y = "Density", title = "Student-t")
# Distribution
cdf1 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pt(x/sd[1], df = df))
cdf2 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pt(x/sd[2], df = df))
cdf3 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pt(x/sd[3], df = df))
cdf4 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pt(x/sd[4], df = df))
cdf5 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pt(x/sd[5], df = df))
# Plot
ggplot()+
  geom_line(aes(grid, cdf1(grid), color = "sd1"))+
  geom_line(aes(grid, cdf2(grid), color = "sd2"))+
  geom_line(aes(grid, cdf3(grid), color = "sd3"))+
  geom_line(aes(grid, cdf4(grid), color = "sd4"))+
  geom_line(aes(grid, cdf5(grid), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Solar radiation", y = "Distribution", title = "Student-t")

# ********************************************* Gaussian mixture *********************************************
# Density
pdf1 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[1]/2, sd[1])))
pdf2 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[2]/2, sd[2])))
pdf3 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[3]/2, sd[3])))
pdf4 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[4]/2, sd[4])))
pdf5 <- function(x) dsolarGHI(x, Ct, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[5]/2, sd[5])))
# Compute expected values
integrate(function(x) x*pdf1(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf2(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf3(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf4(x), lower = min(grid), upper = max(grid))$value
integrate(function(x) x*pdf5(x), lower = min(grid), upper = max(grid))$value
# Plot
ggplot()+
  geom_line(aes(grid, pdf1(grid), color = "sd1"))+
  geom_line(aes(grid, pdf2(grid), color = "sd2"))+
  geom_line(aes(grid, pdf3(grid), color = "sd3"))+
  geom_line(aes(grid, pdf4(grid), color = "sd4"))+
  geom_line(aes(grid, pdf5(grid), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Solar radiation", y = "Density", title = "Gaussian Mixture (-0.1, 0.1, sd/2, sd)")

# Distribution
cdf1 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[1]/2, sd[1])))
cdf2 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[2]/2, sd[2])))
cdf3 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[3]/2, sd[3])))
cdf4 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[4]/2, sd[4])))
cdf5 <- function(x) psolarGHI(x, Ct, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[5]/2, sd[5])))
# Plot
ggplot()+
  geom_line(aes(grid, cdf1(grid), color = "sd1"))+
  geom_line(aes(grid, cdf2(grid), color = "sd2"))+
  geom_line(aes(grid, cdf3(grid), color = "sd3"))+
  geom_line(aes(grid, cdf4(grid), color = "sd4"))+
  geom_line(aes(grid, cdf5(grid), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Solar radiation", y = "Distribution", title = "Gaussian Mixture (-0.1, 0.1, sd/2, sd)")





