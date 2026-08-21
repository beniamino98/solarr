

mom1 <- model_invgumb$moments$conditional[10,]
pdf_Y1 <- function(x) dmixnorm(x, c(mom1$M_Y0, mom1$M_Y1), c(mom1$S_Y0, mom1$S_Y1), c(mom1$p1, 1-mom1$p1))
cdf_Y1 <- function(x) pmixnorm(x, c(mom1$M_Y0, mom1$M_Y1), c(mom1$S_Y0, mom1$S_Y1), c(mom1$p1, 1-mom1$p1))
mom2 <- model_gumb$moments$conditional[10,]
pdf_Y2 <- function(x) dmixnorm(x, c(mom2$M_Y0, mom2$M_Y1), c(mom1$S_Y0, mom2$S_Y1), c(mom2$p1, 1-mom2$p1))
cdf_Y2 <- function(x) pmixnorm(x, c(mom2$M_Y0, mom2$M_Y1), c(mom1$S_Y0, mom2$S_Y1), c(mom2$p1, 1-mom2$p1))
mom3 <- model_logis$moments$conditional[10,]
pdf_Y3 <- function(x) dmixnorm(x, c(mom3$M_Y0, mom3$M_Y1), c(mom3$S_Y0, mom3$S_Y1), c(mom3$p1, 1-mom3$p1))
cdf_Y3 <- function(x) pmixnorm(x, c(mom3$M_Y0, mom3$M_Y1), c(mom3$S_Y0, mom3$S_Y1), c(mom3$p1, 1-mom3$p1))
mom4 <- model_norm$moments$conditional[10,]
pdf_Y4<- function(x) dmixnorm(x, c(mom4$M_Y0, mom4$M_Y1), c(mom4$S_Y0, mom4$S_Y1), c(mom4$p1, 1-mom4$p1))
cdf_Y4<- function(x) pmixnorm(x, c(mom4$M_Y0, mom4$M_Y1), c(mom4$S_Y0, mom4$S_Y1), c(mom4$p1, 1-mom4$p1))
Ct <- mom1$Ct
alpha <- mom1$alpha
beta <- mom1$beta

grid <- mom1$Ct*(1-seq(mom1$alpha, mom1$alpha+mom1$beta, 0.01))

ggplot()+
  geom_line(aes(grid, dsolarGHI(grid, Ct, alpha, beta, pdf_Y1, link = "invgumbel")), color = "blue")+
  geom_line(aes(grid, dsolarGHI(grid, Ct, alpha, beta, pdf_Y2, link = "gumbel")), color = "green")+
  geom_line(aes(grid, dsolarGHI(grid, Ct, alpha, beta, pdf_Y3, link = "logis")))+
  geom_line(aes(grid, dsolarGHI(grid, Ct, alpha, beta, pdf_Y4, link = "norm")), color = "red")

cdf <- ecdf(filter(model$data, Day == mom1$Day & Month == mom1$Month)$GHI)

ggplot()+
  geom_line(aes(grid, psolarGHI(grid, Ct, alpha, beta, cdf_Y1, link = "invgumbel")), color = "blue")+
  #geom_line(aes(grid, psolarGHI(grid, Ct, alpha, beta, cdf_Y2, link = "gumbel")), color = "green")+
  geom_line(aes(grid, psolarGHI(grid, Ct, alpha, beta, cdf_Y3, link = "logis")))+
  geom_line(aes(grid, psolarGHI(grid, Ct, alpha, beta, cdf_Y4, link = "norm")), color = "red")+
  geom_line(aes(grid, cdf(grid)), color = "purple")

qsolarGHI(0.05, Ct, alpha, beta, cdf_Y1, link = "invgumbel")
qsolarGHI(0.05, Ct, alpha, beta, cdf_Y2, link = "gumbel")
qsolarGHI(0.05, Ct, alpha, beta, cdf_Y3, link = "logis")
qsolarGHI(0.05, Ct, alpha, beta, cdf_Y4, link = "norm")

