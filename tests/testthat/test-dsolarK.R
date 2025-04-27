library(solarr)
library(tidyverse)

# Bounds parameters
alpha <- 0.05
beta <- 0.93
# Grid of points
grid <- seq(1-alpha-beta, 1-alpha, length.out = 100)[-c(100)]
# Mixture standard deviations
sd <- c(1, 0.75, 0.5, 0.25, 0.1)


ggplot()+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = sd[1])), color = "sd1"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = sd[2])), color = "sd2"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = sd[3])), color = "sd3"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = sd[4])), color = "sd4"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dnorm(x, sd = sd[5])), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Clearness index", y = "Density", title = "Normal")

# Mixture
ggplot()+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[1]/2, sd[1]))), color = "sd1"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[2]/2, sd[2]))), color = "sd2"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[3]/2, sd[3]))), color = "sd3"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[4]/2, sd[4]))), color = "sd4"))+
  geom_line(aes(grid, dsolarK(grid, alpha, beta, function(x) dmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[5]/2, sd[5]))), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Clearness index", y = "Density", title = "Gaussian Mixture (-0.1, 0.1, sd/2, sd)")

# Distribution function Kt
ggplot()+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = sd[1])), color = "sd1"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = sd[2])), color = "sd2"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = sd[3])), color = "sd3"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = sd[4])), color = "sd4"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pnorm(x, sd = sd[5])), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Clearness index", y = "Density", title = "Normal")

# Mixture
ggplot()+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[1]/2, sd[1]))), color = "sd1"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[2]/2, sd[2]))), color = "sd2"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[3]/2, sd[3]))), color = "sd3"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[4]/2, sd[4]))), color = "sd4"))+
  geom_line(aes(grid, psolarK(grid, alpha, beta, function(x) pmixnorm(x, mean = c(-0.1, 0.1), sd = c(sd[5]/2, sd[5]))), color = "sd5"))+
  scale_color_manual(values = c(sd1 = "red", sd2 = "blue", sd3 = "magenta", sd4 = "green", sd5 = "black"),
                     labels = c(sd1 = sd[1], sd2 = sd[2], sd3 = sd[3], sd4 = sd[4], sd5 = sd[5]))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "Std. Dev", x = "Clearness index", y = "Density", title = "Gaussian Mixture (-0.1, 0.1, sd/2, sd)")

