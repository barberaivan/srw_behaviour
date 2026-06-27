# The attack occurrence models overestimate the attack probabilities
# because there are follows with zero attacks which p_hat are not zero.
# The models have very good coverage in the response scale, but the 
# means suffer underdispersion.
# Is it because the normal random effect is no good? Whould we need a t
# random effect?
# I will fit a simple random effects model with normal effects to check 
# whether this still happens.

library(rstan)
library(tidyverse); theme_set(theme_bw())
# library(arm)
# library(plyr)      # revalue
# library(markovchain)
library(bayestestR)
# library(mgcv)
# library(viridis)
library(posterior)
library(extraDistr) # rtpois
library(logitnorm)

# Functions ---------------------------------------------------------------

softmax <- function(x) exp(x) / sum(exp(x))

normalize <- function(x) x / sum(x)

hdi_lower <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_low"]])
}

hdi_upper <- function(x) {
  samples <- as.numeric(x)
  return(hdi(samples, ci = 0.95)[["CI_high"]])
}

hdmean <- function(x, ci = 0.95, name = "mu") {
  ci <- hdi(x, ci = ci)
  result <- c(ci$CI_low, mean(x), ci$CI_high)
  names(result) <- paste(rep(name, 3), c("lower", "mean", "upper"), sep = "_")
  result
}

get_map <- function(x) {
  x <- as.numeric(x)
  map <- map_estimate(x)
  return(c(map = map))
}

# Load data -----------------------------------------------------------------

bdata <- read.csv("data/srw_behaviour_data.csv",
                  stringsAsFactors = FALSE, sep = ",")
bdata <- bdata[order(bdata$year, bdata$follow), ]
sumlength <- function(x) c("sum" = sum(x), "length" = length(x))
d <- do.call("data.frame", aggregate(fam ~ follow + date + year, data = bdata, FUN = sumlength))
sum(d$fam.sum == 0) / nrow(d) # half the follows have no attacks.
names(d) <- c("follow", "date", "year", "Nwith", "Nint")
str(d)
d$prop <- d$Nwith / d$Nint
sum(d$Nwith > d$Nint)

dagg <- do.call("data.frame", aggregate(prop ~ year, data = d, FUN = mean))
# M1: normal raneff fol + year  --------------------------------------------
sdata <- list(
  N = nrow(d),
  Nfol = length(unique(d$follow)),
  Nyear = length(unique(d$year)),
  N_with = d$Nwith,
  N_int = d$Nint, 
  follow = as.numeric(factor(d$follow)),
  year = as.numeric(factor(d$year))
)
mcode1 <- stan_model("attack/occurrence models trials.stan", verbose = TRUE)
m1 <- sampling(object = mcode1, data = sdata, # seed = 54325, 
               #chains = 1, cores = 1, iter = 10,
               chains = 4, cores = 4, iter = 2000,  
               control = list(adapt_delta = 0.95, max_treedepth = 15))

p_obs <- as.matrix(m1, "p_obs") %>% t
p_obs <- apply(p_obs, 1, hdmean, name = "p") %>% t %>% as.data.frame
p_obs$p_obs <- d$Nwith / d$Nint
p_obs <- p_obs[order(p_obs$p_mean), ]
p_obs$index <- 1:nrow(d)

ggplot(p_obs, aes(x = index, y = p_mean, ymin = p_lower, ymax = p_upper)) +
  geom_line() + geom_ribbon(alpha = 0.2) + 
  geom_point(data = p_obs, mapping = aes(x = index, y = p_obs),
             inherit.aes = F, alpha = 0.5, col = 2)
# sobreestima los que son puro cero, y subestima los que son altos.
# ESO SE LLAMA SHRINKAGE!!!

sd_marginal <- as.matrix(m1, "marginal_sd")
mu_hat <- as.matrix(m1, "alpha")
mean_hat1 <- numeric(length(mu_hat))
for(i in 1:length(mu_hat)) {
  mean_hat1[i] <- momentsLogitnorm(mu_hat[i], sd_marginal[i])["mean"]
}
hdmean(mean_hat1); mean(d$Nwith / d$Nint)
get_map(mean_hat1)

# M2: student t random effects --------------------------------------------

mcode2 <- stan_model("attack/occurrence models trials_t.stan", verbose = TRUE)
m2 <- sampling(object = mcode2, data = sdata, # seed = 54325, 
               #chains = 1, cores = 1, iter = 10,
               chains = 4, cores = 4, iter = 2000,  
               control = list(adapt_delta = 0.95, max_treedepth = 15))
# this one runs slower, but not so much.
sm2 <- summary(m2)[[1]]

p_obs2 <- as.matrix(m2, "p_obs") %>% t
p_obs2 <- apply(p_obs2, 1, hdmean, name = "p") %>% t %>% as.data.frame
p_obs2$p_obs <- d$Nwith / d$Nint
p_obs2 <- p_obs2[order(p_obs2$p_mean), ]
p_obs2$index <- 1:nrow(d)

ggplot(p_obs2, aes(x = index, y = p_mean, ymin = p_lower, ymax = p_upper)) +
  geom_line() + geom_ribbon(alpha = 0.2) + 
  geom_point(data = p_obs2, mapping = aes(x = index, y = p_obs),
             inherit.aes = F, alpha = 0.5, col = 2)


plot(density(as.numeric(as.matrix(m2, "p_obs")), from = 0, to = 1))
lines(density(d$Nwith / d$Nint, from = 0, to = 1), col = 2, lwd = 2)

# Normal and t random effects estimate the same.

sd_marginal <- as.matrix(m2, "marginal_sd")
mu_hat <- as.matrix(m2, "alpha")
mean_hat2 <- numeric(length(mu_hat))
for(i in 1:length(mu_hat)) {
  mean_hat2[i] <- momentsLogitnorm(mu_hat[i], sd_marginal[i])["mean"]
}
hdmean(mean_hat2); mean(d$Nwith / d$Nint)
get_map(mean_hat2)
mean(dagg$prop)


# M3: add day random effect -----------------------------------------------

sdata3 <- list(
  N = nrow(d),
  Nfol = length(unique(d$follow)),
  Ndate = length(unique(d$date)),
  Nyear = length(unique(d$year)),
  N_with = d$Nwith,
  N_int = d$Nint, 
  follow = as.numeric(factor(d$follow)),
  year = as.numeric(factor(d$year)),
  date = as.numeric(factor(d$date))
)
mcode3 <- stan_model("attack/occurrence models trials_date.stan", verbose = TRUE)
m3 <- sampling(object = mcode3, data = sdata3, # seed = 54325, 
               #chains = 1, cores = 1, iter = 10,
               chains = 4, cores = 4, iter = 2000,  
               control = list(adapt_delta = 0.95, max_treedepth = 15))

p_obs3 <- as.matrix(m3, "p_obs") %>% t
p_obs3 <- apply(p_obs3, 1, hdmean, name = "p") %>% t %>% as.data.frame
p_obs3$p_obs <- d$Nwith / d$Nint
p_obs3 <- p_obs3[order(p_obs3$p_mean), ]
p_obs3$index <- 1:nrow(d)

ggplot(p_obs3, aes(x = index, y = p_mean, ymin = p_lower, ymax = p_upper)) +
  geom_line() + geom_ribbon(alpha = 0.2) + 
  geom_point(data = p_obs3, mapping = aes(x = index, y = p_obs),
             inherit.aes = F, alpha = 0.5, col = 2)

# incluyendo date pasa lo mismo.

p_hat3 <- hdmean(rowMeans(as.matrix(m3, "p_obs")))
mean(d$Nwith / d$Nint)
# a la media "empírica" le emboca re bien. quizás el problema sea con la media 
# marginal de la logitnormal.

sd_marginal <- as.matrix(m3, "marginal_sd")
mu_hat <- as.matrix(m3, "alpha")
mean_hat3 <- numeric(length(mu_hat))
for(i in 1:length(mu_hat)) {
  mean_hat3[i] <- momentsLogitnorm(mu_hat[i], sd_marginal[i])["mean"]
}
hdmean(mean_hat3); mean(d$Nwith / d$Nint)
get_map(mean_hat3)
mean(dagg$prop)


# Conclusion: -------------------------------------------------------------

# Tanto la media empírica como la media marginal de la logitnormal le embocan 
# bien al promedio de todos los datos individuales, pero es igual al promedio 
# de los años previamente promediando follows intraaño. 
# El punto es que la var respuesta es la prop por follow, no por año, 
# y las p de cada follow deberían seguir una logitnormal  con la sd marginal
# de year + follow. 
# Chequear si calculando así las obs se acercan más a lo estimado en el modelo
# de ataque.


