library("magrittr")
library("lme4")
library("future.apply")
library("dplyr")
library("ggplot2")

sim_data <- function(N, J) {

  b <- MASS::mvrnorm(N, mu = c(0, 0),
                     Sigma = matrix(c(1.5, -0.35, -0.35, 0.2), nrow = 2))
  
  expand.grid(id = 1:N,
                    time = seq(0, 10, length = J)) %>%
    mutate(time = jitter(time, factor = 0.5),
           y = 30 + 0.76 * time + b[id, 1] + 
             b[id, 2] * time + rnorm(N*J, 0, 1))
}



# Vary N - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
J = 5
Ns = c(2e6, 1e6) %>% unique


times_N <- lapply(rev(Ns), function(N) {
  cat("Running J = ", J, ", N = ", N, "\n")
  
  dat <- sim_data(N = N, J = J)
  
  ft <- system.time({
    fit <- lmer(y ~ time + (time | id), data = dat)
  })["elapsed"]
  
  conv = fit@optinfo$conv$lme4
  if (!is.null(conv$code))
    print(conv$messages)
  
  data.frame(N = N, J = J, nrd = 2, sec = ft, 
             code = ifelse(is.null(conv$code), 0, conv$code)) %>% print()
}) %>%
  do.call(rbind, .)

write.csv(times_N, "lmer_times_N.csv", row.names = FALSE)

ggplot(times_N, aes(x = J*N, y = sec)) +
  geom_point(aes(color = factor(code)))


# Vary J - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

N = 1000
Js = c(9000) %>% unique

times_J <- lapply(rev(Js), function(J) {
  cat("Running J = ", J, ", N = ", N, "\n")
  dat <- sim_data(N = N, J = J)
  
  ft <- system.time({
    fit <- lmer(y ~ time + (time | id), data = dat)
  })["elapsed"]
  
  conv = fit@optinfo$conv$lme4
  if (!is.null(conv$code))
    print(conv$messages)
  
  data.frame(N = N, J = J, nrd = 2, sec = ft, 
             code = ifelse(is.null(conv$code), 0, conv$code))
}) %>%
  do.call(rbind, .)

write.csv(times_J, "lmer_times_J.csv", row.names = FALSE)

ggplot(times_J, aes(x = J*N, y = sec)) +
  geom_point(aes(color = factor(code)))

rbind(
  mutate(times_N, run = "N"),
  mutate(times_J, run = "J")) %>%
  ggplot(aes(x = N*J, y = sec, color = run)) +
  geom_point()



# Random effects - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

N = 50
J = 100


dat <- sim_data(N = N, J = J)

times_rd2 <- lapply(c(12, 13, 14), function (i) {
  
  cat("Running J = ", J, ", N = ", N, " with ", i + 1, " random effects.\n")
  
  ft <- system.time({
    fit <- lmer(y ~ poly(time, i) + (poly(time, i) | id), data = dat)
  })["elapsed"]

  conv = fit@optinfo$conv$lme4
  if (!is.null(conv$code))
    print(conv$messages)
  
  data.frame(N = N, J = J, nrd = i + 1, sec = ft, 
             code = ifelse(is.null(conv$code), 0, conv$code)) %>%
    print()
}) %>%
  do.call(rbind, .)

write.csv(rbind(times_rd1, times_rd2), "lmer_times_rd.csv", row.names = FALSE)

rbind(times_rd2, times_rd1) %>%
  mutate(nobs = N*J,
         beta = 2,
         var = nrd/2 + (nrd^2)/2 + 1,
         b = nrd*N)
ggplot(aes(x = nrd, y = sec)) +
  geom_point(aes(color = factor(N)))



# many covariates

sim_dataX <- function(N, J, p) {
  
  b <- MASS::mvrnorm(N, mu = c(0, 0),
                     Sigma = matrix(c(1.5, -0.35, -0.35, 0.2), nrow = 2))
  
  X <- MASS::mvrnorm(N * J, 
                     mu = rep(0, p),
                     Sigma = diag(p))
  
  beta <- runif(p, -2, 2)
  
  expand.grid(id = 1:N,
              time = seq(0, 10, length = J)) %>%
    cbind(., X) %>%
    mutate(time = jitter(time, factor = 0.5),
           y = -2 + X %*% beta +  0.76 * time + b[id, 1] + 
             b[id, 2] * time + rnorm(N*J, 0, 1))
}


N = 2000
J = 10

times_p = lapply(rev(c(800, 900, 950)), function(p) {
  dat = sim_dataX(N = N, J = J, p = p)
  print(dim(dat))
  
  ft <- system.time({
    fit <- lmer(y ~ . -id + (time | id), data = dat)
  })["elapsed"]
  
  conv = fit@optinfo$conv$lme4
  if (!is.null(conv$code))
    print(conv$messages)
  
  data.frame(N = N, J = J, nrd = 2, p = p, sec = ft, 
             code = ifelse(is.null(conv$code), 0, conv$code))
}) %>%
  do.call(rbind, .)

write.csv(times_p2, "lmer_times_p.csv", row.names = FALSE)

ggplot(times_p, aes(x = p, y = sec)) +
  geom_point(aes(color = factor(code)))




# multivariate models -------------------------------------------------
library("MCMCglmm")


sim_dataMV <- function(N, J, K) {
  
  r = 2
  S = matrix(nrow = r * K, ncol = r * K, data = 0.35)
  diag(S) = 1
  
  betas = runif(r * K, -2, 2)
  print(betas)
  
  b <- MASS::mvrnorm(N, mu = rep(0, r * K),
                     Sigma = S)
  
  d = expand.grid(id = 1:N,
              time = seq(0, 10, length = J)) %>%
    mutate(time = jitter(time, factor = 0.5))
  
  sigs = abs(rnorm(K, 0, 2))
  cat("sigmas: ", sigs, "\n")
  
  for (k in 1:K) {
    d[[paste0("y", k)]] = betas[r * k - 1] + betas[r*k] * d$time + 
      b[d$id, r * k - 1] + b[d$id, r*k] * d$time + rnorm(N*J, 0, sigs[k])
  }
  d
}

d = sim_dataMV(N = 100, J = 10, K = 3)



dlong = reshape2::melt(d, id.vars = c("id", "time")) %>%
  mutate(.by = "variable",
         weight = 1/var(value))


library(nlme)
test1 <- lme(value ~ time*variable, random = ~ time*variable | id, 
             data = dlong, weights = varIdent(form = ~ 1 | variable))
test11 = lme(y1 ~ time, random = ~ time | id, data = d)
test12 = lme(y2 ~ time, random = ~ time | id, data = d)
test13 = lme(y3 ~ time, random = ~ time | id, data = d)

test2 = lmer(value ~ variable * time + (variable * time | id), data = dlong)
test2a = lmer(value ~ variable * time + (variable * time | id), data = dlong,
              weights = weight)

test21 = lmer(y1 ~ time + (time | id), data = d)
test22 = lmer(y2 ~ time + (time | id), data = d)
test23 = lmer(y3 ~ time + (time | id), data = d)


ndf = d
ndflong = dlong
ndflong$fit1 = predict(test1)
ndflong$fit2 = predict(test2)
ndflong$fit2a = predict(test2a)

ndf$fit21 = predict(test21)
ndf$fit22 = predict(test22)
ndf$fit23 = predict(test23)



theid = 5

ggplot(subset(ndflong, id == theid), aes(x = time, y = value, group = variable)) +
  geom_point() +
  geom_line(aes(y = fit2), color = "red") +
  geom_line(aes(y = fit1), color = "green4") +
  geom_line(aes(y = fit2a), color = "magenta") +
  facet_wrap(~ variable) +
  geom_line(data = mutate(subset(ndf, id == theid), variable = "y1"), 
            aes(y = fit21), color = "blue", linetype = 2) +
  geom_line(data = mutate(subset(ndf, id == theid), variable = "y2"), 
            aes(y = fit22), color = "blue", linetype = 2) +
  geom_line(data = mutate(subset(ndf, id == theid), variable = "y3"), 
            aes(y = fit23), color = "blue", linetype = 2)



summary(test1)

library("brms")


# Fit the bivariate model
t0 = Sys.time()
model <- brm(
  bf(y1 ~ time + (1 + time | id)) +
    bf(y2 ~ time + (1 + time | id)) +
    bf(y3 ~ time + (1 + time | id)),
  data = d,
  chains = 3,
  cores = 3,
  iter = 1000
)
t1 = Sys.time()


library(JMbayes)
t0a = Sys.time()
mjmb = mvglmer(
  formulas = list(
  y1 ~ time + (1 + time | id),
  y2 ~ time + (1 + time | id),
  y3 ~ time + (1 + time | id)),
  data = d,
  families = list("gaussian", "gaussian", "gaussian"),
  control = list(
    n.iter = 1000,
    n.burnin = 0,
    # n.chains = 3,
    n.adapt = 5,
    n.thin = 1,
    seed = NULL
))
t1a = Sys.time()

plot(mjmb$mcmc$betas1[, 1])
plot(mjmb)
