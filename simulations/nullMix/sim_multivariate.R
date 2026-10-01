rm(list = ls())

library(lme4)
library(flip)
library(lmerTest)
library(multiwayvcov)
library(lmtest)
library(jointest)

source("simulations/nullMix/createCrossedData.R")
source("simulations/nullMix/nullMix.R")
source("simulations/nullMix/mix_perm.R")

nsim = 1000
R = 20
meanC = 10
sA = 5
sB = 5
sE = 5
L = 10
beta = c(rep(0, L*0.8), rep(1, L*0.2))
B = 1000
alpha = 0.05
seed = 123
rho = 0 #In the paper, rho takes values in {0, 0.2, 0.4, 0.6, 0.8}

set.seed(seed)

sim <- expand.grid(nsim = seq(nsim), R = R, sE = sE,
                   meanC = meanC)

simRC <- expand.grid(R = R, meanC = meanC, KEEP.OUT.ATTRS = FALSE)
Clist <- vector("list", nrow(simRC))
for (j in seq_len(nrow(simRC))) {
  Cj <- stats::rpois(simRC$R[j], simRC$meanC[j])
  Cj <- pmax(Cj, 2L)
  Clist[[j]] <- Cj
}

cor_sim <- c()

pv.eigen_rot <- matrix(NA, ncol = L, nrow = nsim)
pv.eigen_rot1 <- matrix(NA, ncol = L, nrow = nsim)
pv.eigen_rot2 <- matrix(NA, ncol = L, nrow = nsim)

pv.eigen_flip <- matrix(NA, ncol = L, nrow = nsim)
pv.eigen_flip1 <- matrix(NA, ncol = L, nrow = nsim)
pv.eigen_flip2 <- matrix(NA, ncol = L, nrow = nsim)

pv.eigen_perm <- matrix(NA, ncol = L, nrow = nsim)
pv.eigen_perm1 <- matrix(NA, ncol = L, nrow = nsim)
pv.eigen_perm2 <- matrix(NA, ncol = L, nrow = nsim)

pv.lmer <- matrix(NA, ncol = L, nrow= nsim)

observed_outcome_corr <- function(y, use = "pairwise.complete.obs") {
  R <- cor(y, use = use)
  mean_pair <- mean(R[upper.tri(R)])
  list(R = R, mean_pairwise = mean_pair)
}

for(i in seq(nrow(sim))){
  message("simulation ", i, " / ", nrow(sim))

  idx <- which(simRC$R == sim$R[i] & simRC$meanC == sim$meanC[i])
  C <- Clist[[idx]]

  outDB <- createCrossedData(sA = sA, sB = sB, sE = sim$sE[i],
                             R = sim$R[i], C = C,
                             beta = beta,
                             rho = rho,
                             seed = seed + i)

  #rotation
  eigen_mom <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = FALSE,
    variance_method = "mom",
    beta0 = rep(0,L),
    transformation = "rotation"
  )

  eigen_lmer <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = FALSE,
    variance_method = "lmer",
    beta0 = rep(0,L),
    transformation = "rotation"
  )

  eigen_known <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = TRUE,
    sA = sA,
    sB = sB,
    sE = sim$sE[i],
    beta0 = rep(0,L),
    transformation = "rotation"
  )

  pv.eigen_rot[i,] <- eigen_mom$pv
  pv.eigen_rot1[i,] <- eigen_lmer$pv
  pv.eigen_rot2[i,] <- eigen_known$pv

  #signflip
  eigen_mom <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = FALSE,
    variance_method = "mom",
    beta0 = rep(0,L),
    transformation = "signflip"
  )

  eigen_lmer <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = FALSE,
    variance_method = "lmer",
    beta0 = rep(0,L),
    transformation = "signflip"
  )


  eigen_known <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = TRUE,
    sA = sA,
    sB = sB,
    sE = sim$sE[i],
    beta0 = rep(0,L),
    transformation = "signflip"
  )

  pv.eigen_flip[i,] <- eigen_mom$pv
  pv.eigen_flip1[i,] <- eigen_lmer$pv
  pv.eigen_flip2[i,] <- eigen_known$pv

  #permutation
  eigen_mom <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = FALSE,
    variance_method = "mom",
    beta0 = rep(0,L),
    transformation = "permutation"
  )

  eigen_lmer <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = FALSE,
    variance_method = "lmer",
    beta0 = rep(0,L),
    transformation = "permutation"
  )

  eigen_known <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = paste0("y", seq(L)),
    known = TRUE,
    sA = sA,
    sB = sB,
    sE = sim$sE[i],
    beta0 = rep(0,L),
    transformation = "permutation"
  )

  pv.eigen_perm[i,] <- eigen_mom$pv
  pv.eigen_perm1[i,] <- eigen_lmer$pv
  pv.eigen_perm2[i,] <- eigen_known$pv

  #lmer and method ettore

  out_lmer <- sapply(seq(L), function(x){

    idx_y <- which(colnames(outDB$db) == paste0("y", x))
    outDB$db$y <- outDB$db[[idx_y]]

    mod <- lmer(y ~ x + (1|subjects) + (1|items), data = outDB$db)
    sum_mod <- summary(mod)
    sum_mod$coefficients[2,5]
  })

  pv.lmer[i,] <- stats::p.adjust(out_lmer, method = "holm")

  colsY <- grep("^y", names(outDB$db), value = TRUE)
  dbY   <- outDB$db[, colsY]

  cor_sim[i] <- observed_outcome_corr(y = dbY)$mean_pairwise
}


