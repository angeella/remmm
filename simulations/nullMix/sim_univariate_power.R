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
R = 10 #in the paper, R takes values in the set {10, 20, 30}
meanC = 5
sA = 5
sB = 5
sE = 5
beta = 1
B = 1000
alpha = 0.05
seed = 123

set.seed(seed)

sim <- expand.grid(nsim = seq(1000),R = R,sE = sE,meanC = meanC)


simRC <- expand.grid(R =R, meanC = meanC)

Clist <- list()
for(j in seq(nrow(simRC))){
  C <- rpois(simRC$R[j], simRC$meanC[j])
  C <- ifelse(C == 0, 2, C)
  Clist[[j]] <- C
}


for(i in seq(nrow(sim))){
  message("simulation ", i, " / ", nrow(sim))

  idx <- which(simRC$R == sim$R[i] & simRC$meanC == sim$meanC[i])
  C <- Clist[[idx]]

  outDB <- createCrossedData(sA = sA, sB = sB, sE = sim$sE[i],
                             R = sim$R[i], C = C, beta = beta, seed = seed + i)

  #rotation
  eigen_mom <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = "y",
    known = FALSE,
    variance_method = "mom",
    beta0 =  0,
    transformation = "rotation"
  )

  eigen_lmer <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = "y",
    known = FALSE,
    variance_method = "lmer",
    beta0 = 0,
    transformation = "rotation"
  )

  eigen_known <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = "y",
    known = TRUE,
    sA = sA,
    sB = sB,
    sE = sim$sE[i],
    beta0 = 0,
    transformation = "rotation"
  )

  sim$pv.eigen_rot[i] <- eigen_mom$pv
  sim$pv.eigen_rot1[i] <- eigen_lmer$pv
  sim$pv.eigen_rot2[i] <- eigen_known$pv

  #signflip
  eigen_mom <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = "y",
    known = FALSE,
    variance_method = "mom",
    beta0 = 0,
    transformation = "signflip"
  )

  eigen_lmer <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = "y",
    known = FALSE,
    variance_method = "lmer",
    beta0 = 0,
    transformation = "signflip"
  )


  eigen_known <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = "y",
    known = TRUE,
    sA = sA,
    sB = sB,
    sE = sim$sE[i],
    beta0 = 0,
    transformation = "signflip"
  )

  sim$pv.eigen_flip[i] <- eigen_mom$pv
  sim$pv.eigen_flip1[i] <- eigen_lmer$pv
  sim$pv.eigen_flip2[i] <- eigen_known$pv

  #permutation
  eigen_mom <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = "y",
    known = FALSE,
    variance_method = "mom",
    beta0 = 0,
    transformation = "permutation"
  )

  eigen_lmer <- nullMix(
    db = outDB$db,
    B = 1000,
    Zg = outDB$Zg,
    y_cols = "y",
    known = FALSE,
    variance_method = "lmer",
    beta0 = 0,
    transformation = "permutation"
  )

  eigen_known <- nullMix(
    db = outDB$db,
    B = B,
    Zg = outDB$Zg,
    y_cols = "y",
    known = TRUE,
    sA = sA,
    sB = sB,
    sE = sim$sE[i],
    beta0 = 0,
    transformation = "permutation"
  )

  sim$pv.eigen_perm[i] <- eigen_mom$pv
  sim$pv.eigen_perm1[i] <- eigen_lmer$pv
  sim$pv.eigen_perm2[i] <- eigen_known$pv

  #lmer and method ettore

  sim[i, c("pv.lmer", "pv.mix_perm")] = tryCatch({

    out1 <- lmer(y ~ x + (1| subjects) + (1 | items),
                 data = outDB$db)

    a = summary(out1)

    if(!is.null(a$optinfo$conv$lme4$messages)){
      sim$mod.lmer.war1[i] = a$optinfo$conv$lme4$messages
    }

    re <- ranef(out1)
    out_res <- mix_perm(db = outDB$db,
                        mod = out1,
                        rsubj_int = as.vector(re$subjects[,1]),
                        rsubj_slop = rep(0, length(as.vector(re$subjects[,1]))),
                        rstim_int = as.vector(re$items[,1]),
                        rstim_slop = rep(0, length(as.vector(re$items[,1]))))

    c(a$coefficients[2,c(5)], out_res)



  }, error=function(e) NA)




  #LM with Cameron 2011 correction multi-way

  sim$pv.lm[i] = tryCatch({
    out <- lm(y ~ x, data = outDB$db)

    vcov_out <- cluster.vcov(out, cluster = cbind(outDB$db$subjects, outDB$db$items))
    a = lmtest::coeftest(out, vcov_out)
    a[2,4]
  }, error=function(e) NA)


}

