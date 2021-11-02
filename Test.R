## Function Test
library(OpenMx)
devtools::load_all()

### Define Parameter lists
###########################
#### Univariate Development
############################
##### Bilinear spline with a fixed knot
paraFixed <- c("mueta0", "mueta1", "mueta2", "mug",
               paste0("psi", c("00", "01", "02", "11", "12", "22")),
               "residuals")
##### Bilinear spline with a random knot
paraRandom <- c("mueta0", "mueta1", "mueta2", "mug",
                paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")),
                "residuals")
##### Bilinear spline with a fixed knot and two baseline covariates
paraFixedTIC <- c("mueta0", "mueta1", "mueta2", "mug",
                  paste0("psi", c("00", "01", "02", "11", "12", "22")),
                  paste0("beta1", 0:2), paste0("beta2", 0:2),
                  paste0("mux", 1:2), paste0("phi", c("11", "12", "22")),
                  "residuals")
##### Bilinear spline with a random knot and two baseline covariates
paraRandomTIC <- c("mueta0", "mueta1", "mueta2", "mug",
                   paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22",
                                   "2g", "gg")), paste0("beta1", c(0:2, "r")),
                   paste0("beta2", c(0:2, "r")), paste0("mux", 1:2),
                   paste0("phi", c("11", "12", "22")), "residuals")

#### Bivariate Development
##### Bilinear spline with fixed knots
paraBiFixed <- c("mueta0Y", "mueta1Y", "mueta2Y", "mugY", "mueta0Z", "mueta1Z", "mueta2Z", "mugZ",
                 paste0("psi", c("0Y0Y", "0Y1Y", "0Y2Y", "0Y0Z", "0Y1Z", "0Y2Z",
                                 "1Y1Y", "1Y2Y", "1Y0Z", "1Y1Z", "1Y2Z",
                                 "2Y2Y", "2Y0Z", "2Y1Z", "2Y2Z",
                                 "0Z0Z", "0Z1Z", "0Z2Z",
                                 "1Z1Z", "1Z2Z",
                                 "2Z2Z")),
                 "residualsY", "residualsYZ", "residualsZ")

##### Bilinear spline with random knots
paraBiRandom <- c("mueta0Y", "mueta1Y", "mueta2Y", "mugY", "mueta0Z", "mueta1Z", "mueta2Z", "mugZ",
                  paste0("psi", c("0Y0Y", "0Y1Y", "0Y2Y", "0YgY", "0Y0Z", "0Y1Z", "0Y2Z", "0YgZ",
                                  "1Y1Y", "1Y2Y", "1YgY", "1Y0Z", "1Y1Z", "1Y2Z", "1YgZ",
                                  "2Y2Y", "2YgY", "2Y0Z", "2Y1Z", "2Y2Z", "2YgZ",
                                  "gYgY", "gY0Z", "gY1Z", "gY2Z", "gYgZ",
                                  "0Z0Z", "0Z1Z", "0Z2Z", "0ZgZ",
                                  "1Z1Z", "1Z2Z", "1ZgZ",
                                  "2Z2Z", "2ZgZ",
                                  "gZgZ")),
                  "residualsY", "residualsYZ", "residualsZ")

### Test on Generated data sets
################################
#### Paper 1 ###################
################################
data("P1_dat")

BLSGM_F <- getBLSGM_Fixed(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", paraNames = paraFixed)

BLSGM_R <- getBLSGM_Random(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", paraNames = paraRandom)

BLSGM_TIC_F <- getBLSGM_TIC_Fixed(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", x_var = c("x1", "x2"), paraNames = paraFixedTIC)

BLSGM_TIC_R <- getBLSGM_TIC_Random(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", x_var = c("x1", "x2"), paraNames = paraRandomTIC)

################################
#### Paper 2 ###################
################################
data("P2_dat")

BLSGMM_1st <- get_BLSGMM_1st(dat = P2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", paraNames = paraFixed)

BLSGMM_2nd <- get_BLSGMM_2nd(dat = P2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("x1", "x2"),
                             starts = BLSGMM_1st[[3]], beta_starts = matrix(c(0, 0, 0, 0, 1, 1), nrow = 2))

################################
#### Extension 1 ###############
################################
data("E1_dat")

PBLSGM_F <- getPBLSGM_Fixed(dat = E1_dat, nT = rep(10, 2), traj_var = c("Y", "Z"), t_var = rep("T", 2),
                            starts = NA, extraTries = NA, uni_paraNames = paraFixed, paraNames = paraBiFixed)

PBLSGM_R <- getPBLSGM_Random(dat = E1_dat, nT = rep(10, 2), traj_var = c("Y", "Z"), t_var = rep("T", 2),
                             starts = NA, extraTries = NA, uni_paraNames = paraRandom, paraNames = paraBiRandom)

################################
#### Extension 2 ###############
################################
data("E2_dat")

##### Growth Mixture Model
FMM_BLSGM <- getFMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", paraNames = paraFixed)

##### Cluster Predictor Mixture Model
CPMM_BLSGM <- getCPMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("gx1", "gx2"),
                            paraNames = paraFixed)

##### Growth Predictor Mixture Model
GPMM_BLSGM <- getGPMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), paraNames = paraFixedTIC)

##### Full Mixture Model
FullMM_BLSGM <- getFullMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), clus_cov = c("gx1", "gx2"),
                                paraNames = paraFixedTIC)

### Test on real world data set
data("Math_dat")
Math_dat$T1 <- Math_dat$T1 - 60
Math_dat$T2 <- Math_dat$T2 - 60
Math_dat$T3 <- Math_dat$T3 - 60
Math_dat$T4 <- Math_dat$T4 - 60
Math_dat$T5 <- Math_dat$T5 - 60
Math_dat$T6 <- Math_dat$T6 - 60
Math_dat$T7 <- Math_dat$T7 - 60
Math_dat$T8 <- Math_dat$T8 - 60
Math_dat$T9 <- Math_dat$T9 - 60

Math_dat$gx1 <- scale(Math_dat$INCOME)
Math_dat$gx2 <- scale(Math_dat$EDU)
Math_dat$ex1 <- scale(Math_dat$Approach_to_Learning1)
Math_dat$ex2 <- scale(Math_dat$Attention_focus)

################################
#### Paper 1 ###################
################################
BLSGM_F_Math <- getBLSGM_Fixed(dat = Math_dat, nT = 9, traj_var = "M", t_var = "T", paraNames = paraFixed, extraTries = 9)
set.seed(20181022)
BLSGM_R_Math <- getBLSGM_Random(dat = Math_dat, nT = 9, traj_var = "M", t_var = "T", paraNames = paraRandom, extraTries = 29)

BLSGM_TIC_F_Math <- getBLSGM_TIC_Fixed(dat = Math_dat, nT = 9, traj_var = "M", t_var = "T", x_var = c("ex1", "ex2"), paraNames = paraFixedTIC)
set.seed(20181022)
BLSGM_TIC_R_Math <- getBLSGM_TIC_Random(dat = Math_dat, nT = 9, traj_var = "M", t_var = "T", x_var = c("ex1", "ex2"), paraNames = paraRandomTIC, extraTries = 29)

################################
#### Paper 2 ###################
################################
set.seed(20181022)
BLSGMM_1st_Math <- get_BLSGMM_1st(dat = Math_dat, nT = 9, nClass = 3, traj_var = "M", t_var = "T",
                                  res_ratio = c(3, 3, 3), rho = rep(0.1, 3), prop_starts = c(0.23, 0.48, 0.29),
                                  paraNames = paraFixed, extraTries = 9)

BLSGMM_2nd_Math <- get_BLSGMM_2nd(dat = Math_dat, nT = 9, nClass = 3, traj_var = "M", t_var = "T",
                                  clus_cov = c("gx1", "gx2"), starts = BLSGMM_1st_Math[[3]],
                                  beta_starts = matrix(c(0, 0, 0, 0.3, 1, 1, 0.2, 1, 1), nrow = 3))


################################
#### Extension 2 ###############
################################
##### Growth Mixture Model
set.seed(20181022)
FMM_BLSGM_Math <- getFMM_BLSGM(dat = Math_dat, nT = 9, nClass = 3, traj_var = "M", t_var = "T",
                               res_ratio = c(3, 3, 3), rho = rep(0.1, 3), prop_starts = c(0.18, 0.45, 0.37),
                               paraNames = paraFixed, loc = 1, scale = 0.05, extraTries = 9)

##### Cluster Predictor Mixture Model
set.seed(20181022)
CPMM_BLSGM_Math <- getCPMM_BLSGM(dat = Math_dat, nT = 9, nClass = 3, traj_var = "M", t_var = "T", clus_cov = c("gx1", "gx2"),
                                 res_ratio = c(3, 3, 3), rho = rep(0.1, 3), prop_starts = c(0.18, 0.45, 0.37),
                                 paraNames = paraFixed, loc = 1, scale = 0.05, extraTries = 9)

##### Growth Predictor Mixture Model
set.seed(20181022)
GPMM_BLSGM_Math <- getGPMM_BLSGM(dat = Math_dat, nT = 9, nClass = 3, traj_var = "M", t_var = "T", growth_cov = c("ex1", "ex2"),
                                 res_ratio = c(3, 3, 3), rho = rep(0.1, 3), prop_starts = c(0.18, 0.45, 0.37),
                                 paraNames = paraFixedTIC, loc = 1, scale = 0.05, extraTries = 9)

##### Full Mixture Model
set.seed(20181022)
FullMM_BLSGM_Math <- getFullMM_BLSGM(dat = Math_dat, nT = 9, nClass = 3, traj_var = "M", t_var = "T", growth_cov = c("ex1", "ex2"), clus_cov = c("gx1", "gx2"),
                                     res_ratio = c(3, 3, 3), rho = rep(0.1, 3), prop_starts = c(0.18, 0.45, 0.37),
                                     paraNames = paraFixedTIC, loc = 1, scale = 0.05, extraTries = 9)










