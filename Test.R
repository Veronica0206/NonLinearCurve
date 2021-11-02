## Function Test
library(OpenMx)
devtools::load_all()

### Test on Generated data sets
################################
#### Paper 1 ###################
################################
data("P1_dat")
paraFixed <- c("mueta0", "mueta1", "mueta2", "mug",
               paste0("psi", c("00", "01", "02", "11", "12", "22")),
               "residuals")
BLSGM_F <- getBLSGM_Fixed(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", paraNames = paraFixed)
paraRandom <- c("mueta0", "mueta1", "mueta2", "mug",
                paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22", "2g", "gg")),
                "residuals")
BLSGM_R <- getBLSGM_Random(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", paraNames = paraRandom)
paraFixedTIC <- c("mueta0", "mueta1", "mueta2", "mug",
                  paste0("psi", c("00", "01", "02", "11", "12", "22")),
                  paste0("beta1", 0:2), paste0("beta2", 0:2),
                  paste0("mux", 1:2), paste0("phi", c("11", "12", "22")),
                  "residuals")
BLSGM_TIC_F <- getBLSGM_TIC_Fixed(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", x_var = c("x1", "x2"), paraNames = paraFixedTIC)
paraRandomTIC <- c("mueta0", "mueta1", "mueta2", "mug",
                   paste0("psi", c("00", "01", "02", "0g", "11", "12", "1g", "22",
                                   "2g", "gg")), paste0("beta1", c(0:2, "r")),
                   paste0("beta2", c(0:2, "r")), paste0("mux", 1:2),
                   paste0("phi", c("11", "12", "22")), "residuals")
BLSGM_TIC_R <- getBLSGM_TIC_Random(dat = P1_dat, nT = 10, traj_var = "Y", t_var = "T", x_var = c("x1", "x2"), paraNames = paraRandomTIC)

################################
#### Paper 2 ###################
################################
data("E2_dat")
BLSGMM_1st <- get_BLSGMM_1st(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", paraNames = paraFixed)
BLSGMM_2nd <- get_BLSGMM_2nd(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("gx1", "gx2"), starts = BLSGMM_1st[[3]], beta_starts = c(0, 1, 1))

################################
#### Extension 1 ###############
################################
data("E1_dat")
paraBiFixed <- c("mueta0Y", "mueta1Y", "mueta2Y", "mugY", "mueta0Z", "mueta1Z", "mueta2Z", "mugZ",
                 paste0("psi", c("0Y0Y", "0Y1Y", "0Y2Y", "0Y0Z", "0Y1Z", "0Y2Z",
                                 "1Y1Y", "1Y2Y", "1Y0Z", "1Y1Z", "1Y2Z",
                                 "2Y2Y", "2Y0Z", "2Y1Z", "2Y2Z",
                                 "0Z0Z", "0Z1Z", "0Z2Z",
                                 "1Z1Z", "1Z2Z",
                                 "2Z2Z")),
                 "residualsY", "residualsYZ", "residualsZ")
PBLSGM_F <- getPBLSGM_Fixed(dat = E1_dat, nT = rep(10, 2), traj_var = c("Y", "Z"), t_var = rep("T", 2), starts = NA, extraTries = NA, paraNames = paraBiFixed)
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
PBLSGM_R <- getPBLSGM_Random(dat = E1_dat, nT = rep(10, 2), traj_var = c("Y", "Z"), t_var = rep("T", 2), starts = NA, extraTries = NA, paraNames = paraBiRandom)

################################
#### Extension 2 ###############
################################
data("E2_dat")
##### Growth Mixture Model
paraFixed <- c("mueta0", "mueta1", "mueta2", "mug",
               paste0("psi", c("00", "01", "02", "11", "12", "22")),
               "residuals")
paraFixedTIC <- c("mueta0", "mueta1", "mueta2", "mug",
                  paste0("psi", c("00", "01", "02", "11", "12", "22")),
                  paste0("beta1", 0:2), paste0("beta2", 0:2),
                  paste0("mux", 1:2), paste0("phi", c("11", "12", "22")),
                  "residuals")
FMM_BLSGM <- getFMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", paraNames = paraFixed)

##### Cluster Predictor Mixture Model
CPMM_BLSGM <- getCPMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("gx1", "gx2"),
                            beta_starts = c(0, 1, 1), paraNames = paraFixed)

##### Growth Predictor Mixture Model
GPMM_BLSGM <- getGPMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), paraNames = paraFixedTIC)

##### Full Mixture Model
FullMM_BLSGM <- getFullMM_BLSGM(dat = E2_dat, nT = 10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), clus_cov = c("gx1", "gx2"),
                                beta_starts = c(0, 1, 1), paraNames = paraFixedTIC)

### Test on real world data set
################################
#### Paper 1 ###################
################################
data("Math")
Math$T1 <- Math$T1 - 60
Math$T2 <- Math$T2 - 60
Math$T3 <- Math$T3 - 60
Math$T4 <- Math$T4 - 60
Math$T5 <- Math$T5 - 60
Math$T6 <- Math$T6 - 60
Math$T7 <- Math$T7 - 60
Math$T8 <- Math$T8 - 60
Math$T9 <- Math$T9 - 60


BLSGM_F_math <- getBLSGM_Fixed(dat = Math, nT = 9, traj_var = "M", t_var = "T", extraTries = 29)
set.seed(20181022)
BLSGM_R_math <- getBLSGM_Random(dat = Math, nT = 9, traj_var = "M", t_var = "T", extraTries = 29)

Math$x1 <- scale(Math$Approach_to_Learning1)
Math$x2 <- scale(Math$Attention_focus)

BLSGM_TIC_F_math <- getBLSGM_TIC_Fixed(dat = Math, nT = 9, traj_var = "M", t_var = "T", x_var = c("x1", "x2"), extraTries = 29)
set.seed(20181022)
BLSGM_TIC_R_math <- getBLSGM_TIC_Random(dat = Math, nT = 9, traj_var = "M", t_var = "T", x_var = c("x1", "x2"), extraTries = 499)
















