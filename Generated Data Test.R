## Function Test
library(OpenMx)
mxOption(model = NULL, key = "Default optimizer", "CSOLNP", reset = FALSE)
devtools::load_all()
# NPSOL | SLSQP | CSOLNP
##########################
#### BLSGM ###############
##########################
### Test on Generated data sets

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
paraBiFixed <- c("mueta0Y", "mueta1Y", "mueta2Y", "mugY", paste0("psi", c("0Y0Y", "0Y1Y", "0Y2Y", "1Y1Y", "1Y2Y", "2Y2Y")),
                 "mueta0Z", "mueta1Z", "mueta2Z", "mugZ", paste0("psi", c("0Z0Z", "0Z1Z", "0Z2Z", "1Z1Z", "1Z2Z", "2Z2Z")),
                 paste0("psi", c("0Y0Z", "1Y0Z", "2Y0Z", "0Y1Z", "1Y1Z", "2Y1Z", "0Y2Z", "1Y2Z", "2Y2Z")),
                 "residualsY", "residualsYZ", "residualsZ")

##### Bilinear spline with random knots
paraBiRandom <- c("mueta0Y", "mueta1Y", "mueta2Y", "mugY", paste0("psi", c("0Y0Y", "0Y1Y", "0Y2Y", "0YgY", "1Y1Y", "1Y2Y", "1YgY", "2Y2Y", "2YgY", "gYgY")),
                  "mueta0Z", "mueta1Z", "mueta2Z", "mugZ", paste0("psi", c("0Z0Z", "0Z1Z", "0Z2Z", "0ZgZ", "1Z1Z", "1Z2Z", "1ZgZ", "2Z2Z", "2ZgZ", "gZgZ")),
                  paste0("psi", c("0Y0Z", "1Y0Z", "2Y0Z", "gY0Z",
                                  "0Y1Z", "1Y1Z", "2Y1Z", "gY1Z",
                                  "0Y2Z", "1Y2Z", "2Y2Z", "gY2Z",
                                  "0YgZ", "1YgZ", "2YgZ", "gYgZ")),
                  "residualsY", "residualsYZ", "residualsZ")

#### Longitudinal Mediation
paraMed2 <- c("muetaY1", "muetaYr", "muetaY2", "mugY", "muetaM1", "muetaMr", "muetaM2", "mugM",
              paste0("psi", c("Y1Y1", "Y1Yr", "Y1Y2",
                              "YrYr", "YrY2",
                              "Y2Y2",
                              "M1M1", "M1Mr", "M1M2",
                              "MrMr", "MrM2",
                              "M2M2")),
              "muX", "phi11", paste0("Mediator", c("11", "1r", "rr", "12", "r2", "22")),
              paste0("total", c("1", "r", "2")),
              paste0("beta", rep(c("Y", "M"), each = 3), rep(c(1, "r", 2), 2)),
              paste0("beta", c("M1Y1", "M1Yr", "M1Y2", "MrYr", "MrY2", "M2Y2")),
              "residualsY", "residualsM", "residualsYM")

paraMed3 <- c("muetaY1", "muetaYr", "muetaY2", "mugY", "muetaM1", "muetaMr", "muetaM2", "mugM", "muetaX1", "muetaXr", "muetaX2", "mugX",
              paste0("psi", c("Y1Y1", "Y1Yr", "Y1Y2",
                              "YrYr", "YrY2",
                              "Y2Y2",
                              "M1M1", "M1Mr", "M1M2",
                              "MrMr", "MrM2",
                              "M2M2",
                              "X1X1", "X1Xr", "X1X2",
                              "XrXr", "XrX2",
                              "X2X2")),
              paste0("beta", c("X1Y1", "X1Yr", "X1Y2", "XrYr", "XrY2", "X2Y2",
                               "X1M1", "X1Mr", "X1M2", "XrMr", "XrM2", "X2M2",
                               "M1Y1", "M1Yr", "M1Y2", "MrYr", "MrY2", "M2Y2")),
              paste0("mediator", c("111", "rrr", "222", "11r", "1rr", "112", "1r2", "122", "rr2", "r22")),
              paste0("total", c("11", "1r", "12", "rr", "r2", "22")),
              "residualsY", "residualsM", "residualsX", "residualsYM", "residualsYX", "residualsMX")

################################
#### Paper 1 ###################
################################
data("BLSGM_uni_dat")

BLSGM_F <- getBLSGM_Fixed(dat = BLSGM_uni_dat, T_records = 1:10, traj_var = "Y", t_var = "T", paraNames = paraFixed)

BLSGM_R <- getBLSGM_Random(dat = BLSGM_uni_dat, T_records = 1:10, traj_var = "Y", t_var = "T", paraNames = paraRandom)

BLSGM_TIC_F <- getBLSGM_TIC_Fixed(dat = BLSGM_uni_dat, T_records = 1:10, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), paraNames = paraFixedTIC)

BLSGM_TIC_R <- getBLSGM_TIC_Random(dat = BLSGM_uni_dat, T_records = 1:10, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), paraNames = paraRandomTIC)

################################
#### Paper 2 ###################
################################
data("BLSGM_uni_sub_dat")

BLSGMM_1st <- get_BLSGMM_1st(dat = BLSGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", paraNames = paraFixed)

BLSGMM_2nd <- get_BLSGMM_2nd(dat = BLSGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("gx1", "gx2"),
                             starts = BLSGMM_1st[[2]], beta_starts = matrix(c(0, 0, 0, 0, 1, 1), nrow = 2))

################################
#### Extension 1 ###############
################################
data("BLSGM_bi_dat")

PBLSGM_F <- getPBLSGM_Fixed(dat = BLSGM_bi_dat, T_records = list(1:10, 1:10), traj_var = c("Y", "Z"), t_var = rep("T", 2), paraNames = paraBiFixed)

PBLSGM_R <- getPBLSGM_Random(dat = BLSGM_bi_dat, T_records = list(1:10, 1:10), traj_var = c("Y", "Z"), t_var = rep("T", 2), paraNames = paraBiRandom)

################################
#### Extension 2 ###############
################################
##### Growth Mixture Model
FMM_BLSGM <- getFMM_BLSGM(dat = BLSGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", paraNames = paraFixed)

##### Cluster Predictor Mixture Model
CPMM_BLSGM <- getCPMM_BLSGM(dat = BLSGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("gx1", "gx2"),
                            paraNames = paraFixed)

##### Growth Predictor Mixture Model
GPMM_BLSGM <- getGPMM_BLSGM(dat = BLSGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"),
                            paraNames = paraFixedTIC)

##### Full Mixture Model
FullMM_BLSGM <- getFullMM_BLSGM(dat = BLSGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"),
                                clus_cov = c("gx1", "gx2"), paraNames = paraFixedTIC)

################################
#### Extension 3 ###############
################################
data("BLSGM_bi_sub_dat")
FMM_PBLSGM <- getFMM_PBLSGM(dat = BLSGM_bi_sub_dat, T_records = list(1:10, 1:10), nClass = 2, traj_var = c("Y", "Z"), t_var = rep("T", 2), paraNames = paraBiFixed)

CPMM_PBLSGM <- getCPMM_PBLSGM(dat = BLSGM_bi_sub_dat, T_records = list(1:10, 1:10), nClass = 2, traj_var = c("Y", "Z"), t_var = rep("T", 2), clus_cov = c("x1", "x2"),
                              paraNames = paraBiFixed)

################################
#### Extension 4 ###############
################################
data("BLSGM_med2_dat")

Med2_BLSGM <- getBLSGM_Med2(dat = BLSGM_med2_dat, M_records = 1:10, Y_records = 1:10, X_var = "X", M_var = "M", Y_var = "Y",
                            t_var = rep("T", 2), paraNames = paraMed2)

data("BLSGM_med3_dat")

Med3_BLSGM <- getBLSGM_Med3(dat = BLSGM_med3_dat, X_records = 1:10, M_records = 1:10, Y_records = 1:10, X_var = "X", M_var = "M", Y_var = "Y",
                            t_var = rep("T", 3), paraNames = paraMed3)

#########################
#### LBGM ###############
#########################
### Define Parameter lists
###########################
#### Univariate Development
############################
##### LBGM
paraLBGM <- c("mueta0", "mueta1", paste0("psi", c("00", "01", "11")), paste0("rel_rate", 1:9), paste0("abs_rate", 1:9),
              "residuals")

paraLBGM_LCSM <- c("mueta0", "mueta1", paste0("psi", c("00", "01", "11")),
                   paste0("rel_rate", 1:9), paste0("abs_rate", 1:9), paste0("abs_rate_se", 1:9),
                   paste0("change_in_interval", 1:9), paste0("change_from_baseline", 1:9), "residuals")

paraLBGMTIC <- c("mueta0", "mueta1", paste0("psi", c("00", "01", "11")), paste0("beta1", 0:1), paste0("beta2", 0:1),
                 paste0("mux", 1:2), paste0("phi", c("11", "12", "22")), paste0("rel_rate", 1:9), paste0("abs_rate", 1:9),
                 "residuals")

paraBiLBGM <- c("mueta0Y", "mueta1Y", paste0("psi", c("0Y0Y", "0Y1Y", "1Y1Y")), paste0("rel_rateY", 1:9), paste0("rel_rateZ", 1:9),
                "mueta0Z", "mueta1Z", paste0("psi", c("0Z0Z", "0Z1Z", "1Z1Z")), paste0("abs_rateY", 1:9), paste0("abs_rateZ", 1:9),
                paste0("psi", c("0Y0Z", "1Y0Z", "0Y1Z", "1Y1Z")),
                "residualsY", "residualsYZ", "residualsZ")

data("LBGM_uni_dat")

LBGM <- getLBGM(dat = LBGM_uni_dat, T_records = 1:10, traj_var = "Y", t_var = "T", paraNames = paraLBGM_LCSM)

LBGM_TIC <- getLBGM_TIC(dat = LBGM_uni_dat, T_records = 1:10, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"), paraNames = paraLBGMTIC)

data("LBGM_uni_sub_dat")
FMM_LBGM <- getFMM_LBGM(dat = LBGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T",  paraNames = paraLBGM)

CPMM_LBGM <- getCPMM_LBGM(dat = LBGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", clus_cov = c("gx1", "gx2"),
                          paraNames = paraLBGM)

GPMM_LBGM <- getGPMM_LBGM(dat = LBGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"),
                          paraNames = paraLBGMTIC)

FullMM_LBGM <- getFullMM_LBGM(dat = LBGM_uni_sub_dat, T_records = 1:10, nClass = 2, traj_var = "Y", t_var = "T", growth_cov = c("ex1", "ex2"),
                              clus_cov = c("gx1", "gx2"), paraNames = paraLBGMTIC)

data("LBGM_bi_dat")

PLBGM <- getPLBGM(dat = LBGM_bi_dat, T_records = list(1:10, 1:10), traj_var = c("Y", "Z"), t_var = rep("T", 2), paraNames = paraBiLBGM)

#########################
#### QUAD ###############
#########################
paraQUAD_LGCM <- c("mueta0", "mueta1", "mueta2", paste0("psi", c("00", "01", "02", "11", "12", "22")), "residuals")

paraQUAD_LCSM <- c("mueta0", "mueta1", "mueta2", paste0("psi", c("00", "01", "02", "11", "12", "22")), paste0("instant_rate_est", 1:9),
                   paste0("instant_rate_var", 1:9), paste0("change_in_interval", 1:9), paste0("change_from_baseline", 1:9), "residuals")


data("QUAD_dat")

QUAD_LGCM <- getLGCM_QUAD(dat = QUAD_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                          paraNames = paraQUAD_LGCM)
QUAD_LCSM <- getLCSM_QUAD(dat = QUAD_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                          paraNames = paraQUAD_LCSM)

#########################
#### EXP ################
#########################
paraEXP_LGCM <- c("mueta0", "mueta1", paste0("psi", c("00", "01", "11")), "gamma", "residuals")

paraEXP_LCSM <- c("mueta0", "mueta1", paste0("psi", c("00", "01", "11")), "gamma", paste0("instant_rate_est", 1:9),
                  paste0("instant_rate_var", 1:9), paste0("change_in_interval", 1:9), paste0("change_from_baseline", 1:9), "residuals")

data("EXP_dat")

EXP_LGCM <- getLGCM_EXP(dat = EXP_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                        paraNames = paraEXP_LGCM)
EXP_LCSM <- getLCSM_EXP(dat = EXP_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                        paraNames = paraEXP_LCSM)

#########################
#### JB #################
#########################
paraJB_LGCM <- c("mueta0", "mueta1", "mueta2", paste0("psi", c("00", "01", "02", "11", "12", "22")),
                 "gamma", "residuals")

paraJB_LCSM <- c("mueta0", "mueta1", "mueta2", paste0("psi", c("00", "01", "02", "11", "12", "22")), "gamma", paste0("instant_rate_est", 1:9),
                 paste0("instant_rate_var", 1:9), paste0("change_in_interval", 1:9), paste0("change_from_baseline", 1:9), "residuals")


data("JB_dat")

JB_LGCM <- getLGCM_JB(dat = JB_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                      paraNames = paraJB_LGCM)
JB_LCSM <- getLCSM_JB(dat = JB_dat, T_records = 1:10, traj_var = "Y", t_var = "T",
                      paraNames = paraJB_LCSM)

save.image("../../test/simu_dat.Rdata")
