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

