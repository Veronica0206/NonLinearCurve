getBLSGM_Fixed <- function(dat, nT, traj_var, t_var, res_ratio = 4, rho = 0.1, starts = NA, extraTries = NA, paraNames, optimizer = "CSOLNP"){
  mxOption(model = NULL, key = "Default optimizer", optimizer, reset = FALSE)
  if (any(is.na(starts))){
    dat_traj <- dat[, paste0(traj_var, 1:nT)]
    dat_time <- dat[, paste0(t_var, 1:nT)]
    #### Decide the initial value of the intercept mean
    eta0 <- as.numeric(lm(traj_bl ~ ., data = data.frame(traj_bl = dat_traj[, 1]),
                          na.action = na.exclude)$coefficients)
    #### Decide the initial value of the intercept variance
    psi00 <- var(dat_traj[, 1])
    slp <- rep(0, nT - 1)
    for (j in 1:(nT - 1)){
      slp[j] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[, j + 1] - dat_traj[, j])/
                                                                  (time_delta = dat_time[, j + 1] - dat_time[, j])),
                              na.action = na.exclude)$coefficients)
    }
    slp_diff <- diff(slp)
    for (ratio in seq(2, 0, -0.05)){
      knot_candidate <- which(abs(slp_diff) > ratio * abs(mean((dat_traj[, nT] - dat_traj[, 1])/(dat_time[, nT] - dat_time[, 1]), na.rm = T)))
      if (length(knot_candidate) == 2){
        #### Decide the initial value of the knot mean
        gamma <- as.numeric(lm(knot ~ ., data = data.frame(knot = apply(dat_time[, knot_candidate + 1], 1, mean, na.rm = T)),
                               na.action = na.exclude)$coefficients[1])
        break
      }
    }
    for (j in 1:(nT - 1)){
      if (I(gamma > mean(dat_time[, j]) & (gamma < mean(dat_time[, j + 1])))){
        gamma_f <- j
        gamma_c <- j + 1
      }
    }
    #### Decide the initial value of the first slope mean
    eta1 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[, gamma_f] - dat_traj[, 1])/
                                                              (dat_time[, gamma_f] - dat_time[, 1])),
                          na.action = na.exclude)$coefficients[1])
    #### Decide the initial value of the first slope variance
    psi11 <- var((dat_traj[, gamma_f] - dat_traj[, 1])/(dat_time[, gamma_f] - dat_time[, 1]))
    #### Decide the initial value of the second slope mean
    eta2 <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[, nT] - dat_traj[, gamma_c])/
                                                              (dat_time[, nT] - dat_time[, gamma_c])),
                          na.action = na.exclude)$coefficients)
    #### Decide the initial value of the second slope variance
    psi22 <- var((dat_traj[, nT] - dat_traj[, gamma_c])/(dat_time[, nT] - dat_time[, gamma_c]))
    #### Decide the initial value of the residual variance
    residuals <- var(dat_traj[, 1])/res_ratio
    ### Initial values in the original parameter space
    mean0 <- c(eta0, eta1, eta2, gamma)
    psi0 <- matrix(c(psi00, rho * sqrt(psi00 * psi11), rho * sqrt(psi00 * psi22),
                     rho * sqrt(psi00 * psi11), psi11, rho * sqrt(psi11 * psi22),
                     rho * sqrt(psi00 * psi22), rho * sqrt(psi11 * psi22), psi22), nrow = 3)
    ### Transformed matrices obtained by multivariate Delta method
    #### For mean vector
    func0 <- matrix(c(1, mean0[4], 0,
                      0, 0.5, 0.5,
                      0, -0.5, 0.5), nrow = 3, byrow = T)
    #### For var-cov matrix
    grad0 <- matrix(c(1, mean0[4], 0,
                      0, 0.5, 0.5,
                      0, -0.5, 0.5), nrow = 3, byrow = T)

    mean0.s <- func0 %*% mean0[1:3]
    psi0.s <- grad0 %*% psi0 %*% t(grad0)
    starts <- c(mean0.s[1:3], mean0[4], psi0.s[row(psi0.s) >= col(psi0.s)], residuals)
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, 1:nT)
  ### Define latent variables
  latents <- c("eta0s", "eta1s", "eta2s")
  outDef <- list(); outLoads1 <- list(); outLoads2 <- list()
  for(j in 1:nT){
    outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var, j), name = paste0("t", j))
    outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - mug"), name = paste0("L1", j))
    outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - mug)"), name = paste0("L2", j))
  }
  ### Create a mxModel object
  model_mx <- mxModel("Estimate a fixed knot", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loading from latent variables to manifests
                      mxPath(from = "eta0s", to = manifests, arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:nT, "[1,1]")),
                      mxPath(from = "eta2s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 1:nT, "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[11],
                             labels = "residuals"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents, arrows = 1, free = T, values = starts[1:3],
                             labels = c("mueta0s", "mueta1s", "mueta2s")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[5:10],
                             labels = c("psi0s0s", "psi0s1s", "psi0s2s",
                                        "psi1s1s", "psi1s2s", "psi2s2s")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[4],
                               labels = "muknot", name = "mug"),
                      outDef, outLoads1, outLoads2,
                      mxAlgebra(rbind(mueta0s, mueta1s, mueta2s), name = "mean_s"),
                      mxAlgebra(rbind(cbind(psi0s0s, psi0s1s, psi0s2s),
                                      cbind(psi0s1s, psi1s1s, psi1s2s),
                                      cbind(psi0s2s, psi1s2s, psi2s2s)), name = "psi_s"),
                      mxAlgebra(rbind(cbind(1, -mug, mug),
                                      cbind(0, 1, -1),
                                      cbind(0, 1, 1)), name = "func"),
                      mxAlgebra(rbind(cbind(1, -mug, mug),
                                      cbind(0, 1, -1),
                                      cbind(0, 1, 1)), name = "grad"),
                      mxAlgebra(func %*% mean_s, name = "mean"),
                      mxAlgebra(grad %*% psi_s %*% t(grad), name = "psi"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0)
  }
  else{
    model <- mxRun(model_mx)
  }
  model.para <- summary(model)$parameters[, c(1, 5, 6)]
  model.est <- c(model$mean$result, model.para[model.para$name == "muknot", 2],
                 model$psi$result[row(model$psi$result) >= col(model$psi$result)],
                 model.para[model.para$name == "residuals", 2])
  mean.se <- mxSE(mean, model)
  psi.se <- mxSE(psi, model)
  model.se <- c(mean.se, model.para[model.para$name == "muknot", 3],
                psi.se[row(psi.se) >= col(psi.se)],
                model.para[model.para$name == "residuals", 3])
  estimate_out <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
  return(list(model, estimate_out))
}


