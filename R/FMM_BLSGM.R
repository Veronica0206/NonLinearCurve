getFMM_BLSGM <- function(dat, nT, nClass, traj_var, t_var, res_ratio = 4, rho = 0.1, starts = NA, w_starts = 1, extraTries = NA, paraNames, optimizer = "CSOLNP"){
  mxOption(model = NULL, key = "Default optimizer", optimizer, reset = FALSE)
  if (any(is.na(starts))){
    starts <- list()
    ID <- 1:nrow(dat)
    dat_order <- dat[ID[order(dat[, paste0(traj_var, 1)])], ]
    for (k in 1:nClass){
      dat_traj <- dat_order[as.integer(nrow(dat)/nClass * (k - 1) + 1):as.integer(nrow(dat)/nClass * k), paste0(traj_var, 1:nT)]
      dat_time <- dat_order[as.integer(nrow(dat)/nClass * (k - 1) + 1):as.integer(nrow(dat)/nClass * k), paste0(t_var, 1:nT)]
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
        if (length(knot_candidate) == 4){
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
      starts[[length(starts) + 1]] <- c(mean0.s[1:3], mean0[4], psi0.s[row(psi0.s) >= col(psi0.s)], residuals)
    }
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, 1:nT)
  ### Define latent variables
  latents <- c("eta0s", "eta1s", "eta2s")
  outDef <- list(); outLoads1 <- list(); outLoads2 <- list()
  class.list <- list()
  for (k in 1:nClass){
    for(j in 1:nT){
      outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.T", j),
                              name = paste0("t", j))
      outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " -c", k, "mug"), name = paste0("c", k, "L1", j))
      outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " -c", k, "mug)"),
                                            name = paste0("c", k, "L2", j))
    }
    ### Create a mxModel object
    class.list[[k]] <- mxModel(name = paste0("Class", k), type = "RAM",
                               manifestVars = manifests, latentVars = latents,
                               #### Define factor loadings from latent variables to manifests
                               mxPath(from = "eta0s", to = manifests, arrows = 1, free = F, values = 1),
                               mxPath(from = "eta1s", to = manifests, arrows = 1, free = F, values = 0,
                                      labels = paste0("c", k, "L1", 1:nT, "[1,1]")),
                               mxPath(from = "eta2s", to = manifests, arrows = 1, free = F, values = 0,
                                      labels = paste0("c", k, "L2", 1:nT, "[1,1]")),
                               #### Define the variances of residuals
                               mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[[k]][11],
                                      labels = paste0("c", k, "residuals")),
                               #### Define means of latent variables
                               mxPath(from = "one", to = latents[1:3], arrows = 1, free = T, values = starts[[k]][1:3],
                                      labels = paste0("c", k, c("mueta0s", "mueta1s", "mueta2s"))),
                               #### Define var-cov matrix of latent variables
                               mxPath(from = latents, to = latents, arrows = 2,
                                      connect = "unique.pairs", free = T,
                                      values = starts[[k]][c(5:10)],
                                      labels = paste0("c", k, c("psi0s0s", "psi0s1s", "psi0s2s",
                                                                "psi1s1s", "psi1s2s", "psi2s2s"))),
                               #### Add additional parameter and constraints
                               mxMatrix("Full", 1, 1, free = T, values = starts[[k]][4],
                                        labels = paste0("c", k, "muknot"), name = paste0("c", k, "mug")),
                               outDef, outLoads1, outLoads2,
                               mxAlgebraFromString(paste0("rbind(c", k, "mueta0s, c", k, "mueta1s, c", k, "mueta2s)"),
                                                   name = paste0("c", k, "mean_s")),
                               mxAlgebraFromString(paste0("rbind(cbind(c", k, "psi0s0s, c", k, "psi0s1s, c", k, "psi0s2s),",
                                                          "cbind(c", k, "psi0s1s, c", k, "psi1s1s, c", k, "psi1s2s),",
                                                          "cbind(c", k, "psi0s2s, c", k, "psi1s2s, c", k, "psi2s2s))"),
                                                   name = paste0("c", k, "psi_s")),
                               mxAlgebraFromString(paste0("rbind(cbind(", "1,", "-c", k, "muknot,", "c", k, "muknot),",
                                                          "cbind(0, 1, -1), cbind(0, 1, 1))"),
                                                   name = paste0("c", k, "func")),
                               mxAlgebraFromString(paste0("rbind(cbind(", "1,", "-c", k, "muknot,", "c", k, "muknot),",
                                                          "cbind(0, 1, -1), cbind(0, 1, 1))"),
                                                   name = paste0("c", k, "grad")),
                               mxAlgebraFromString(paste0("c", k, "func %*% c", k, "mean_s"), name = paste0("c", k, "mean")),
                               mxAlgebraFromString(paste0("c", k, "grad %*% c", k, "psi_s %*% t(c", k, "grad)"),
                                                   name = paste0("c", k, "psi")),
                               mxFitFunctionML(vector = T))
  }
  ### Make the class proportion matrix, fixing one parameter at a non-zero constant (one)
  classP <- mxMatrix("Full", nClass, 1, free = c(F, rep(T, nClass - 1)), values = w_starts,
                     labels = paste0("w", 1:nClass), name = "weights")
  algebraObjective <- mxExpectationMixture(paste0("Class", 1:nClass), weights = "weights", scale = "softmax")
  objective <- mxFitFunctionML()
  model_mx <- mxModel("Finite Mixture Model, BLSGM with Fixed Knots",
                      mxData(observed = dat, type = "raw"), class.list, classP, algebraObjective, objective)
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0)
  }
  else{
    model <- mxRun(model_mx)
  }
  model.para <- summary(model)$parameters[, c(1, 5, 6)]
  ### Estimates of parameters in each latent class
  model.est <- model.se <- est <- list()
  for (k in 1:nClass){
    model.est[[k]] <- c(mxEvalByName(paste0("c", k, "mean"), model = model@submodels[[k]]),
                        model.para[model.para$name == paste0("c", k, "muknot"), 2],
                        mxEvalByName(paste0("c", k, "psi"), model = model@submodels[[k]])[
                          row(mxEvalByName(paste0("c", k, "psi"), model = model@submodels[[k]])) >=
                            col(mxEvalByName(paste0("c", k, "psi"), model = model@submodels[[k]]))],
                        model.para[model.para$name == paste0("c", k, "residuals"), 2])
    model.se[[k]] <- c(mxSE(paste0("Class", k, ".c", k, "mean"), model, forceName = T),
                       model.para[model.para$name == paste0("c", k, "muknot"), 3],
                       mxSE(paste0("Class", k, ".c", k, "psi"), model, forceName = T)[
                         row(mxSE(paste0("Class", k, ".c", k, "psi"), model, forceName = T)) >=
                           col(mxSE(paste0("Class", k, ".c", k, "psi"), model, forceName = T))],
                       model.para[model.para$name == paste0("c", k, "residuals"), 3])
    est[[k]] <- data.frame(Name = paste0("c", k, paraNames),
                           Estimate = model.est[[k]], SE = model.se[[k]])
  }
  est.weights <- data.frame(Name = paste0("p", 2:nClass),
                            Estimate = model.para[grep("w", model.para$name), 2],
                            SE = model.para[grep("w", model.para$name), 3])
  estimate_out <- rbind(do.call(rbind.data.frame, est), est.weights)
  return(list(model, estimate_out))
}
