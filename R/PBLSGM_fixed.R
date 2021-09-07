getPBLSGM_Fixed <- function(dat, nT, traj_var, t_var, res_ratio = rep(4, 2), rho = rep(0.1, 2), btw_rho = rep(0.1, 1), starts = NA, extraTries = NA, paraNames, optimizer = "CSOLNP"){
  mxOption(model = NULL, key = "Default optimizer", optimizer, reset = FALSE)
  if (any(is.na(starts))){
    dat_traj <- dat_time <- list();
    eta0 <- eta1 <- eta2 <- gamma <- psi00 <- psi11 <- psi22 <- uni_residual <- list();
    uni_mean0 <- uni_psi0 <- uni_func0 <- uni_grad0 <- uni_mean0.s <- list()
    for (traj in 1:length(traj_var)){
      dat_traj[[length(dat_traj) + 1]] <- dat[, paste0(traj_var[traj], 1:nT[traj])]
      dat_time[[length(dat_time) + 1]] <- dat[, paste0(t_var[traj], 1:nT[traj])]

      #### Decide the initial value of the intercept mean
      eta0[[length(eta0) + 1]] <- as.numeric(lm(traj_bl ~ ., data = data.frame(traj_bl = dat_traj[[traj]][, 1]),
                                                na.action = na.exclude)$coefficients)
      #### Decide the initial value of the intercept variance
      psi00[[length(psi00) + 1]] <- var(dat_traj[[traj]][, 1])
      slp <- rep(0, nT[traj] - 1)
      for (j in 1:(nT[traj] - 1)){
        slp[j] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[[traj]][, j + 1] - dat_traj[[traj]][, j])/
                                                                    (time_delta = dat_time[[traj]][, j + 1] - dat_time[[traj]][, j])),
                                na.action = na.exclude)$coefficients)
      }
      slp_diff <- diff(slp)
      for (ratio in seq(2, 0, -0.05)){
        knot_candidate <- which(abs(slp_diff) > ratio * abs(mean((dat_traj[[traj]][, nT[traj]] - dat_traj[[traj]][, 1])/
                                                                   (dat_time[[traj]][, nT[traj]] - dat_time[[traj]][, 1]), na.rm = T)))
        if (length(knot_candidate) == 2){
          #### Decide the initial value of the knot mean
          gamma[[length(gamma) + 1]] <- as.numeric(lm(knot ~ ., data = data.frame(knot = apply(dat_time[[traj]][, knot_candidate + 1], 1, mean, na.rm = T)),
                                                      na.action = na.exclude)$coefficients[1])
          break
        }
      }
      for (j in 1:(nT[traj] - 1)){
        if (I(gamma[[traj]] > mean(dat_time[[traj]][, j]) & (gamma[[traj]] < mean(dat_time[[traj]][, j + 1])))){
          gamma_f <- j
          gamma_c <- j + 1
          break
        }
      }
      #### Decide the initial value of the first slope mean
      eta1[[length(eta1) + 1]] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[[traj]][, gamma_f] - dat_traj[[traj]][, 1])/
                                                                                    (dat_time[[traj]][, gamma_f] - dat_time[[traj]][, 1])),
                                                na.action = na.exclude)$coefficients[1])
      #### Decide the initial value of the first slope variance
      psi11[[length(psi11) + 1]] <- var((dat_traj[[traj]][, gamma_f] - dat_traj[[traj]][, 1])/(dat_time[[traj]][, gamma_f] - dat_time[[traj]][, 1]))
      #### Decide the initial value of the second slope mean
      eta2[[length(eta2) + 1]] <- as.numeric(lm(traj_delta ~ ., data = data.frame(traj_delta = (dat_traj[[traj]][, nT[traj]] - dat_traj[[traj]][, gamma_c])/
                                                                                    (dat_time[[traj]][, nT[traj]] - dat_time[[traj]][, gamma_c])),
                                                na.action = na.exclude)$coefficients)
      #### Decide the initial value of the second slope variance
      psi22[[length(psi22) + 1]] <- var((dat_traj[[traj]][, nT[traj]] - dat_traj[[traj]][, gamma_c])/
                                          (dat_time[[traj]][, nT[traj]] - dat_time[[traj]][, gamma_c]))
      #### Decide the initial value of the residual variance
      uni_residual[[length(uni_residual) + 1]] <- var(dat_traj[[traj]][, 1])/res_ratio[traj]
      ### Initial values in the original parameter space
      uni_mean0[[length(uni_mean0) + 1]] <- c(eta0[[traj]], eta1[[traj]], eta2[[traj]], gamma[[traj]])
      uni_psi0[[length(uni_psi0) + 1]] <- matrix(c(psi00[[traj]], rho[traj] * sqrt(psi00[[traj]] * psi11[[traj]]), rho[traj] * sqrt(psi00[[traj]] * psi22[[traj]]),
                                                   rho[traj] * sqrt(psi00[[traj]] * psi11[[traj]]), psi11[[traj]], rho[traj] * sqrt(psi11[[traj]] * psi22[[traj]]),
                                                   rho[traj] * sqrt(psi00[[traj]] * psi22[[traj]]), rho[traj] * sqrt(psi11[[traj]] * psi22[[traj]]), psi22[[traj]]), nrow = 3)
      ### Transformed matrices obtained by multivariate Delta method
      #### For mean vector
      uni_func0[[length(uni_func0) + 1]] <- matrix(c(1, uni_mean0[[traj]][4], 0,
                                                     0, 0.5, 0.5,
                                                     0, -0.5, 0.5), nrow = 3, byrow = T)
      ### For var-cov matrix
      #######################
      uni_grad0[[length(uni_grad0) + 1]] <- matrix(c(1, uni_mean0[[traj]][4], 0,
                                                     0, 0.5, 0.5,
                                                     0, -0.5, 0.5), nrow = 3, byrow = T)
      uni_mean0.s[[length(uni_mean0.s) + 1]] <- c(uni_func0[[traj]] %*% uni_mean0[[traj]][1:3], uni_mean0[[traj]][4])
    }
    psi0 <- as.matrix(Matrix::bdiag(uni_psi0))
    grad0 <- as.matrix(Matrix::bdiag(uni_grad0))
    residuals <- as.matrix(Matrix::bdiag(uni_residual))
    for (traj in 1:(length(traj_var) - 1)){
      psi0[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), (traj * 3 + 1):(traj * 3 + 3)] <-
        matrix(c(btw_rho[traj] * sqrt(psi00[[traj]] * psi00[[traj + 1]]), btw_rho[traj] * sqrt(psi00[[traj]] * psi11[[traj + 1]]), btw_rho[traj] * sqrt(psi00[[traj]] * psi22[[traj + 1]]),
                 btw_rho[traj] * sqrt(psi11[[traj]] * psi00[[traj + 1]]), btw_rho[traj] * sqrt(psi11[[traj]] * psi11[[traj + 1]]), btw_rho[traj] * sqrt(psi11[[traj]] * psi22[[traj + 1]]),
                 btw_rho[traj] * sqrt(psi22[[traj]] * psi00[[traj + 1]]), btw_rho[traj] * sqrt(psi22[[traj]] * psi11[[traj + 1]]), btw_rho[traj] * sqrt(psi22[[traj]] * psi22[[traj + 1]])), nrow = 3)
      psi0[(traj * 3 + 1):(traj * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)] <-
        t(matrix(c(btw_rho[traj] * sqrt(psi00[[traj]] * psi00[[traj + 1]]), btw_rho[traj] * sqrt(psi00[[traj]] * psi11[[traj + 1]]), btw_rho[traj] * sqrt(psi00[[traj]] * psi22[[traj + 1]]),
                   btw_rho[traj] * sqrt(psi11[[traj]] * psi00[[traj + 1]]), btw_rho[traj] * sqrt(psi11[[traj]] * psi11[[traj + 1]]), btw_rho[traj] * sqrt(psi11[[traj]] * psi22[[traj + 1]]),
                   btw_rho[traj] * sqrt(psi22[[traj]] * psi00[[traj + 1]]), btw_rho[traj] * sqrt(psi22[[traj]] * psi11[[traj + 1]]), btw_rho[traj] * sqrt(psi22[[traj]] * psi22[[traj + 1]])), nrow = 3))
      residuals[traj, (traj + 1)] <- residuals[(traj + 1), traj] <- 0.3 * sqrt(uni_residual[[traj]] * uni_residual[[traj + 1]])
    }
    psi0.s <- grad0 %*% psi0 %*% t(grad0)
    starts <- c(unlist(uni_mean0.s), psi0.s[row(psi0.s) >= col(psi0.s)], residuals[row(residuals) >= col(residuals)])
  }
  ### Define manifest variables
  manifests <- c(paste0(traj_var[1], 1:nT[1]), paste0(traj_var[2], 1:nT[2]))
  ### Define latent variables
  latents <- paste0(rep(c("eta0s", "eta1s", "eta2s"), 2), rep(c("Y", "Z"), each = 3))
  outDefY <- outDefZ <- list(); outLoadsY1 <- outLoadsY2 <- outLoadsZ1 <- outLoadsZ2 <- list()
  for(j in 1:nT[1]){
    outDefY[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j), name = paste0("yt", j))
    outLoadsY1[[j]] <- mxAlgebraFromString(paste0("yt", j, " - mugY"), name = paste0("L1", j, "Y"))
    outLoadsY2[[j]] <- mxAlgebraFromString(paste0("abs(yt", j, " - mugY)"), name = paste0("L2", j, "Y"))
  }
  for(j in 1:nT[2]){
    outDefZ[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[2], j), name = paste0("zt", j))
    outLoadsZ1[[j]] <- mxAlgebraFromString(paste0("zt", j, " - mugZ"), name = paste0("L1", j, "Z"))
    outLoadsZ2[[j]] <- mxAlgebraFromString(paste0("abs(zt", j, " - mugZ)"), name = paste0("L2", j, "Z"))
  }
  ### Create a mxModel object
  model_mx <- mxModel("Bivariate Nonlinear Growth Model with Fixed Knots", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loadings from latent variables to manifests
                      mxPath(from = "eta0sY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1sY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:nT[1], "Y", "[1,1]")),
                      mxPath(from = "eta2sY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 1:nT[1], "Y", "[1,1]")),
                      mxPath(from = "eta0sZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1sZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:nT[2], "Z", "[1,1]")),
                      mxPath(from = "eta2sZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 1:nT[2], "Z", "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests[1:nT[1]], to = manifests[1:nT[1]],
                             arrows = 2, free = T, values = starts[30],
                             labels = "residualsY"),
                      mxPath(from = manifests[(nT[1] + 1):(nT[1] + nT[2])], to = manifests[(nT[1] + 1):(nT[1] + nT[2])],
                             arrows = 2, free = T, values = starts[32],
                             labels = "residualsZ"),
                      #### Define the covariances of residuals
                      mxPath(from = manifests[1:nT[1]], to = manifests[(nT[1] + 1):(nT[1] + nT[2])],
                             arrows = 2, free = T, values = starts[31],
                             labels = "residualsYZ"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents, arrows = 1, free = T, values = starts[c(1:3, 5:7)],
                             labels = c("mueta0sY", "mueta1sY", "mueta2sY",
                                        "mueta0sZ", "mueta1sZ", "mueta2sZ")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[c(9:29)],
                             labels = c("psi0sY0sY", "psi0sY1sY", "psi0sY2sY", "psi0sY0sZ", "psi0sY1sZ", "psi0sY2sZ",
                                        "psi1sY1sY", "psi1sY2sY", "psi1sY0sZ", "psi1sY1sZ", "psi1sY2sZ",
                                        "psi2sY2sY", "psi2sY0sZ", "psi2sY1sZ", "psi2sY2sZ",
                                        "psi0sZ0sZ", "psi0sZ1sZ", "psi0sZ2sZ",
                                        "psi1sZ1sZ", "psi1sZ2sZ",
                                        "psi2sZ2sZ")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[4],
                               labels = "muknot_Y", name = "mugY"),
                      mxMatrix("Full", 1, 1, free = T, values = starts[8],
                               labels = "muknot_Z", name = "mugZ"),
                      outDefY, outDefZ, outLoadsY1, outLoadsY2, outLoadsZ1, outLoadsZ2,
                      mxAlgebra(rbind(mueta0sY, mueta1sY, mueta2sY, mueta0sZ, mueta1sZ, mueta2sZ), name = "mean_s"),

                      mxAlgebra(rbind(cbind(psi0sY0sY, psi0sY1sY, psi0sY2sY, psi0sY0sZ, psi0sY1sZ, psi0sY2sZ),
                                      cbind(psi0sY1sY, psi1sY1sY, psi1sY2sY, psi1sY0sZ, psi1sY1sZ, psi1sY2sZ),
                                      cbind(psi0sY2sY, psi1sY2sY, psi2sY2sY, psi2sY0sZ, psi2sY1sZ, psi2sY2sZ),
                                      cbind(psi0sY0sZ, psi1sY0sZ, psi2sY0sZ, psi0sZ0sZ, psi0sZ1sZ, psi0sZ2sZ),
                                      cbind(psi0sY1sZ, psi1sY1sZ, psi2sY1sZ, psi0sZ1sZ, psi1sZ1sZ, psi1sZ2sZ),
                                      cbind(psi0sY2sZ, psi1sY2sZ, psi2sY2sZ, psi0sZ2sZ, psi1sZ2sZ, psi2sZ2sZ)), name = "psi_s"),

                      mxAlgebra(rbind(cbind(1, -mugY, mugY, 0, 0, 0),
                                      cbind(0, 1, -1, 0, 0, 0),
                                      cbind(0, 1, 1, 0, 0, 0),
                                      cbind(0, 0, 0, 1, -mugZ, mugZ),
                                      cbind(0, 0, 0, 0, 1, -1),
                                      cbind(0, 0, 0, 0, 1, 1)), name = "func"),

                      mxAlgebra(rbind(cbind(1, -mugY, mugY, 0, 0, 0),
                                      cbind(0, 1, -1, 0, 0, 0),
                                      cbind(0, 1, 1, 0, 0, 0),
                                      cbind(0, 0, 0, 1, -mugZ, mugZ),
                                      cbind(0, 0, 0, 0, 1, -1),
                                      cbind(0, 0, 0, 0, 1, 1)), name = "grad"),
                      mxAlgebra(func %*% mean_s, name = "mean"),
                      mxAlgebra(grad %*% psi_s %*% t(grad), name = "psi"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0)
  }
  else{
    model <- mxRun(model_mx)
  }
  model.para <- summary(model)$parameters[, c(1, 5, 6)]
  model.est <- c(model$mean$result[1:3], model.para[model.para$name == "muknot_Y", 2],
                 model$mean$result[4:6], model.para[model.para$name == "muknot_Z", 2],
                 model$psi$result[row(model$psi$result) >= col(model$psi$result)],
                 model.para[grep("residuals", model.para$name), 2])
  mean.se <- mxSE(mean, model)
  psi.se <- mxSE(psi, model)
  model.se <- c(mean.se[1:3], model.para[model.para$name == "muknot_Y", 3],
                mean.se[4:6], model.para[model.para$name == "muknot_Z", 3],
                psi.se[row(psi.se) >= col(psi.se)],
                model.para[grep("residuals", model.para$name), 3])
  estimate_out <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
  return(list(model, estimate_out))
}


