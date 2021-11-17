getPBLSGM_Random <- function(dat, nT, traj_var, t_var, res_ratio = rep(4, 2), rho = rep(0.1, 2), btw_rho = rep(0.1, 1), btw_res = rep(0.3, 1),
                             rho_gamma = list(c(rep(0, 3)), c(rep(0, 3))), btw_rho_gamma = list(c(rep(0, 3))), btw_gamma = rep(0, 1), starts = NA,
                             loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA, optimizer = "CSOLNP"){
  mxOption(model = NULL, key = "Default optimizer", optimizer, reset = FALSE)
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  if (any(is.na(starts))){
    mug <- psigg <- rep(0, length(traj_var))
    uni_residual <- list();
    uni_mean0.s <- uni_psi0 <- uni_grad0 <- list()
    ### Obtain estimates from the reduced model
    reduced_model <- getPBLSGM_Fixed(dat = dat, nT = nT, traj_var = traj_var, t_var = t_var, res_ratio = res_ratio, rho = rho,
                                     btw_rho = btw_rho, btw_res = btw_res, original = F)
    mug <- summary(reduced_model)$parameters$Estimate[c(grep("muknot", summary(reduced_model)$parameters$name))]
    reduced_psi0 <- reduced_model$psi$result
    for (traj in 1:length(traj_var)){
      ### Mean vector (in the reparameterized framework)
      uni_mean0.s[[length(uni_mean0.s) + 1]] <- c(reduced_model$mean_s$result[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)], mug[traj])
      dat_time <- dat[, paste0(t_var, 1:nT[traj])]
      mean_time <- apply(dat_time, 2, mean)
      for (t in 1:length(mean_time)){
        if (mean_time[t] <= mug[traj] & mean_time[t + 1] >= mug[traj]){
          psigg[traj] <- var(apply(dat_time[, c(t, t + 1)], 1, mean))
          stop
        }
      }
      ### var-cov matrix (in the original framework)
      ###################
      uni_psi0[[length(uni_psi0) + 1]] <- matrix(0, nrow = 4, ncol = 4)
      uni_psi0[[length(uni_psi0)]][1:3, 1:3] <- reduced_psi0[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3), ((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)]
      uni_psi0[[length(uni_psi0)]][4, 4] <- psigg[traj]
      uni_psi0[[length(uni_psi0)]][1, 4] <- uni_psi0[[length(uni_psi0)]][4, 1] <-
        rho_gamma[[traj]][1] * sqrt(psigg[traj] * uni_psi0[[length(uni_psi0)]][1, 1])
      uni_psi0[[length(uni_psi0)]][2, 4] <- uni_psi0[[length(uni_psi0)]][4, 2] <-
        rho_gamma[[traj]][2] * sqrt(psigg[traj] * uni_psi0[[length(uni_psi0)]][2, 2])
      uni_psi0[[length(uni_psi0)]][3, 4] <- uni_psi0[[length(uni_psi0)]][4, 3] <-
        rho_gamma[[traj]][3] * sqrt(psigg[traj] * uni_psi0[[length(uni_psi0)]][3, 3])
      uni_grad0[[length(uni_grad0) + 1]] <- matrix(c(1, mug[traj], 0, reduced_model$mean$result[((traj - 1) * 3 + 2)],
                                                     0, 0.5, 0.5, 0,
                                                     0, -0.5, 0.5, 0,
                                                     0, 0, 0, 1), nrow = 4, byrow = T)
      uni_residual[[length(uni_residual) + 1]] <- reduced_model$S$values[(traj - 1) * nT[traj] + 1, (traj - 1) * nT[traj] + 1]
    }
  }
  psi0 <- as.matrix(Matrix::bdiag(uni_psi0))
  grad0 <- as.matrix(Matrix::bdiag(uni_grad0))
  residuals <- as.matrix(Matrix::bdiag(uni_residual))
  for (traj_i in 1:(length(traj_var) - 1)){
    for (traj_j in traj_i:(length(traj_var) - 1)){
      reduced_cor <- reduced_psi0[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)]
      psi0[((traj_i - 1) * 4 + 1):((traj_i - 1) * 4 + 4), (traj_j * 4 + 1):(traj_j * 4 + 4)] <-
        matrix(c(reduced_cor[1, ], btw_rho_gamma[[traj_i]][1] * sqrt(psi0[(traj_i - 1) * 4 + 1, (traj_i - 1) * 4 + 1] * psigg[traj_j + 1]),
                 reduced_cor[2, ], btw_rho_gamma[[traj_i]][2] * sqrt(psi0[(traj_i - 1) * 4 + 2, (traj_i - 1) * 4 + 2] * psigg[traj_j + 1]),
                 reduced_cor[3, ], btw_rho_gamma[[traj_i]][3] * sqrt(psi0[(traj_i - 1) * 4 + 3, (traj_i - 1) * 4 + 3] * psigg[traj_j + 1]),
                 btw_rho_gamma[[traj_i]][1] * sqrt(psi0[traj_j * 4 + 1, traj_j * 4 + 1] * psigg[traj_i]),
                 btw_rho_gamma[[traj_i]][2] * sqrt(psi0[traj_j * 4 + 2, traj_j * 4 + 2] * psigg[traj_i]),
                 btw_rho_gamma[[traj_i]][3] * sqrt(psi0[traj_j * 4 + 3, traj_j * 4 + 3] * psigg[traj_i]),
                 btw_gamma[traj_i] * sqrt(psigg[traj_i] * psigg[traj_j + 1])), byrow = T, nrow = 4)
      residuals[(traj_i + 1), traj_j] <- reduced_model$S$values[((traj_i - 1) * nT[traj_i] + 1), (traj_j * nT[traj_j] + 1)]
    }
  }
  psi0.s <- grad0 %*% t(psi0) %*% t(grad0)
  starts <- c(unlist(uni_mean0.s), psi0.s[row(psi0.s) >= col(psi0.s)], residuals[row(residuals) >= col(residuals)])
  ### Define manifest variables
  manifests <- c(paste0(traj_var[1], 1:nT[1]), paste0(traj_var[2], 1:nT[2]))
  ### Define latent variables
  latents <- paste0(rep(c("eta0s", "eta1s", "eta2s", "delta"), 2), rep(c("Y", "Z"), each = 4))
  outDefY <- outDefZ <- list(); outLoadsY1 <- outLoadsY2 <- outLoadsY3 <- outLoadsZ1 <- outLoadsZ2 <- outLoadsZ3 <- list()
  for(j in 1:nT[1]){
    outDefY[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[1], j), name = paste0("yt", j))
    outLoadsY1[[j]] <- mxAlgebraFromString(paste0("yt", j, " - mugY"), name = paste0("L1", j, "Y"))
    outLoadsY2[[j]] <- mxAlgebraFromString(paste0("abs(yt", j, " - mugY)"), name = paste0("L2", j, "Y"))
    outLoadsY3[[j]] <- mxAlgebraFromString(paste0("-mueta2sY * (yt", j, " - mugY)/abs(yt", j,
                                                  " - mugY) - mueta2sY"), name = paste0("L3", j, "Y"))
  }
  for(j in 1:nT[2]){
    outDefZ[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var[2], j), name = paste0("zt", j))
    outLoadsZ1[[j]] <- mxAlgebraFromString(paste0("zt", j, " - mugZ"), name = paste0("L1", j, "Z"))
    outLoadsZ2[[j]] <- mxAlgebraFromString(paste0("abs(zt", j, " - mugZ)"), name = paste0("L2", j, "Z"))
    outLoadsZ3[[j]] <- mxAlgebraFromString(paste0("-mueta2sZ * (zt", j, " - mugZ)/abs(zt", j,
                                                  " - mugZ) - mueta2sZ"), name = paste0("L3", j, "Z"))
  }
  ### Create a mxModel object
  model_mx <- mxModel("Bivariate Nonlinear Growth Model with Random Knots", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loadings from latent variables to manifests
                      mxPath(from = "eta0sY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1sY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:nT[1], "Y", "[1,1]")),
                      mxPath(from = "eta2sY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 1:nT[1], "Y", "[1,1]")),
                      mxPath(from = "deltaY", to = manifests[1:nT[1]], arrows = 1, free = F, values = 0,
                             labels = paste0("L3", 1:nT[1], "Y", "[1,1]")),
                      mxPath(from = "eta0sZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1sZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:nT[2], "Z", "[1,1]")),
                      mxPath(from = "eta2sZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 1:nT[2], "Z", "[1,1]")),
                      mxPath(from = "deltaZ", to = manifests[(nT[1] + 1):(nT[1] + nT[2])], arrows = 1, free = F, values = 0,
                             labels = paste0("L3", 1:nT[2], "Z", "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests[1:nT[1]], to = manifests[1:nT[1]],
                             arrows = 2, free = T, values = starts[45],
                             labels = "residualsY"),
                      mxPath(from = manifests[(nT[1] + 1):(nT[1] + nT[2])], to = manifests[(nT[1] + 1):(nT[1] + nT[2])],
                             arrows = 2, free = T, values = starts[47],
                             labels = "residualsZ"),
                      #### Define the covariances of residuals
                      mxPath(from = manifests[1:nT[1]], to = manifests[(nT[1] + 1):(nT[1] + nT[2])],
                             arrows = 2, free = T, values = starts[46],
                             labels = "residualsYZ"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents[c(1:3, 5:7)], arrows = 1, free = T, values = starts[c(1:3, 5:7)],
                             labels = c("mueta0sY", "mueta1sY", "mueta2sY",
                                        "mueta0sZ", "mueta1sZ", "mueta2sZ")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2, connect = "unique.pairs",
                             free = T, values = starts[c(9:44)],
                             labels = c("psi0sY0sY", "psi0sY1sY", "psi0sY2sY", "psi0sYgY", "psi0sY0sZ", "psi0sY1sZ", "psi0sY2sZ", "psi0sYgZ",
                                        "psi1sY1sY", "psi1sY2sY", "psi1sYgY", "psi1sY0sZ", "psi1sY1sZ", "psi1sY2sZ", "psi1sYgZ",
                                        "psi2sY2sY", "psi2sYgY", "psi2sY0sZ", "psi2sY1sZ", "psi2sY2sZ", "psi2sYgZ",
                                        "psigYgY", "psigY0sZ", "psigY1sZ", "psigY2sZ", "psigYgZ",
                                        "psi0sZ0sZ", "psi0sZ1sZ", "psi0sZ2sZ", "psi0sZgZ",
                                        "psi1sZ1sZ", "psi1sZ2sZ", "psi1sZgZ",
                                        "psi2sZ2sZ", "psi2sZgZ",
                                        "psigZgZ")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[4],
                               labels = "muknot_Y", name = "mugY"),
                      mxMatrix("Full", 1, 1, free = T, values = starts[8],
                               labels = "muknot_Z", name = "mugZ"),
                      outDefY, outDefZ, outLoadsY1, outLoadsY2, outLoadsY3, outLoadsZ1, outLoadsZ2, outLoadsZ3,
                      mxAlgebra(rbind(mueta0sY, mueta1sY, mueta2sY, mueta0sZ, mueta1sZ, mueta2sZ), name = "mean_s"),

                      mxAlgebra(rbind(cbind(psi0sY0sY, psi0sY1sY, psi0sY2sY, psi0sYgY, psi0sY0sZ, psi0sY1sZ, psi0sY2sZ, psi0sYgZ),
                                      cbind(psi0sY1sY, psi1sY1sY, psi1sY2sY, psi1sYgY, psi1sY0sZ, psi1sY1sZ, psi1sY2sZ, psi1sYgZ),
                                      cbind(psi0sY2sY, psi1sY2sY, psi2sY2sY, psi2sYgY, psi2sY0sZ, psi2sY1sZ, psi2sY2sZ, psi2sYgZ),
                                      cbind(psi0sYgY, psi1sYgY, psi2sYgY, psigYgY, psigY0sZ, psigY1sZ, psigY1sZ, psigYgZ),
                                      cbind(psi0sY0sZ, psi1sY0sZ, psi2sY0sZ, psigY0sZ, psi0sZ0sZ, psi0sZ1sZ, psi0sZ2sZ, psi0sZgZ),
                                      cbind(psi0sY1sZ, psi1sY1sZ, psi2sY1sZ, psigY1sZ, psi0sZ1sZ, psi1sZ1sZ, psi1sZ2sZ, psi1sZgZ),
                                      cbind(psi0sY2sZ, psi1sY2sZ, psi2sY2sZ, psigY2sZ, psi0sZ2sZ, psi1sZ2sZ, psi2sZ2sZ, psi2sZgZ),
                                      cbind(psi0sYgZ, psi1sYgZ, psi2sYgZ, psigYgZ, psi0sZgZ, psi1sZgZ, psi2sZgZ, psigZgZ)), name = "psi_s"),

                      mxAlgebra(rbind(cbind(1, -mugY, mugY, 0, 0, 0),
                                      cbind(0, 1, -1, 0, 0, 0),
                                      cbind(0, 1, 1, 0, 0, 0),
                                      cbind(0, 0, 0, 1, -mugZ, mugZ),
                                      cbind(0, 0, 0, 0, 1, -1),
                                      cbind(0, 0, 0, 0, 1, 1)), name = "func"),

                      mxAlgebra(rbind(cbind(1, -mugY, mugY, 0, 0, 0, 0, 0),
                                      cbind(0, 1, -1, 0, 0, 0, 0, 0),
                                      cbind(0, 1, 1, 0, 0, 0, 0, 0),
                                      cbind(0, 0, 0, 1, 0, 0, 0, 0),
                                      cbind(0, 0, 0, 0, 1, -mugZ, mugZ, 0),
                                      cbind(0, 0, 0, 0, 0, 1, -1, 0),
                                      cbind(0, 0, 0, 0, 0, 1, 1, 0),
                                      cbind(0, 0, 0, 0, 0, 0, 0, 1)), name = "grad"),
                      mxAlgebra(func %*% mean_s, name = "mean"),
                      mxAlgebra(grad %*% psi_s %*% t(grad), name = "psi"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
    model.para <- summary(model)$parameters[, c(1, 5, 6)]
    mean_est <- mean_se <- list()
    mug <- model.para$Estimate[c(grep("muknot", model.para$name))]
    psigg <- model.para$Std.Error[c(grep("muknot", model.para$name))]
    mean.se <- mxSE(mean, model)
    for (traj in 1:length(traj_var)){
      mean_est[[length(mean_est) + 1]] <- c(model$mean$result[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)], mug[traj])
      mean_se[[length(mean_se) + 1]] <- c(mean.se[((traj - 1) * 3 + 1):((traj - 1) * 3 + 3)], psigg[traj])
    }
    model.est <- c(unlist(mean_est), model$psi$result[row(model$psi$result) >= col(model$psi$result)],
                   model.para[grep("residuals", model.para$name), 2])
    psi.se <- mxSE(psi, model)
    model.se <- c(unlist(mean_se), psi.se[row(psi.se) >= col(psi.se)],
                  model.para[grep("residuals", model.para$name), 3])
    estimate_out <- data.frame(Name = paraNames, Estimate = model.est, SE = model.se)
    return(list(model, estimate_out))
  }
  else{
    return(model)
  }
}


