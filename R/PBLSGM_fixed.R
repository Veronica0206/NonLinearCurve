getPBLSGM_Fixed <- function(dat, nT, traj_var, t_var, res_ratio = rep(4, 2), rho = rep(0.1, 2), btw_rho = rep(0.1, 1),
                            btw_res = rep(0.3, 1), starts = NA, loc = 1, scale = 0.25, extraTries = NA,
                            uni_paraNames, paraNames, optimizer = "CSOLNP"){
  mxOption(model = NULL, key = "Default optimizer", optimizer, reset = FALSE)
  if (any(is.na(starts))){
    uni_residual <- list();
    uni_mean0.s <- uni_psi0 <-uni_grad0 <- list()
    for (traj in 1:length(traj_var)){
      BLSGM_F <- getBLSGM_Fixed(dat = dat, nT = nT[traj], traj_var = traj_var[traj], t_var = t_var[traj], paraNames = uni_paraNames)
      ### Mean vector (in the reparameterized framework)
      uni_mean0.s[[length(uni_mean0.s) + 1]] <- c(BLSGM_F[[1]]$mean_s$result, BLSGM_F[[1]]$mug$values)
      ### var-cov matrix (in the original framework)
      ###################
      uni_psi0[[length(uni_psi0) + 1]] <- BLSGM_F[[1]]$psi$result
      uni_grad0[[length(uni_grad0) + 1]] <- matrix(c(1, BLSGM_F[[1]]$mug$values, 0,
                                                     0, 0.5, 0.5,
                                                     0, -0.5, 0.5), nrow = 3, byrow = T)
      uni_residual[[length(uni_residual) + 1]] <- BLSGM_F[[1]]$output$estimate[1]

    }
    psi0 <- as.matrix(Matrix::bdiag(uni_psi0))
    grad0 <- as.matrix(Matrix::bdiag(uni_grad0))
    residuals <- as.matrix(Matrix::bdiag(uni_residual))
    for (traj_i in 1:(length(traj_var) - 1)){
      for (traj_j in traj_i:(length(traj_var) - 1)){
        psi0[((traj_i - 1) * 3 + 1):((traj_i - 1) * 3 + 3), (traj_j * 3 + 1):(traj_j * 3 + 3)] <-
          btw_rho[traj_i] * sqrt(diag(uni_psi0[[traj_i]]) %*% t(diag(uni_psi0[[traj_j + 1]])))
        residuals[traj_i, (traj_j + 1)] <- btw_res * sqrt(uni_residual[[traj_i]] * uni_residual[[traj_j + 1]])
      }
    }
    psi0.s <- grad0 %*% t(psi0) %*% t(grad0)
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
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
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


