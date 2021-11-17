getBLSGM_Random <- function(dat, nT, traj_var, t_var, res_ratio = 4, rho = 0.1, rho_gamma = 0, starts = NA,
                            loc = 1, scale = 0.25, extraTries = NA, original = T, paraNames = NA, optimizer = "CSOLNP"){
  mxOption(model = NULL, key = "Default optimizer", optimizer, reset = FALSE)
  if (I(original & any(is.na(paraNames)))){
    print("Please enter the original parameters if want to obtain them!")
    break
  }
  if (any(is.na(starts))){
    dat_traj <- dat[, paste0(traj_var, 1:nT)]
    dat_time <- dat[, paste0(t_var, 1:nT)]
    ### Obtain estimates from the reduced model
    reduced_model <- getBLSGM_Fixed(dat = dat, nT = nT, traj_var = traj_var, t_var = t_var, res_ratio = res_ratio, rho = rho,
                                    loc = loc, scale = scale, extraTries = extraTries, original = F)
    ### Obtain knot variance
    mean_time <- apply(dat_time, 2, mean)
    mug <- reduced_model$mug$values
    for (t in 1:length(mean_time)){
      if (mean_time[t] <= mug & mean_time[t + 1] >= mug){
        psigg <- var(apply(dat_time[, c(t, t + 1)], 1, mean))
        stop
      }
    }
    ### Initial values in the original parameter space
    mean0 <- c(reduced_model$mean$result, mug)
    reduced_psi0 <- reduced_model$psi$result
    psi0 <- matrix(c(reduced_psi0[1, ], rho_gamma * sqrt(reduced_psi0[1, 1] * psigg),
                     reduced_psi0[2, ], rho_gamma * sqrt(reduced_psi0[2, 2] * psigg),
                     reduced_psi0[3, ], rho_gamma * sqrt(reduced_psi0[3, 3] * psigg),
                     rho_gamma * sqrt(reduced_psi0[1, 1] * psigg), rho_gamma * sqrt(reduced_psi0[2, 2] * psigg),
                     rho_gamma * sqrt(reduced_psi0[3, 3] * psigg), psigg), nrow = 4)
    ### Transformed matrices obtained by multivariate Delta method
    #### For mean vector
    func0 <- matrix(c(1, mean0[4], 0, 0,
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = T)
    #### For var-cov matrix
    grad0 <- matrix(c(1, mean0[4], 0, mean0[2],
                      0, 0.5, 0.5, 0,
                      0, -0.5, 0.5, 0,
                      0, 0, 0, 1), nrow = 4, byrow = T)

    mean0.s <- func0[1:3, 1:3] %*% mean0[1:3]
    psi0.s <- grad0 %*% psi0 %*% t(grad0)
    true.s <- c(mean0.s[1:3], mean0[4], psi0.s[row(psi0.s) >= col(psi0.s)])
    starts <- c(true.s, reduced_model$S$values[1, 1])
  }
  ### Define manifest variables
  manifests <- paste0(traj_var, 1:nT)
  ### Define latent variables
  latents <- c("eta0s", "eta1s", "eta2s", "delta")
  outDef <- list(); outLoads1 <- list(); outLoads2 <- list(); outLoads3 <- list()
  for(j in 1:nT){
    outDef[[j]] <- mxMatrix("Full", 1, 1, free = F, labels = paste0("data.", t_var, j), name = paste0("t", j))
    outLoads1[[j]] <- mxAlgebraFromString(paste0("t", j, " - mug"), name = paste0("L1", j))
    outLoads2[[j]] <- mxAlgebraFromString(paste0("abs(t", j, " - mug)"), name = paste0("L2", j))
    outLoads3[[j]] <- mxAlgebraFromString(paste0("-mueta2s * (t", j, " - mug)/abs(t", j,
                                                 " - mug) - mueta1s"), name = paste0("L3", j))
  }
  ### Create a mxModel object
  model_mx <- mxModel("Estimate a random knot", type = "RAM",
                      manifestVars = manifests, latentVars = latents,
                      mxData(observed = dat, type = "raw"),
                      #### Define factor loadings from latent variables to manifests
                      mxPath(from = "eta0s", to = manifests, arrows = 1, free = F, values = 1),
                      mxPath(from = "eta1s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L1", 1:nT, "[1,1]")),
                      mxPath(from = "eta2s", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L2", 1:nT, "[1,1]")),
                      mxPath(from = "delta", to = manifests, arrows = 1, free = F, values = 0,
                             labels = paste0("L3", 1:nT, "[1,1]")),
                      #### Define the variances of residuals
                      mxPath(from = manifests, to = manifests, arrows = 2, free = T, values = starts[15],
                             labels = "residuals"),
                      #### Define means of latent variables
                      mxPath(from = "one", to = latents[1:3], arrows = 1, free = T, values = starts[1:3],
                             labels = c("mueta0s", "mueta1s", "mueta2s")),
                      #### Define var-cov matrix of latent variables
                      mxPath(from = latents, to = latents, arrows = 2,
                             connect = "unique.pairs", free = T,
                             values = starts[5:14],
                             labels = c("psi0s0s", "psi0s1s", "psi0s2s", "psi0sg", "psi1s1s",
                                        "psi1s2s", "psi1sg", "psi2s2s", "psi2sg", "psigg")),
                      #### Add additional parameter and constraints
                      mxMatrix("Full", 1, 1, free = T, values = starts[4],
                               labels = "muknot", name = "mug"),
                      outDef, outLoads1, outLoads2, outLoads3,
                      mxAlgebra(rbind(mueta0s, mueta1s, mueta2s), name = "mean_s"),
                      mxAlgebra(rbind(cbind(psi0s0s, psi0s1s, psi0s2s, psi0sg),
                                      cbind(psi0s1s, psi1s1s, psi1s2s, psi1sg),
                                      cbind(psi0s2s, psi1s2s, psi2s2s, psi2sg),
                                      cbind(psi0sg, psi1sg, psi2sg, psigg)), name = "psi_s"),
                      mxAlgebra(rbind(cbind(1, -mug, mug, 0),
                                      cbind(0, 1, -1, 0),
                                      cbind(0, 1, 1, 0),
                                      cbind(0, 0, 0, 1)), name = "func"),
                      mxAlgebra(rbind(cbind(1, -mug, mug, mueta2s - mueta1s),
                                      cbind(0, 1, -1, 0),
                                      cbind(0, 1, 1, 0),
                                      cbind(0, 0, 0, 1)), name = "grad"),
                      mxAlgebra(func[1:3, 1:3] %*% mean_s, name = "mean"),
                      mxAlgebra(grad %*% psi_s %*% t(grad), name = "psi"))
  if (!is.na(extraTries)){
    model <- mxTryHard(model_mx, extraTries = extraTries, OKstatuscodes = 0, loc = loc, scale = scale)
  }
  else{
    model <- mxRun(model_mx)
  }
  if(original){
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
  else{
    return(model)
  }
}