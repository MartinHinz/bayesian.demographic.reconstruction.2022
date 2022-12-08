sweep_run <- function(model_constants_list, this_model_data, model_code) {

  batches <- 100000 # run in batches of this length
  thinning <- 10 # save only each ... of the chains

  nrow_data <- nrow(this_model_data)

  sweep_results <- list()

  for (i in seq_along(model_constants_list)) {

    this_model_constants <- model_constants_list[[i]]

    this_model_data$constraint_lambda_lower <- 1
    this_model_data$constraint_lambda_upper <- 1

    r <- log(1 - (10^(1/this_model_constants$nYears) - 1))

    model_inits <- list(
      lambda = rep(r,nrow_data),
      PopDens = rep(this_model_constants$nEnd,nrow_data),
      nSites = rep(50,nrow_data)
    )

    model_monitors <- c("PopDens", "p")

    # Spare 4 cores from the calculation
    ncore <- detectCores() - 4
    ncore <- max(ncore, 1)

    # Register the cluster
    cl <- makeCluster(ncore)
    registerDoParallel(cl)

    # Export common values for the cluster
    clusterExport(cl,
                  c("model_code",
                    "model_inits",
                    "this_model_data",
                    "this_model_constants",
                    "model_monitors",
                    "batches",
                    "thinning"
                  ),
                  envir=environment()
    )

    # Set random seeds
    for (j in seq_along(cl)) {
      set.seed(j)
    }

    start_time <- Sys.time()

    mcmcSamples <- clusterEvalQ(cl, {

      # Load necessary libraries
      library(nimble)
      library(coda)

      # initiate the model with the parameters and dates
      model <- nimbleModel(code=model_code,
                           data=this_model_data,
                           constants = this_model_constants,
                           inits = model_inits)

      # Compile the model
      Cmodel <- compileNimble(model)

      # Configure the model
      modelConf <- configureMCMC(model, thin = thinning)

      # Add the monitor(s)
      modelConf$addMonitors(model_monitors)

      # Build the mcmc
      modelMCMC <- buildMCMC(modelConf)

      # Compile the final model
      CmodelMCMC <- compileNimble(modelMCMC, project = model)

      # Run the model for the number of iterations specified in batches
      CmodelMCMC$run(batches, reset = FALSE)
      return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
    })


    end_time <- Sys.time()
    end_time - start_time

    # Initialize with non-convergence
    converged <- F

    # minimum psrf until the model is considered converged
    min_convergence <- 1.1

    start_time <- Sys.time()

    # As long as not converged
    while (!converged) {

      # run the model for another batch, with resetting the values (burn-in)
      mcmcSamples <- clusterEvalQ(cl, {
        CmodelMCMC$run(batches, reset = FALSE, resetMV = TRUE)
        return(as.mcmc(as.matrix(CmodelMCMC$mvSamples)))
      })

      # convert the resulting list to an mcmc.list object
      mcmcSamples <- as.mcmc.list(mcmcSamples)

      # Check the psrf values from the Gelman and Rubin's convergence diagnostic
      gelman <- gelman.diag(mcmcSamples, multivariate = F)
      psrf <- gelman$psrf[,1]

      # The model is considered converged when
      # all psrf are lower than the minimum convergence criterion
      converged <- all(psrf<min_convergence, na.rm=T)

      # optional: create an mcmc traceplot for visual inspection
      #MCMCtrace(window(mcmcSamples, start = (batches/thinning)-min((batches/thinning),5000)), params = "PopDens", pdf = F)
    }

    end_time <- Sys.time()
    end_time - start_time

    result <- MCMCsummary(mcmcSamples, HPD = T, params = "PopDens")
    sweep_results[[i]] <- result
    stopCluster(cl)
  }

  return(sweep_results)
}
