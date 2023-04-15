library(brms)
library(ggplot2)

load(file="data/nestdata.RData")

nX <- length(Xvars.all)

fits <- as.list(1:(2*nX))
loos <- as.list(1:(2*nX))
dim(fits) <- c(nX,2)
dim(loos) <- c(nX,2)

row.names(fits) <- Xvars.all
row.names(loos) <- Xvars.all


prior <- c(set_prior("normal(0,2)", class = "Intercept"),
            set_prior("normal(0,2)", class = "b"))
niter <- 4000
seed <- 2371
cores <- 1
silent <- 2
refresh <- 0

sink("results/nestAnalyses_log.txt", type="output")

for (da in 1:2) {
  if (da==1) {
    data <- subdata
    
    cat("FULL DATA \n")
    cat("\n")
  } else if (da==2) {
    data <- subdata_noAL
    
    cat("WITHOUT ÅLAND \n")
    cat("\n")
  }

  for (n in 1:nX) {
    X <- Xvars.all[n]
    cat("\n")
    cat(X)
    cat("\n")
    cat("\n")
    model <- paste0("hybrid ~ ", X)
    fits[[n,da]] <- brm(model,data=data, family = bernoulli, prior=prior, seed = seed,
                        iter = niter, cores = cores, silent = silent,
                        refresh = refresh, file = paste0("models/Nest_",da,"_",X))
    
    print(summary(fits[[n,da]]))
    cat("\n")
    fits[[n,da]] <- add_criterion(fits[[n,da]],"loo")
    loos[[n,da]] <- fits[[n,da]]$criteria$loo
    print(loo(fits[[n,da]]))
    cat("\n")
  }
}

sink()
save(fits,loos,file="results/nestResults.Rdata")



