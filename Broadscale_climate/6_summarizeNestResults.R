library(ggplot2)
library(brms)
library(dplyr)

load(file="data/nestdata.RData")
load(file="results/nestResults.Rdata")

nX <- length(Xvars.all)
res.table <- data.frame()


datasets <- c("all","noAl")

for (da in 1:2) {

  result.table  <- as.data.frame((loo_compare(loos[,da])))
  result.table <- result.table[c("elpd_diff","se_diff")]
  result.table$model <- as.factor(row.names(result.table))
  result.table$pprob = NULL
  result.table$effect = NULL
  result.table$data <- datasets[da]
  result.table$bvar <- NULL
  result.table$year <- NULL
  
  for (n in 1:nX) {
    X <- Xvars.all[n]
    fitx <- fits[[n,da]]
    a <- hypothesis(fitx,paste0(X," > 0"))
    rowx <- rownames(result.table)==X
    pprob <- a$hypothesis$Post.Prob
    result.table$pprob[rowx] = (max(pprob,1-pprob))
    result.table$effect[rowx] = (summary(fitx)$fixed$Estimate[2])
    result.table$lb[rowx] = (summary(fitx)$fixed$`l-95% CI`[2])
    result.table$ub[rowx] = (summary(fitx)$fixed$`u-95% CI`[2])
    
    
    xvrow <- Xvars.all.table$var==X
    result.table$bvar[rowx] = Xvars.all.table$var.base[xvrow]
    result.table$year[rowx] = Xvars.all.table$year[xvrow]
  }
  res.table <- rbind(res.table,result.table)
}



res.table <- res.table %>%
  mutate(postprob = factor((pprob>0.975)+(pprob>0.995)+(pprob>0.9995)))  %>%
  mutate(model=factor(model, levels=Xvars.all))

saveRDS(res.table,file="results/restable.Rdata")  
