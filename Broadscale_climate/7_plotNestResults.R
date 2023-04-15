library(ggplot2)

pproblabels <- c('-', '>0.975', '>0.995','>0.9995')

load(file="data/nestdata.RData")
load(file="results/nestResults.Rdata")

res.table <- readRDS(file="results/restable.Rdata")  

png(filename = "figures/all_effects.png",width = 1200, height = 800)
p <- res.table %>%
  filter(data=="all") %>%
  ggplot(aes(x=effect,y=model,color = postprob, shape = postprob)) + geom_point(size=3) + 
  geom_errorbar(aes(xmin = lb,xmax=ub)) + geom_vline(xintercept=0,linetype="dashed") + 
  scale_color_discrete(labels=pproblabels) + 
  scale_shape_discrete(labels=pproblabels)
print(p)  
dev.off()


png(filename = "figures/noAl_effects.png",width = 1200, height = 800)
p <- res.table %>%
  filter(data=="noAl") %>%
  ggplot(aes(x=effect,y=model,color = postprob, shape = postprob)) + geom_point(size=3) + 
  geom_errorbar(aes(xmin = lb,xmax=ub)) + geom_vline(xintercept=0,linetype="dashed") + 
  scale_color_discrete(labels=pproblabels) + 
  scale_shape_discrete(labels=pproblabels)
print(p)  
dev.off()


png(filename = "figures/effects.png",width = 1200, height = 800)
p <- res.table %>%
  ggplot(aes(x=effect,y=model,color = postprob, shape = postprob)) + geom_point(size=3) + 
  geom_errorbar(aes(xmin = lb,xmax=ub)) + geom_vline(xintercept=0,linetype="dashed") + 
  scale_color_discrete(labels=pproblabels) + 
  scale_shape_discrete(labels=pproblabels) +
  facet_grid(cols = vars(data))
print(p)  
dev.off()

png(filename = "figures/effects_20yr.png",width = 1200, height = 800)
p <- res.table %>%
  filter(year %in% c("20yr","Lat","Long")) %>%
  ggplot(aes(x=effect,y=model,color = postprob, shape = postprob)) + geom_point(size=3) + 
  geom_errorbar(aes(xmin = lb,xmax=ub)) + geom_vline(xintercept=0,linetype="dashed") + 
  scale_color_discrete(labels=pproblabels) + 
  scale_shape_discrete(labels=pproblabels) +
  facet_grid(cols = vars(data))
print(p)  
dev.off()

png(filename = "figures/all_effects_20yr.png",width = 1200, height = 800)
p <- res.table %>%
  filter(data=="all" & year %in% c("20yr","Lat","Long")) %>%
  ggplot(aes(x=effect,y=model,color = postprob, shape = postprob)) + geom_point(size=3) + 
  geom_errorbar(aes(xmin = lb,xmax=ub)) + geom_vline(xintercept=0,linetype="dashed") + 
  scale_color_discrete(labels=pproblabels) + 
  scale_shape_discrete(labels=pproblabels) 
print(p)  
dev.off()


png(filename = "figures/elpd_effect.png",width = 1200, height = 800)
p <- res.table %>%
  ggplot(aes(x=elpd_diff,y=abs(effect),shape=data,color=postprob)) +
  geom_point(size=5) +
  scale_color_discrete(labels=pproblabels) + 
  geom_line(aes(group=model),show.legend = F, color="black", linetype="dotted")
print(p)
dev.off()
