antdata <- read.csv2("data/Winter_Spring_Temp_Prec_33antnests_for_modelling_roundedCoord.csv", header = TRUE, sep = ";")

antdata$hybrid <- 1*(antdata$Genomic_assignment=="Admixed1")

Xvars.all <- c()

Xvars.all <- c(Xvars.all,"ETRS.TM35FIN_Y","ETRS.TM35FIN_X")

tempyears <- c("1yr","5yr","20yr")
Climvars <- c("spring_meanT","spring_Nfrostdays","spring_perc25T",  "spring_meanP",
              "winter_meanP", "winter_meanT", "winter_Ndays0plus","winter_perc75T")

for (cvar in Climvars) {
  for (yr in tempyears) {
    Xvars.all <- c(Xvars.all,paste0(cvar,"_",yr))
  }
}

Xvars.all.table <- data.frame(var = Xvars.all,
                              var.base = c("Loc","Loc",rep(Climvars,each=3)),
                              year = c("ETRS.TM35FIN_Y","ETRS.TM35FIN_X",rep(tempyears,8)))

# Data without Åland nests
mainlandNests <- which(antdata$ETRS.TM35FIN_X>150000)
subdata_noAL <- antdata[mainlandNests,]
subdata_noAL[Xvars.all] <- scale(subdata_noAL[Xvars.all])

save(antdata, subdata_noAL, Xvars.all,Xvars.all.table,file="data/nestdata.RData")
