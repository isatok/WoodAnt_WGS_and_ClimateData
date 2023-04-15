setwd()
setwd("/Users/marijone/Documents/Mirkka Documents BAU/BAU projects/Kulmuni Ants/R scripts for GitHub and input data broadscale/")
### Calculation of Winter and Spring seasonal temperature and precipitation per ant nest ###

# First, import annual winter temperature data and calculate variables from this for ant modelling.
# Calculated variables are Dec-Feb daily mean temperatures and N days >0C calculated for the ant observation year itself 
# and averaged over temporal windows extending to 5, 10, 15 and 20 years before sampling.

winter_data = read.csv2("output data/WinterTemperatures_antnests_1986-2019.csv")

varnames = c("winter_meanT_", "winter_Ndays0plus_", "winter_perc75T_")
varnames_short = c("_mean", "_Ndays_0plus", "_perc75")
windows = c(1, 5, 10, 15, 20)

winterTemps = data.frame(matrix(NA, nrow = nrow(winter_data), ncol = 15))
names(winterTemps) = paste0(sort(rep(varnames,5)), rep(paste0(windows,"yr"), 3))

for(i in 1:length(windows))
{
  window = windows[i]
  
  for(n in 1:nrow(winter_data))
  {
    endyear = winter_data$Sampling_year[n]
    
    for(j in 1:3)
    {
      var = varnames_short[j]
      varname = varnames[j]
      
      subset = winter_data[n, grep(var, names(winter_data))]
      subset = subset[, (grep(endyear, names(subset))-(window-1)):grep(endyear, names(subset))]
      data_column = paste0(varname, window, "yr")
      winterTemps[n, grep(data_column, names(winterTemps))] = mean(as.numeric(subset))
    }
  }
}

winterTemps = data.frame(winter_data[, 1:7], winterTemps)


# Then, import annual spring temperature data and extract variables from this for ant modelling.
# Calculated variables are Apr-May daily mean temp and N days <0C calculated for the ant observation year itself 
# and averaged over temporal windows extending to 5, 10, 15 and 20 years before sampling.

spring_data = read.csv2("output data/SpringTemperatures_antnests_1986-2019.csv")

varnames = c("spring_meanT_", "spring_Nfrostdays_", "spring_perc25T_")
varnames_short = c("_mean", "_Nfrostdays", "_perc25")
windows = c(1, 5, 10, 15, 20)

springTemps = data.frame(matrix(NA, nrow = nrow(spring_data), ncol = 15))
names(springTemps) = paste0(sort(rep(varnames,5)), rep(paste0(windows,"yr"), 3))

for(i in 1:length(windows))
{
  window = windows[i]
  
  for(n in 1:nrow(spring_data))
  {
    endyear = spring_data$Sampling_year[n]
    
    for(j in 1:3)
    {
      var = varnames_short[j]
      varname = varnames[j]
      
      subset = spring_data[n, grep(var, names(spring_data))]
      subset = subset[, (grep(endyear, names(subset))-(window-1)):grep(endyear, names(subset))]
      data_column = paste0(varname, window, "yr")
      springTemps[n, grep(data_column, names(springTemps))] = mean(as.numeric(subset))
    }
  }
}


### Calculation of Winter and Spring mean precipitation (mean of daily sums) per ant nest ###

# Next, import annual winter precipitation data and means of Dec-Feb daily precipitation sums for the ant observation year itself 
# and averaged over temporal windows extending to 5, 10, 15 and 20 years before sampling.

winter_data = read.csv2("output data/WinterPrecipitation_antnests_1986-2019.csv")

varname = c("winter_meanP_")
varname_short = c("_mean")
windows = c(1, 5, 10, 15, 20)

winterPrec = data.frame(matrix(NA, nrow = nrow(winter_data), ncol = 5))
names(winterPrec) = paste0(rep(varname,5), paste0(windows,"yr"))

for(i in 1:length(windows))
{
  window = windows[i]
  
  for(n in 1:nrow(winter_data))
  {
    endyear = winter_data$Sampling_year[n]
      
      subset = winter_data[n, grep(varname_short, names(winter_data))]
      subset = subset[, (grep(endyear, names(subset))-(window-1)):grep(endyear, names(subset))]
      data_column = paste0(varname, window, "yr")
      winterPrec[n, grep(data_column, names(winterPrec))] = mean(as.numeric(subset))
  }
}


# Finally, import annual spring precipitation data and means of April-May daily precipitation sums for the ant observation year itself 
# and averaged over temporal windows extending to 5, 10, 15 and 20 years before sampling.

spring_data = read.csv2("output data/SpringPrecipitation_antnests_1986-2019.csv")

varname = c("spring_meanP_")
varname_short = c("_mean")
windows = c(1, 5, 10, 15, 20)

springPrec = data.frame(matrix(NA, nrow = nrow(spring_data), ncol = 5))
names(springPrec) = paste0(rep(varname,5), paste0(windows,"yr"))

for(i in 1:length(windows))
{
  window = windows[i]
  
  for(n in 1:nrow(spring_data))
  {
    endyear = spring_data$Sampling_year[n]
    
    subset = spring_data[n, grep(varname_short, names(spring_data))]
    subset = subset[, (grep(endyear, names(subset))-(window-1)):grep(endyear, names(subset))]
    data_column = paste0(varname, window, "yr")
    springPrec[n, grep(data_column, names(springPrec))] = mean(as.numeric(subset))
  }
}

nests_springwinterT_P = data.frame(winterTemps, winterPrec, springTemps, springPrec)

write.csv2(nests_springwinterT_P, "output data/Winter_Spring_Temp_Prec_33antnests_for_modelling_roundedCoord.csv", row.names = F)
