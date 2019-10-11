###### FILAMO Ecological Forecasting October 2019
##### Forecasting and quantifying species range shifts
#### Forecasting changes in overlap between two species
### Gemma Carroll 

install.packages("mgcv")
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("maps")
install.packages("dismo")

require(mgcv)
require(ggplot2)
require(gridExtra)
require(maps)
require(dismo)


### read in species observation data set
#### arrowtooth flounder Eastern Bering Sea 1982 - 2017
### NOAA Alaska Fisheries Science Center groundtrawl survey
data <- read.csv('flounder.csv')

#### turn abundance data into presence/absence to model occurrence
data$pres <- as.numeric(data$CPUE_NOHA > 0)
names(data) <- c("vessel", "haul", "year", "cruise", "lat", "lon", "station_id", "gear_depth", "bottom_depth", "temp_surface", "temp_gear", "cpue", "pres")

#### remove NA values 
data <- data[complete.cases(data),]

### in the file 'data', we now have locations (latitude and longitude) with presence(1) and absence(0) observations for arrowtooth flounder in the Eastern Bering Sea
world <- map_data("world")

ggplot() + 
  geom_point(data = data, aes(x = lon, y = lat, colour = as.factor(pres))) + 
  scale_colour_manual(values = c("#400b0b", "#ffdb00"), name = "pres_abs") + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year) + 
  theme_classic()

### and take a take a look at the spatial distribution of bottom temperature (temp_gear) during the time series
#### could also plot surface temperature and bathymetry if interested
tempCols <- colorRampPalette(c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))

ggplot() + 
  geom_point(data = data, aes(x = lon, y = lat, colour = temp_gear)) + 
  scale_colour_gradientn(colors = tempCols(20)) + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)  +
  theme_classic() 


### BUILDING A SPECIES DISTRIBUTION MODEL

##### quickly check out species ~ environment relationships
#### histogram of arrowtooth flounder presence as a function of:
#### temperature at gear
ggplot()+
  geom_histogram(data = data[data$pres == 1,], aes(temp_gear),bins = 20, fill = "#F21A00", col = "grey20") + 
  theme_classic()

#### surface temp
ggplot()+
  geom_histogram(data = data[data$pres == 1,], aes(temp_surface),bins = 20, fill = "#E1AF00", col = "grey20") + 
  theme_classic()

#### bottom depth
ggplot()+
  geom_histogram(data = data[data$pres == 1,], aes(bottom_depth),bins = 20, fill = "#78B7C5", col = "grey20") + 
  theme_classic()

### check for correlation between covariates
cor(data$temp_gear, data$temp_surface)
cor(data$temp_gear, data$bottom_depth)

### fit a binomial generalised additive model (GAM) to occurrence data using the package "mgcv"
### "s(x)" defines a smoothing term that allows the relationship between occurrence and the covariate to be flexible and non-linear
gam1 <- gam(pres ~ s(lon, lat) + s(year), family = "binomial", data = data)
gam2 <- gam(pres ~ s(lon, lat)  + s(year) + s(temp_gear), family = "binomial", data = data)
gam3 <- gam(pres ~ s(lon, lat)  + s(year) + s(temp_gear) + s(temp_surface), family = "binomial", data = data)
gam4 <- gam(pres ~ s(lon, lat) + s(year) + s(temp_gear) + s(temp_surface) + s(bottom_depth), family = "binomial", data = data)

#### Akaike's Information Criterion for model selection
#### lowest value is best model (most accurate relative to complexity)
AIC(gam1, gam2, gam3, gam4)

### take a look at the model summary - looks pretty good! 
### temperature explains around 80% of the variability in species occurrence
### (in the real world you would spend a lot more time on model fitting, model selection, validation, testing etc) 
summary(gam4)

#### look at partial plots for each covariate in this model, they describe species ~ environment relationships
plot(gam4)

### save best model formula 
formula <- pres ~ s(lon, lat) + s(year) + s(temp_gear) + s(temp_surface) + s(bottom_depth)

#### perform k-folds cross validation
#### randomly partition the data into 10 parts, train the data on 9 of those parts and test on the 10th
eval_kfold <- function(data, formula){
  data = data
  data$Kset <- kfold(data,10) #randomly allocate k groups
  kfold_gam <- as.data.frame(matrix(data=0,nrow=10,ncol=3))
  colnames(kfold_gam) <- c("k","AUC","TSS")
  counter=1
  for (k in 1:10){
  print(k)
  data_train <- data[data$Kset!=k,]
  data_test <- data[data$Kset==k,]
  data.kfolds <- mgcv::gam(formula, family = "binomial", data = data_train)
  preds <- as.data.frame(predict(data.kfolds, data_test, se=TRUE, type="response"))
  d <- cbind(data_test$pres, preds)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  kfold_gam[counter,1] <- k
  kfold_gam[counter,2] <- e@auc
  kfold_gam[counter,3] <- max(e@TPR + e@TNR-1)
  counter=counter+1 
  }
  return(kfold_gam)}

eval_kfold(data, formula)

#### looks like the model performs pretty well! 


### let's take the model fitted values to see what they look like geographically 
data$gam.predict <- predict(gam4, type = "response", newdata = data)

occurrenceCols <- colorRampPalette(c("#400b0b", "#a12424","#ee7b06", "#ffa904", "#ffdb00"))

### instead of presence and absence, we now have a surface of probability of occurrence values (0-1) for each location, derived from the relationship between occurrence and environment
ggplot() + 
  geom_point(data = data, aes(x = lon, y = lat, colour = gam.predict)) + 
  scale_colour_gradientn(colours = occurrenceCols(20), name = "P(occurrence)") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

### we can also plot the spatially-explicit uncertainty in our model results from the standard error around the fitted values estimated by the GAM
#### perhaps unsurprisingly tells us that the greatest uncertainty exists along the margin of the distribution
gam.se <- predict(gam4, type = "response", se.fit = "TRUE", newdata = data)
data$gam.se = gam.se$se.fit
  
ggplot() + 
  geom_point(data = data, aes(x = lon, y = lat, colour = gam.se)) + 
  scale_colour_gradientn(colours = occurrenceCols(20), name = "error") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

#### hold the most recent year in this data set (2017) out of training and see how the model does
#### "hindcast" approach on 2017 as an out-of-sample test data set
data_train <- data[data$year < 2017,]
data_test <- data[data$year == 2017,]

gam4_lo2017 <- gam(formula, family = "binomial", data = data_train)

##### "hindcasting" to determine model performance on years left out of the training process
data_test$gam.predict <- predict(gam4_lo2017, type = "response", newdata = data_test)

##### calculate the difference between observations and predictions 
data_test$diff <- data_test$pres - data_test$gam.predict

p1 <- ggplot() + 
  geom_point(data = data_test, aes(x = lon, y = lat, colour = gam.predict)) + 
  scale_colour_gradientn(colours = occurrenceCols(20), name = "P(occur)") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

p2 <- ggplot() + 
  geom_point(data = data_test, aes(x = lon, y = lat, colour = as.factor(pres))) + 
  scale_colour_manual(values = c("#400b0b", "#ffdb00"), name = "pres_abs") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

p3 <- ggplot() + 
  geom_point(data = data_test, aes(x = lon, y = lat, colour = diff)) + 
  scale_colour_gradient2(low = "red", mid = "grey90",
  high = "blue", midpoint = 0, name = "difference") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

grid.arrange(p2, p1, p3, ncol = 1)

#### hold another anomalous year out of training and see how the model does (2012 was a very cold year in this data set and the range of arrowtooth flounder was more restricted than usual)
data_train <- data[!data$year == 2012,]
data_test <- data[data$year == 2012,]

gam4_lo2012 <- gam(formula, family = "binomial", data = data_train)

##### "hindcasting" to determine model performance on 2012
data_test$gam.predict <- predict(gam4_lo2012, type = "response", newdata = data_test)
data_test$diff <- data_test$pres - data_test$gam.predict

p1 <- ggplot() + 
  geom_point(data = data_test, aes(x = lon, y = lat, colour = gam.predict)) + 
  scale_colour_gradientn(colours = occurrenceCols(20), name = "P(occur)") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

p2 <- ggplot() + 
  geom_point(data = data_test, aes(x = lon, y = lat, colour = as.factor(pres))) + 
  scale_colour_manual(values = c("#400b0b", "#ffdb00"), name = "pres_abs") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

p3 <- ggplot() + 
  geom_point(data = data_test, aes(x = lon, y = lat, colour = diff)) + 
  scale_colour_gradient2(low = "red", mid = "grey95",
                         high = "blue", midpoint = 0, name = "difference") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) + 
  facet_wrap(~year)

grid.arrange(p2, p1, p3, ncol = 1)


#### model seems to do a pretty good job of recapturing spatial patterns in out-of-sample test cases, giving us some confidence to apply the model to new years

#### iteratively test how well the model would forecast the following year starting from 1982
### refit the model each time on increasingly more data
#### turn year into a linear parameter to reduce degrees of freedom so we can start early in the time series
formula <- pres ~ s(lon, lat) + s(temp_gear) + s(temp_surface) + s(bottom_depth) + year

forecast_iter <- function(data, formula){
  data = data
  perform <- as.data.frame(matrix(data=0,nrow=length(seq(1983:2016)),ncol=3))
  colnames(perform) <- c("year","AUC","TSS")
  counter = 0 
for (y in 1983:2017){
  y1 = 1982
  yi = seq(y1, y1 + counter,1)
  print (yi)
  data_train <- data[data$year %in% yi,]
  data_test <- data[data$year == max(yi) + 1,]
  gam <- gam(formula, family = "binomial", data = data_train)
  preds <- predict(gam, type = "response", newdata = data_test)
  d <- cbind(data_test$pres, preds)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  perform[counter+1,1] <- max(yi) + 1
  perform[counter+1,2] <- e@auc
  perform[counter+1,3] <- max(e@TPR + e@TNR-1)
  counter = counter + 1
}
  return(perform)}

perform <- forecast_iter(data, formula)

#### plot the AUC as a function of the amount of information that goes into the forecast
#### reduce process error?
ggplot(data = perform, aes(year, AUC)) + 
  geom_point() + 
  theme_classic()

ggplot(data = perform, aes(year, TSS)) + 
  geom_point() + 
  theme_classic()

##### so at this point we can either forecast species distributions onto existing environmental data where we don't have species observations (e.g. transfer to a new location), or onto forecast environmental data - i.e. we could use physical model output and estimate where species will be under forecast conditions. "Nowcasts" are becoming common, where models are produced continuously for every day, and assimilate observed data - Climate projections (e.g. 2050 or 2100) of species distributions are also common, and seasonal to decadal scale forecasting is getting better skill

#### read in example forecast data layers for 2019 (the variables have been derived and renamed to fit the model)
forecast2019 <- read.csv('forecast.csv')

### predict the gam onto forecast data
forecast2019$gam.predict <- predict(gam4, type = "response", newdata = forecast2019)

ggplot() + 
  geom_point(data = forecast2019, aes(x = lon, y = lat, colour = gam.predict)) + 
  scale_colour_gradientn(colours = occurrenceCols(20), name = "P(occur)") + 
  theme_classic() + 
  geom_path(data = world, aes(long, lat, group = group)) + 
  ylim(53, 64) + 
  xlim(-181, -156) 


