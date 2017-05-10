# Office PC directory
setwd("C:/Users/tyatabe/OneDrive/Docs/Grants/Data incubator")
# Laptop
setwd("C:/Users/Tadaishi/SkyDrive/Docs/Grants/Data incubator")

# read in and prepare asthma attacks data
asthma <- read.csv("Asthma_Emergency_Department_Visit_Rates_by_ZIP_Code.csv")

# clean ZIP code data (separate zip from coordinates)
library(zipcode)
asthma$zip <- as.factor(clean.zipcodes(substr(asthma$ZIP.code, 1,5 )))


# Get total emergency visits by zipcode
cases <- asthma[asthma$Age.Group=="All Ages",]
cases <- cases[,-2]

# zipcode coordinates

data(zipcode)
zipcode$zip <- as.factor(zipcode$zip)
cases <- merge(cases, zipcode, by="zip")


# Read in and prepare air quality data
air <- read.csv("openaq.csv")
# each contaminant on a separate column
bc <- air[air$parameter=="bc", c("location", "utc", "value", "latitude", "longitude")]
co <- air[air$parameter=="co", c("location", "utc", "value", "latitude", "longitude")]
no2 <- air[air$parameter=="no2", c("location", "utc", "value", "latitude", "longitude")]
o3 <- air[air$parameter=="o3", c("location", "utc", "value", "latitude", "longitude")]
pm10 <- air[air$parameter=="pm10", c("location", "utc", "value", "latitude", "longitude")]
pm25 <- air[air$parameter=="pm25", c("location", "utc", "value", "latitude", "longitude")]
so2 <- air[air$parameter=="so2", c("location", "utc", "value", "latitude", "longitude")]

# Getting the maximum value at each location (could be other, like the mean)
library(doBy)
bc <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = bc)
co <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = co)
no2 <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = no2)
o3 <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = o3)
pm10 <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = pm10)
pm25 <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = pm25)
so2 <- summaryBy(value ~ location + latitude + longitude, FUN=max, data = so2)




# Getting this into a map
library(raster)
library(rgdal)
# US shapefile
usa <- getData('GADM', country='USA', level=2)
usa <- usa[! usa$NAME_1 %in% c('Alaska','Hawaii'), ]
# California shapefile (county level)
caco <- shapefile("C:/Users/tyatabe/OneDrive/Docs/PhD Epi/Spring 14'/GEO200CN/HW/data/data/California/counties_2000_TA.shp")
# Projecting usa to ca
usa <-spTransform(usa, CRS(projection(caco)))# Not a good projection for display


# Create spatial points data frame for all contaminants
bc.pts <- SpatialPoints(data.frame(bc$longitude, bc$latitude))
bc.pts <- SpatialPointsDataFrame(bc.pts, bc)
projection(bc.pts) <- "+proj=longlat +datum=NAD83"
co.pts <- SpatialPoints(data.frame(co$longitude, co$latitude))
co.pts <- SpatialPointsDataFrame(co.pts, co)
projection(co.pts) <- "+proj=longlat +datum=NAD83"
no2.pts <- SpatialPoints(data.frame(no2$longitude, no2$latitude))
no2.pts <- SpatialPointsDataFrame(no2.pts, no2)
projection(no2.pts) <- "+proj=longlat +datum=NAD83"
o3.pts <- SpatialPoints(data.frame(o3$longitude, o3$latitude))
o3.pts <- SpatialPointsDataFrame(o3.pts, o3)
projection(o3.pts) <- "+proj=longlat +datum=NAD83"
pm10.pts <- SpatialPoints(data.frame(pm10$longitude, pm10$latitude))
pm10.pts <- SpatialPointsDataFrame(pm10.pts, pm10)
projection(pm10.pts) <- "+proj=longlat +datum=NAD83"
pm25.pts <- SpatialPoints(data.frame(pm25$longitude, pm25$latitude))
pm25.pts <- SpatialPointsDataFrame(pm25.pts, pm25)
projection(pm25.pts) <- "+proj=longlat +datum=NAD83"
so2.pts <- SpatialPoints(data.frame(so2$longitude, so2$latitude))
so2.pts <- SpatialPointsDataFrame(so2.pts, so2)
projection(so2.pts) <- "+proj=longlat +datum=NAD83"

## Transform the data to Teale Albers (note that I am using km here!)
TA <- CRS(" +proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0")
bc.pts <- spTransform(bc.pts, TA)
co.pts <- spTransform(co.pts, TA)
no2.pts <- spTransform(no2.pts, TA)
o3.pts <- spTransform(o3.pts, TA)
pm10.pts <- spTransform(pm10.pts, TA)
pm25.pts <- spTransform(pm25.pts, TA)
so2.pts <- spTransform(so2.pts, TA)
# CA and USA to TA
taca <- spTransform(caco, TA)
usa <- spTransform(usa, TA)
# Interpolation (only for o3, pm10 and pm25) through thin plate spline (or universal kriging)
# Add covariates to use for interpolation
## Get elevation data:
elv <- getData('worldclim', res=2.5, var='alt')
r <- raster(usa, res=5)
celv <- projectRaster(elv, r)
celv <- mask(celv, usa)

# Wind speed
wind <- raster("wc2.0_5m_wind_04.tif")
cwind <- projectRaster(wind, r)
cwind <- mask(cwind, usa)

# Mean temp
temp <- raster("wc2.0_2.5m_tavg_04.tif")
ctemp <- projectRaster(temp, r)
ctemp <- mask(ctemp, usa)

# Precipitation
prec <- raster("wc2.0_2.5m_prec_04.tif")
cprec <- projectRaster(prec, r)
cprec <- mask(cprec, usa)

# Radiation
rad <- raster("wc2.0_2.5m_srad_04.tif")
crad <- projectRaster(rad, r)
crad <- mask(crad, usa)

# Pop density...maybe later


# convert cases to spatialpoints data frame
cases.pts <- SpatialPoints(data.frame(cases$longitude, cases$latitude))
cases.pts <- SpatialPointsDataFrame(cases.pts, cases)
projection(cases.pts) <- "+proj=longlat +datum=NAD83"
#Project
cases.pts <- spTransform(cases.pts, TA)

# OZONE
## Will do k-fold cross-validation. Using 20% of the data
set.seed(0)
kf <- sample(1:5, nrow(o3.pts), replace=TRUE)
d <- data.frame(coordinates(o3.pts), o3=o3.pts$value.max)
colnames(d)[1:2] = c('x', 'y')
d$kf <- kf

## Now add elevation. Note that I use the original (high
## spatial resolution) elevation data and lon/lat points. The warning message
## about transformation is given because of the difference in datum between the
## two sources (WGS84 and NAD83)
d$elev <- extract(elv, o3.pts)
d$wind <- extract(wind, o3.pts)
d$temp <- extract(temp, o3.pts)
d$prec <- extract(prec, o3.pts)
d$rad <- extract(rad, o3.pts)


## Interpolation with a thin plate spline with elevation, wind, and temp as predictor variables
library(fields)
rr <- vector(length=5)
for (kf in 1:5) {
  test <- d[d$kf == kf, c('x', 'y', 'elev', 'wind', 'temp', 'o3')]
  train <- d[d$kf != kf, c('x', 'y', 'elev', 'wind', 'temp', 'o3')]
  tps <- Tps(train[,1:5], train$o3)
  p <- predict(tps, test[,1:5])
  rr[kf] <- cor(test$o3, p)^2
}
mean(rr)


## This looks better (not so great though). A prediction	
tps <- Tps(d[, c('x', 'y', 'elev', 'wind', 'temp')], d$o3)
rstack <- stack(celv, cwind, ctemp)
rt <- interpolate(rstack, tps, xyOnly=FALSE)
plot(rt)
plot(o3.pts, add=T, pch=1, cex=o3.pts$value.max*30, col="blue")

#### PM10
set.seed(0)
kf <- sample(1:5, nrow(pm10.pts), replace=TRUE)
d <- data.frame(coordinates(pm10.pts), pm10=pm10.pts$value.max)
colnames(d)[1:2] = c('x', 'y')
d$kf <- kf

d$elev <- extract(elv, pm10.pts)
d$wind <- extract(wind, pm10.pts)
d$temp <- extract(temp, pm10.pts)

## Interpolation...perhaps other variables later
rr <- vector(length=5)
for (kf in 1:5) {
  test <- d[d$kf == kf, c('x', 'y', 'elev', 'wind', 'temp', 'pm10')]
  train <- d[d$kf != kf, c('x', 'y', 'elev', 'wind', 'temp', 'pm10')]
  tps <- Tps(train[,1:5], train$pm10)
  p <- predict(tps, test[,1:5])
  rr[kf] <- cor(test$pm10, p)^2
}
mean(rr)

# A prediction
tps <- Tps(d[, c('x', 'y', 'elev', 'wind', 'temp')], d$pm10)
rstack <- stack(celv, cwind, ctemp)
rt.pm10 <- interpolate(rstack, tps, xyOnly=FALSE)
plot(rt.pm10)
plot(pm10.pts, add=T, pch=1, cex=pm10.pts$value.max/max(pm10.pts$value.max), col="blue")


#### PM25
set.seed(0)
kf <- sample(1:5, nrow(pm25.pts), replace=TRUE)
d <- data.frame(coordinates(pm25.pts), pm25=pm25.pts$value.max)
colnames(d)[1:2] = c('x', 'y')
d$kf <- kf

d$elev <- extract(elv, pm25.pts)
d$wind <- extract(wind, pm25.pts)
d$temp <- extract(temp, pm25.pts)

## Interpolation...perhaps other variables later
rr <- vector(length=5)
for (kf in 1:5) {
  test <- d[d$kf == kf, c('x', 'y', 'elev', 'wind', 'temp', 'pm25')]
  train <- d[d$kf != kf, c('x', 'y', 'elev', 'wind', 'temp', 'pm25')]
  tps <- Tps(train[,1:5], train$pm25)
  p <- predict(tps, test[,1:5])
  rr[kf] <- cor(test$pm25, p)^2
}
mean(rr)

# A prediction
tps <- Tps(d[, c('x', 'y', 'elev', 'wind', 'temp')], d$pm25)
rstack <- stack(celv, cwind, ctemp)
rt.pm25 <- interpolate(rstack, tps, xyOnly=FALSE)
plot(rt.pm25)
plot(pm25.pts, add=T, pch=1, cex=pm25.pts$value.max/max(pm25.pts$value.max), col="blue")


# Extract o3, pm10, amd pm25 values at each zip code where cases data was collected
cases$o3 <- extract(rt, cases.pts)
cases$pm10 <- extract(rt.pm10, cases.pts)
cases$pm25 <- extract(rt.pm25, cases.pts)

plot(cases$o3, cases$Number.of.Visits)
plot(cases$pm10, cases$Number.of.Visits)
plot(cases$pm25, cases$Number.of.Visits)

cases <- cases[!is.na(cases$o3),]

# Modeling cases of asthma
# split the  data into two groups (take out 20% for model evaluation).

## Sampling the data
set.seed(1)
random <- sample(1:nrow(cases), round(nrow(cases)*0.8))
training <-cases[random,]
testing <- cases[-random,]

# trying a random forest
library(randomForest)
rfbf = randomForest(Number.of.Visits ~ o3+pm10+pm25, data=training)
rfbf
#Most important variables
varImpPlot(rfbf)

# Predicting with the testing data
pred <- predict(rfbf, testing[,c("o3", "pm10", "pm25")], type="response",
        norm.votes=F, predict.all=FALSE, proximity=FALSE, nodes=FALSE)

#Model looks good
plot(pred, testing$Number.of.Visits)
cor(pred, testing$Number.of.Visits)^2

# Extract o3, pm10, and pm25 data for all of the USA's zip codes
data(zipcode)
zipcode <- zipcode[!is.na(zipcode$latitude),]
zip.pts <- SpatialPoints(data.frame(zipcode$longitude, zipcode$latitude))
zip.pts <- SpatialPointsDataFrame(zip.pts, zipcode)
projection(zip.pts) <- "+proj=longlat +datum=NAD83"
zip.pts <- spTransform(zip.pts, TA)


zipcode$o3 <- extract(rt, zip.pts)
zipcode$pm10 <- extract(rt.pm10, zip.pts)
zipcode$pm25 <- extract(rt.pm25, zip.pts)
# removing NA's
zipcode <- zipcode[complete.cases(zipcode),]

# Predicting number of cases with the USA estimated contaminats data
zipcode$pred <- predict(rfbf, zipcode[,c("o3", "pm10", "pm25")], type="response",
                norm.votes=F, predict.all=FALSE, proximity=FALSE, nodes=FALSE)

# Projecting for plotting
pred.pts <- SpatialPoints(data.frame(zipcode$longitude, zipcode$latitude))
pred.pts <- SpatialPointsDataFrame(pred.pts, zipcode)
projection(pred.pts) <- "+proj=longlat +datum=NAD83"
pred.pts <- spTransform(pred.pts, TA)

# Plotting
library(choroplethrZip)
colnames(zipcode)[c(1,9)] <- c("region","value")
zip.dat <- zipcode[,c(1,9)]
# US mainland
data(zip.regions)
zip_choropleth(zip.dat, state_zoom = unique(zip.regions$state.name)[-c(18,49)])
# California and New York
zip_choropleth(zip.dat, state_zoom = "california")
zip_choropleth(zip.dat, state_zoom = "new york")
