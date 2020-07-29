#
# read in LST temperature dbf output files
# clean by RCAID, sum length, change column names
# normalize by adding YEAR and BASIN values
# output as csv
#
library(zoo)
library(foreign)
library(tidyverse)
library(dplyr)
library(pals)
library(ggplot2)

path <- "/Users/jonathanarmstrong/Dropbox/For_JA"

basin <- c("ASO", "ENT", "JDR", "LEM", "LOC", "LOL", "MET", "PAH", "PAN", "SFS", "TUC", "UGR", "USL", "WEN")
year <- c("2011", "2012", "2013", "2014", "2015", "2016")

winter <- c(1:10, 39:46)
spring <- c(11:20)
summer <- c(21:31)
hot <- c(24:28)
fall <- c(32:38)
weeks <- c(1:46)
fish.range <- 500
temp_bins <- c(12, 18, 18, 24)

name.dbf <- "8D_Mn.dbf"
name.csv <- "8D_Mn.csv"


#order vs. width
#w = 0.542· e^0.824*order from Downing et al Inland Waters 2012
ordw<-as.data.frame(cbind(1:8,0.54*exp(0.824*1:8)))
names(ordw)<-c("order","width")

# clean all the year x basin mean temp model dbf files downloaded from SFR wiki site

for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    filename <- file.path(path, paste(basin[ii], year[jj], name.dbf, sep="_"))
    if (file.exists(filename)){
      print(filename)
      tmp.dat <- read.csv(filename) #read in dbf
      tmp.dat <- select(tmp.dat, matches("^RCA|^TM")) #use only RCAID and TM. columns
      names(tmp.dat)[names(tmp.dat)=="RCA_ID"] <- "RCAID" #deal with variation in RCAID column name
      names(tmp.dat)[names(tmp.dat)=="rca_id"] <- "RCAID" #deal with variation in RCAID column name
      tmp.dat <- unique(tmp.dat) #remove duplicate rows
      tmp.dat <- arrange(tmp.dat, RCAID) #sort by RCAID
      tmp.dat[, grep("^T", colnames(tmp.dat))] <- round(tmp.dat[, grep("^T", colnames(tmp.dat))], digits=2) #trunc temps to 2 sig fig
      tmp.dat[, grep("^T", colnames(tmp.dat))][tmp.dat[, grep("^T", colnames(tmp.dat))]<0] <- 0 #set negative temp = 0
      newname <- unlist(map(unlist(map(colnames(tmp.dat), str_split, pattern="[1-9][1-9]", n=2), recursive=FALSE), str_c, collapse="")) #drop year reference in TM. column name
      newname <- unlist(map(unlist(map(newname, str_split, pattern="i", n=2), recursive=FALSE), str_c, collapse="")) #fix the TMin v TMn naming
      names(tmp.dat) <- newname #update names
      tmp.dat <- add_column(tmp.dat, BASIN = basin[ii], .after=1) #add BASIN attribute
      tmp.dat <- add_column(tmp.dat, YEAR = year[jj], .after=2) #add YEAR attribute
      outfile <- file.path(path, paste(basin[ii], year[jj], name.csv, sep="_")) #where to write
      write.csv(tmp.dat, outfile, row.names = FALSE) #write out csv
    }
  }
}


# build a single data structure (matrix), highly normalized, with all temp by year and 'shed

for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    filename <- file.path(path, paste(basin[ii], year[jj], name.csv, sep="_"))
    if (file.exists(filename)){
      print(filename)
      tmp.dat <- read.csv(filename) #read in the clean temperature data files year x basin
      if((ii+jj)==2)
        all.dat <- tmp.dat
      else
        all.dat <- rbind(all.dat, tmp.dat)
    }
  }
}

# adjust TMn column names to something slightly more sensible

names(all.dat)[grep("^T", names(all.dat))] <- paste("TMn.", 1:46, sep="")

# Hack to replace three missing weeks of data with adjacent value
#all.dat[(all.dat$BASIN=="JDR" & all.dat$YEAR=="2016"), 'TMn.46'] <- all.dat[(all.dat$BASIN=="JDR" & all.dat$YEAR=="2016"), 'TMn.43']
#all.dat[(all.dat$BASIN=="JDR" & all.dat$YEAR=="2016"), 'TMn.45'] <- all.dat[(all.dat$BASIN=="JDR" & all.dat$YEAR=="2016"), 'TMn.43']
#all.dat[(all.dat$BASIN=="JDR" & all.dat$YEAR=="2016"), 'TMn.44'] <- all.dat[(all.dat$BASIN=="JDR" & all.dat$YEAR=="2016"), 'TMn.43']
#all.dat[(all.dat$BASIN=="UGR" & all.dat$YEAR=="2011"), 'TMn.46'] <- all.dat[(all.dat$BASIN=="UGR" & all.dat$YEAR=="2011"), 'TMn.45']


# save temperature data set - later use readRDS(filename) to load

saveRDS(all.dat, file.path(path, "all.temp.RData"))

# all.dat <- readRDS(file.path(path, "all.temp.RData"))

#build models to estimate weekly growth increment by Pvalue

#Growth [g/g/d] - From WIBioE model runs (params from FB4) for P <- c(0.4, 1.0) 
growth.rate.file <- file.path(path, "SC_OM_growthrate.csv")
grow.data <- read.csv(growth.rate.file, header=TRUE)

#NASA week (8d) growth for a 10g fish
out_grow_OM_040 <- 10.0 * 8 * grow.data$OM_P040
out_grow_OM_100 <- 10.0 * 8 * grow.data$OM_P100
out_grow_SC_040 <- 10.0 * 8 * grow.data$SC_P040
out_grow_SC_100 <- 10.0 * 8 * grow.data$SC_P100

#fit 6th order polynomial to growth(T) data by P-value
in_temp <- grow.data$Temp
est_grow_OM_040 <- lm(out_grow_OM_040 ~ poly(in_temp, 8))
est_grow_OM_100 <- lm(out_grow_OM_100 ~ poly(in_temp, 8))
est_grow_SC_040 <- lm(out_grow_SC_040 ~ poly(in_temp, 8))
est_grow_SC_100 <- lm(out_grow_SC_100 ~ poly(in_temp, 8))


#predict growth for all TMn values over all years and basins
OMP040.names <- sub("TMn", "OMP040", grep("^T", names(all.dat)))
OMP100.names <- sub("TMn", "OMP100", grep("^T", names(all.dat)))
SCP040.names <- sub("TMn", "SCP040", grep("^T", names(all.dat)))
SCP100.names <- sub("TMn", "SCP100", grep("^T", names(all.dat)))

tmp.dat <- select(all.dat, grep("^T", names(all.dat), value=TRUE))

all.dat <- cbind(all.dat, OMP040.names = apply(tmp.dat, 2, function(x){predict.lm(est_grow_OM_040,data.frame(in_temp=x))}))
all.dat <- cbind(all.dat, OMP100.names = apply(tmp.dat, 2, function(x){predict.lm(est_grow_OM_100,data.frame(in_temp=x))}))
all.dat <- cbind(all.dat, SCP040.names = apply(tmp.dat, 2, function(x){predict.lm(est_grow_SC_040,data.frame(in_temp=x))}))
all.dat <- cbind(all.dat, SCP100.names = apply(tmp.dat, 2, function(x){predict.lm(est_grow_SC_100,data.frame(in_temp=x))}))

names(all.dat) <- sub(".names.TMn", "", names(all.dat))

# more.all.dat <- cbind(all.dat, apply(all.dat, grep("^T", names(all.dat), values=TRUE),
#                                       function(x){predict.lm(est_grow_40, data.frame(in_temp = x))}))
# more.all.dat <- cbind(all.dat, apply(all.dat, 2, function(x){predict.lm(est_grow_40,data.frame(in_temp=x))}))

#predict growth for all TMn values over all years and basins
#for (ii in 1:46){
#  all.dat[,49+ii] <- predict.lm(est_grow_40, data.frame(in_temp = all.dat[, 3+ii]))
#  names(all.dat)[49+ii] <- paste("OMP40_", ii, sep="")
#}
#for (ii in 1:46){
#  all.dat[,95+ii] <- predict.lm(est_grow_45, data.frame(in_temp = all.dat[, 3+ii]))
#  names(all.dat)[95+ii] <- paste("OMP45_", ii, sep="")
#}

# save temperature / growth data set - later use readRDS(filename) to load

saveRDS(all.dat, file.path(path, "all.temp.growth.RData"))

#all.dat <- readRDS(file.path(path, "all.temp.growth.RData"))

#import spatial reference data
space.dat.name <- file.path(path, "Basin_Geo_Ref_20191216.csv")
space.dat <- read.csv(space.dat.name, header = TRUE)

all.dat <- merge(all.dat, space.dat, by=c("BASIN", "RCAID"))

space.dat.name <- file.path(path, "Basin_Geo_Ref_Add_SO_20200405.csv")
space.dat <- read.csv(space.dat.name, header = TRUE)

all.dat <- merge(all.dat, space.dat, by=c("BASIN", "RCAID"))

saveRDS(all.dat, file.path(path, "all.temp.growth.spatial.RData"))

#all.dat <- readRDS(file.path(path, "all.temp.growth.spatial.RData"))

#
#
## THIS IS ADDING UP LENGTH ATTRIBUTE
TOT_LEN <- matrix(, nrow=length(basin), ncol=length(year))

for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    TOT_LEN[ii, jj] <- 0
    TOT_LEN[ii, jj] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), "LENGTH")))
  }
}  

##
TOT_Count <- matrix(, nrow=length(basin), ncol=length(year))
row.names(TOT_Count)<-basin
colnames(TOT_Count)<-year
for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    TOT_Count[ii, jj] <- 0
    TOT_Count[ii, jj] <- length(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), "LENGTH")))
  }
}  

#new code from CJ 5-5-20
TOT_LEN_NORM <- data.frame(BASIN=character(), YEAR=character(), LENGTH=double(), stringsAsFactors = FALSE)
cntr=1
for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    TOT_LEN_NORM[cntr,] <- c(basin[ii], year[jj], sum(as.numeric(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), "LENGTH")))))
    cntr=cntr+1
  }
}  
#
#

#summary metrics - growth potential [g] && production [g/m] by reach summed over all weeks


all.dat$OMP040_TOT_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", weeks, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$OMP040_WIN_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", winter, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$OMP040_SPR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", spring, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$OMP040_SUM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", summer, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$OMP040_FAL_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", fall, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$OMP040_TOT_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", weeks, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$OMP040_WIN_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", winter, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$OMP040_SPR_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", spring, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$OMP040_SUM_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", summer, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$OMP040_FAL_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", fall, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$OMP040_TOT_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", weeks, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$OMP040_WIN_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", winter, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$OMP040_SPR_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", spring, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$OMP040_SUM_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", summer, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$OMP040_FAL_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP040.", fall, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$OMP040_TOT_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", weeks, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$OMP040_WIN_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", winter, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$OMP040_SPR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", spring, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$OMP040_SUM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", summer, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$OMP040_FAL_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", fall, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$OMP040_TOT_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", weeks, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$OMP040_WIN_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", winter, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$OMP040_SPR_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", spring, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$OMP040_SUM_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", summer, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$OMP040_FAL_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", fall, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$OMP040_TOT_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", weeks, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$OMP040_WIN_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", winter, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$OMP040_SPR_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", spring, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$OMP040_SUM_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", summer, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$OMP040_FAL_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP040.", fall, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$OMP100_TOT_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", weeks, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$OMP100_WIN_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", winter, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$OMP100_SPR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", spring, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$OMP100_SUM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", summer, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$OMP100_FAL_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", fall, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$OMP100_TOT_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", weeks, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$OMP100_WIN_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", winter, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$OMP100_SPR_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", spring, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$OMP100_SUM_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", summer, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$OMP100_FAL_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", fall, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$OMP100_TOT_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", weeks, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$OMP100_WIN_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", winter, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$OMP100_SPR_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", spring, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$OMP100_SUM_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", summer, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$OMP100_FAL_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("OMP100.", fall, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$OMP100_TOT_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", weeks, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$OMP100_WIN_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", winter, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$OMP100_SPR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", spring, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$OMP100_SUM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", summer, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$OMP100_FAL_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", fall, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$OMP100_TOT_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", weeks, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$OMP100_WIN_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", winter, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$OMP100_SPR_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", spring, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$OMP100_SUM_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", summer, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$OMP100_FAL_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", fall, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$OMP100_TOT_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", weeks, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$OMP100_WIN_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", winter, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$OMP100_SPR_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", spring, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$OMP100_SUM_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", summer, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$OMP100_FAL_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("OMP100.", fall, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$SCP040_TOT_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", weeks, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$SCP040_WIN_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", winter, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$SCP040_SPR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", spring, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$SCP040_SUM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", summer, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$SCP040_FAL_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", fall, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$SCP040_TOT_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", weeks, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$SCP040_WIN_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", winter, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$SCP040_SPR_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", spring, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$SCP040_SUM_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", summer, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$SCP040_FAL_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", fall, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$SCP040_TOT_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", weeks, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$SCP040_WIN_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", winter, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$SCP040_SPR_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", spring, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$SCP040_SUM_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", summer, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$SCP040_FAL_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP040.", fall, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$SCP040_TOT_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", weeks, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$SCP040_WIN_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", winter, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$SCP040_SPR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", spring, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$SCP040_SUM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", summer, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$SCP040_FAL_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", fall, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$SCP040_TOT_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", weeks, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$SCP040_WIN_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", winter, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$SCP040_SPR_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", spring, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$SCP040_SUM_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", summer, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$SCP040_FAL_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", fall, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$SCP040_TOT_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", weeks, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$SCP040_WIN_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", winter, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$SCP040_SPR_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", spring, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$SCP040_SUM_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", summer, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$SCP040_FAL_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP040.", fall, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$SCP100_TOT_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", weeks, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$SCP100_WIN_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", winter, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$SCP100_SPR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", spring, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$SCP100_SUM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", summer, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$SCP100_FAL_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", fall, sep="")])) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$SCP100_TOT_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", weeks, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$SCP100_WIN_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", winter, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$SCP100_SPR_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", spring, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$SCP100_SUM_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", summer, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$SCP100_FAL_SR_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", fall, sep="")])) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$SCP100_TOT_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", weeks, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$SCP100_WIN_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", winter, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$SCP100_SPR_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", spring, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$SCP100_SUM_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", summer, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$SCP100_FAL_SRM_GP <- apply(all.dat, 1, function(x){sum(as.numeric(x[paste("SCP100.", fall, sep="")])) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})

all.dat$SCP100_TOT_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", weeks, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(weeks))})
all.dat$SCP100_WIN_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", winter, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(winter))})
all.dat$SCP100_SPR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", spring, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(spring))})
all.dat$SCP100_SUM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", summer, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(summer))})
all.dat$SCP100_FAL_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", fall, sep="")], na.rm=TRUE)) *
                                                      as.numeric(x["LENGTH"])/(fish.range*length(fall))})

all.dat$SCP100_TOT_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", weeks, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(weeks))})
all.dat$SCP100_WIN_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", winter, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(winter))})
all.dat$SCP100_SPR_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", spring, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(spring))})
all.dat$SCP100_SUM_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", summer, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(summer))})
all.dat$SCP100_FAL_SR_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", fall, sep="")], na.rm=TRUE)) *
                                                          as.numeric(x["OM_SR"])/(fish.range*length(fall))})

all.dat$SCP100_TOT_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", weeks, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(weeks))})
all.dat$SCP100_WIN_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", winter, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(winter))})
all.dat$SCP100_SPR_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", spring, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(spring))})
all.dat$SCP100_SUM_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", summer, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(summer))})
all.dat$SCP100_FAL_SRM_PR <- apply(all.dat, 1, function(x){median(as.numeric(x[paste("SCP100.", fall, sep="")], na.rm=TRUE)) *
                                                          (as.numeric(x["OM_SR"]) + as.numeric(x["OM_RM"]))/(fish.range*length(fall))})



# these lines generate a 0/1 binary variable of whether reaches exceed a temperature threshold averaged across specific weeks
all.dat$MAX_T_13 <- apply(all.dat, 1, function(x){ifelse(mean(as.numeric(x[27:34]))>13, 1, 0)})
all.dat$MAX_T_17 <- apply(all.dat, 1, function(x){ifelse(mean(as.numeric(x[27:34]))>17, 1, 0)})
all.dat$MAX_T_16 <- apply(all.dat, 1, function(x){ifelse(mean(as.numeric(x[27:34]))>16, 1, 0)})
all.dat$MAX_T_18 <- apply(all.dat, 1, function(x){ifelse(mean(as.numeric(x[27:34]))>18, 1, 0)})

#new code from CJ, start...
#Asks what fraction of watershed is bimodal growth
#normalized length metric
junk <- apply(all.dat, 1, function(x){as.numeric(TOT_LEN_NORM$LENGTH[TOT_LEN_NORM$BASIN==as.character(x["BASIN"]) & 
                                                                       TOT_LEN_NORM$YEAR==as.character(x["YEAR"])])})

all.dat$NORMLEN <- all.dat$LENGTH/junk #length of each reach divided by total lenght of watershed, so what fraction is this reach


OMP040.UvB <- t(apply(all.dat, 1, function(x){ifelse((as.numeric(x[paste("OMP040.", weeks, sep="")]) / 
                                                        mean(as.numeric(x[paste("OMP040.", hot, sep="")])) > 1), 
                                                     x["NORMLEN"], 0)}))

OMP040.UvB <- apply(OMP040.UvB, 2, as.numeric)
OMP040.UvB.names <- paste("OMP040.UvB.", 1:46, sep="")
colnames(OMP040.UvB) <- OMP040.UvB.names
all.dat <- cbind(all.dat, OMP040.UvB)

OMP100.UvB <- t(apply(all.dat, 1, function(x){ifelse((as.numeric(x[paste("OMP100.", weeks, sep="")]) / 
                                                        mean(as.numeric(x[paste("OMP100.", hot, sep="")])) > 1), 
                                                     x["NORMLEN"], 0)}))
OMP100.UvB <- apply(OMP100.UvB, 2, as.numeric)
OMP100.UvB.names <- paste("OMP100.UvB.", 1:46, sep="")
colnames(OMP100.UvB) <- OMP100.UvB.names
all.dat <- cbind(all.dat, OMP100.UvB)

SCP040.UvB <- t(apply(all.dat, 1, function(x){ifelse((as.numeric(x[paste("SCP040.", weeks, sep="")]) / 
                                                        mean(as.numeric(x[paste("SCP040.", hot, sep="")])) > 1), 
                                                     x["NORMLEN"], 0)}))
SCP040.UvB <- apply(SCP040.UvB, 2, as.numeric)
SCP040.UvB.names <- paste("SCP040.UvB.", 1:46, sep="")
colnames(SCP040.UvB) <- SCP040.UvB.names
all.dat <- cbind(all.dat, SCP040.UvB)

SCP100.UvB <- t(apply(all.dat, 1, function(x){ifelse((as.numeric(x[paste("SCP100.", weeks, sep="")]) / 
                                                        mean(as.numeric(x[paste("SCP100.", hot, sep="")])) > 1), 
                                                     x["NORMLEN"], 0)}))
SCP100.UvB <- apply(SCP100.UvB, 2, as.numeric)
SCP100.UvB.names <- paste("SCP100.UvB.", 1:46, sep="")
colnames(SCP100.UvB) <- SCP100.UvB.names
all.dat <- cbind(all.dat, SCP100.UvB)



#saveRDS(all.dat, file.path(path, "all.temp.growth.spatial.GPPR.RData"))

#all.dat <- readRDS(file.path(path, "all.temp.growth.spatial.GPPR.RData"))


#
#summarize data over reaches, by week, basin and year
#JA:This aggregates across reaches


Tot_GP_PR_TR_colnames <- c("WEEK", "YEAR", "BASIN", 
                          "OM040GP", "OM100GP", "SC040GP", "SC100GP", 
                          "OM040PR", "OM100PR", "SC040PR", "SC100PR",
                          "OM040GP_OMSR", "OM100GP_OMSR", "SC040GP_OMSR", "SC100GP_OMSR", 
                          "OM040PR_OMSR", "OM100PR_OMSR", "SC040PR_OMSR", "SC100PR_OMSR",
                          "OM040GP_OMSRM", "OM100GP_OMSRM", "SC040GP_OMSRM", "SC100GP_OMSRM", 
                          "OM040PR_OMSRM", "OM100PR_OMSRM", "SC040PR_OMSRM", "SC100PR_OMSRM",
                          "TRSU", "TRHT", "TRGR", "TRMD",
                          "OM040GPSU13", "OM100GPSU13", "SC040GPSU13", "SC100GPSU13", 
                          "OM040GPSU17", "OM100GPSU17", "SC040GPSU17", "SC100GPSU17",
                          "OM040GPSU13_OMSR", "OM100GPSU13_OMSR", "SC040GPSU13_OMSR", "SC100GPSU13_OMSR", 
                          "OM040GPSU17_OMSR", "OM100GPSU17_OMSR", "SC040GPSU17_OMSR", "SC100GPSU17_OMSR",
                          "OM040GPSU13_OMSRM", "OM100GPSU13_OMSRM", "SC040GPSU13_OMSRM", "SC100GPSU13_OMSRM", 
                          "OM040GPSU17_OMSRM", "OM100GPSU17_OMSRM", "SC040GPSU17_OMSRM", "SC100GPSU17_OMSRM",
                          "OM040GPposSU13", "OM100GPposSU13", "SC040GPposSU13", "SC100GPposSU13",
                          "OM040GPposSU17", "OM100GPposSU17", "SC040GPposSU17", "SC100GPposSU17",
                          "OM040GPposSU13_OMSR", "OM100GPposSU13_OMSR", "SC040GPposSU13_OMSR", "SC100GPposSU13_OMSR",
                          "OM040GPposSU17_OMSR", "OM100GPposSU17_OMSR", "SC040GPposSU17_OMSR", "SC100GPposSU17_OMSR",
                          "OM040GPposSU13_OMSRM", "OM100GPposSU13_OMSRM", "SC040GPposSU13_OMSRM", "SC100GPposSU13_OMSRM",
                          "OM040GPposSU17_OMSRM", "OM100GPposSU17_OMSRM", "SC040GPposSU17_OMSRM", "SC100GPposSU17_OMSRM",
                          "OPTLEN_10_16", "OPTLEN_12_18", 
                          "OPTLEN_10_16_OMSR", "OPTLEN_12_18_OMSR",
                          "OPTLEN_10_16_OMSRM", "OPTLEN_12_18_OMSRM",
                          "OPTLEN_10_16_SO<4", "OPTLEN_12_18_SO<4", 
                          "OPTLEN_10_16_OMSR_SO<4", "OPTLEN_12_18_OMSR_SO<4",
                          "OPTLEN_10_16_OMSRM_SO<4", "OPTLEN_12_18_OMSRM_SO<4",
                          "OPTLEN_10_16_SO>4", "OPTLEN_12_18_SO>4", 
                          "OPTLEN_10_16_OMSR_SO>4", "OPTLEN_12_18_OMSR_SO>4",
                          "OPTLEN_10_16_OMSRM_SO>4", "OPTLEN_12_18_OMSRM_SO>4",
                          "OMP040.UvB", "OMP100.UvB", "SCP040.UvB", "SC100.UvB", #last two rows of names are new from CJ
                          "OM040GP_m2", "OM100GP_m2", "SC040GP_m2", "SC100GP_m2")

Tot_GP_PR_TR <- matrix(, nrow=length(basin)*length(year)*length(weeks), ncol=length(Tot_GP_PR_TR_colnames))

colnames(Tot_GP_PR_TR) <- Tot_GP_PR_TR_colnames

cntr = 0


for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    if(dim(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]))[1]>0){
      print(c(basin[ii], year[jj]))
      for (kk in 1:length(weeks)){
        #print(kk)
        cntr = cntr+1
# "WEEK", "YEAR", "BASIN"        
        Tot_GP_PR_TR[cntr,1] <- kk
        Tot_GP_PR_TR[cntr,2] <- jj
        Tot_GP_PR_TR[cntr,3] <- ii
# "OM040GP", "OM100GP", "SC040GP", "SC100GP"
        Tot_GP_PR_TR[cntr,4] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) *
                                    unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range) 
        Tot_GP_PR_TR[cntr,5] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) *
                                    unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range) 
        Tot_GP_PR_TR[cntr,6] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) *
                                    unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range) 
        Tot_GP_PR_TR[cntr,7] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) *
                                    unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
# "OM040PR", "OM100PR", "SC040PR", "SC100PR"
        Tot_GP_PR_TR[cntr,8] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) *
                                       unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")/fish.range))        
        Tot_GP_PR_TR[cntr,9] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) *
                                       unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")/fish.range))        
        Tot_GP_PR_TR[cntr,10] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) *
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")/fish.range))        
        Tot_GP_PR_TR[cntr,11] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) *
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")/fish.range))        
# "OM040GP_OMSR", "OM100GP_OMSR", "SC040GP_OMSR", "SC100GP_OMSR"        
        Tot_GP_PR_TR[cntr,12] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range) 
        Tot_GP_PR_TR[cntr,13] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range) 
        Tot_GP_PR_TR[cntr,14] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range) 
        Tot_GP_PR_TR[cntr,15] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
# "OM040PR_OMSR", "OM100PR_OMSR", "SC040PR_OMSR", "SC100PR_OMSR"               
        Tot_GP_PR_TR[cntr,16] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) *
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")/fish.range))        
        Tot_GP_PR_TR[cntr,17] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) *
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")/fish.range))        
        Tot_GP_PR_TR[cntr,18] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) *
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")/fish.range))        
        Tot_GP_PR_TR[cntr,19] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) *
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")/fish.range))        
# "OM040GP_OMSRM", "OM100GP_OMSRM", "SC040GP_OMSRM", "SC100GP_OMSRM"
        Tot_GP_PR_TR[cntr,20] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
        Tot_GP_PR_TR[cntr,21] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
        Tot_GP_PR_TR[cntr,22] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
        Tot_GP_PR_TR[cntr,23] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
# "OM040PR_OMSRM", "OM100PR_OMSRM", "SC040PR_OMSRM", "SC100PR_OMSRM"        
        Tot_GP_PR_TR[cntr,24] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) *
                                       (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
        Tot_GP_PR_TR[cntr,25] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) *
                                       (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
        Tot_GP_PR_TR[cntr,26] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) *
                                       (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
        Tot_GP_PR_TR[cntr,27] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) *
                                       (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                        unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range) 
# "TRSU", "TRHT", "TRGR", "TRMD"        
        Tot_GP_PR_TR[cntr,28] <- sum(ifelse(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("TMn.",kk,sep=""))) > temp_bins[1],
                                            ifelse(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("TMn.",kk,sep=""))) < temp_bins[2],
                                                   unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/
                                                     sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))), 0), 0))
        Tot_GP_PR_TR[cntr,29] <- sum(ifelse(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("TMn.",kk,sep=""))) > temp_bins[3],
                                            ifelse(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("TMn.",kk,sep=""))) < temp_bins[4],
                                                   unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/
                                                     sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))), 0), 0))
        Tot_GP_PR_TR[cntr,30] <- sum(ifelse(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("TMn.",kk,sep=""))) > temp_bins[1],
                                            unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/
                                            sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))), 0))
        Tot_GP_PR_TR[cntr,31] <- median(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("TMn.",kk,sep=""))))
# "OM040GPSU13", "OM100GPSU13", "SC040GPSU13", "SC100GPSU13"
        Tot_GP_PR_TR[cntr,32] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("OMP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,33] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("OMP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,34] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("SCP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,35] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("SCP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
# "OM040GPSU17", "OM100GPSU17", "SC040GPSU17", "SC100GPSU17"        
        Tot_GP_PR_TR[cntr,36] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("OMP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,37] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("OMP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,38] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("SCP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,39] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("SCP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
# "OM040GPSU13_OMSR", "OM100GPSU13_OMSR", "SC040GPSU13_OMSR", "SC100GPSU13_OMSR"        
        Tot_GP_PR_TR[cntr,40] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("OMP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,41] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("OMP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,42] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("SCP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,43] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("SCP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
# "OM040GPSU17_OMSR", "OM100GPSU17_OMSR", "SC040GPSU17_OMSR", "SC100GPSU17_OMSR"        
        Tot_GP_PR_TR[cntr,44] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("OMP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,45] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("OMP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,46] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("SCP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,47] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("SCP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
# "OM040GPSU13_OMSRM", "OM100GPSU13_OMSRM", "SC040GPSU13_OMSRM", "SC100GPSU13_OMSRM", 
        Tot_GP_PR_TR[cntr,48] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("OMP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,49] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("OMP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,50] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("SCP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,51] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1)), paste("SCP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
# "OM040GPSU17_OMSRM", "OM100GPSU17_OMSRM", "SC040GPSU17_OMSRM", "SC100GPSU17_OMSRM",
        Tot_GP_PR_TR[cntr,52] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("OMP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,53] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("OMP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,54] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("SCP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,55] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1)), paste("SCP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
# "OM040GPposSU13", "OM100GPposSU13", "SC040GPposSU13", "SC100GPposSU13",
        Tot_GP_PR_TR[cntr,56] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("OMP040.",kk,sep=""))>0)), 
                                                           paste("OMP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,57] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("OMP100.",kk,sep=""))>0)), 
                                                           paste("OMP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,58] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("SCP040.",kk,sep=""))>0)), 
                                                           paste("SCP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,59] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("SCP100.",kk,sep=""))>0)), 
                                                           paste("SCP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
# "OM040GPposSU17", "OM100GPposSU17", "SC040GPposSU17", "SC100GPposSU17",
        Tot_GP_PR_TR[cntr,60] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("OMP040.",kk,sep=""))>0)), 
                                                           paste("OMP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,61] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("OMP100.",kk,sep=""))>0)), 
                                                           paste("OMP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,62] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("SCP040.",kk,sep=""))>0)), 
                                                           paste("SCP040.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
        Tot_GP_PR_TR[cntr,63] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("SCP100.",kk,sep=""))>0)), 
                                                           paste("SCP100.",kk,sep=""))) * 
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH"))/fish.range)
#"OM040GPposSU13_OMSR", "OM100GPposSU13_OMSR", "SC040GPposSU13_OMSR", "SC100GPposSU13_OMSR",
        Tot_GP_PR_TR[cntr,64] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("OMP040.",kk,sep=""))>0)), 
                                                   paste("OMP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,65] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("OMP100.",kk,sep=""))>0)), 
                                                   paste("OMP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,66] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("SCP040.",kk,sep=""))>0)), 
                                                   paste("SCP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,67] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("SCP100.",kk,sep=""))>0)), 
                                                   paste("SCP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
# "OM040GPposSU17_OMSR", "OM100GPposSU17_OMSR", "SC040GPposSU17_OMSR", "SC100GPposSU17_OMSR"        
        Tot_GP_PR_TR[cntr,68] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("OMP040.",kk,sep=""))>0)),
                                                   paste("OMP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,69] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("OMP100.",kk,sep=""))>0)), 
                                                   paste("OMP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,70] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("SCP040.",kk,sep=""))>0)), 
                                                   paste("SCP040.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
        Tot_GP_PR_TR[cntr,71] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("SCP100.",kk,sep=""))>0)), 
                                                   paste("SCP100.",kk,sep=""))) *
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR"))/fish.range)
# "OM040GPposSU13_OMSRM", "OM100GPposSU13_OMSRM", "SC040GPposSU13_OMSRM", "SC100GPposSU13_OMSRM",
        Tot_GP_PR_TR[cntr,72] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("OMP040.",kk,sep=""))>0)), 
                                                   paste("OMP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,73] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("OMP100.",kk,sep=""))>0)), 
                                                   paste("OMP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,74] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("SCP040.",kk,sep=""))>0)), 
                                                   paste("SCP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,75] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_13==1 & !!sym(paste("SCP100.",kk,sep=""))>0)), 
                                                   paste("SCP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
# "OM040GPposSU17_OMSRM", "OM100GPposSU17_OMSRM", "SC040GPposSU17_OMSRM", "SC100GPposSU17_OMSRM",
        Tot_GP_PR_TR[cntr,76] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("OMP040.",kk,sep=""))>0)),
                                                   paste("OMP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,77] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("OMP100.",kk,sep=""))>0)), 
                                                   paste("OMP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,78] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("SCP040.",kk,sep=""))>0)), 
                                                   paste("SCP040.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
        Tot_GP_PR_TR[cntr,79] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & MAX_T_17==1 & !!sym(paste("SCP100.",kk,sep=""))>0)), 
                                                   paste("SCP100.",kk,sep=""))) *
                                    (unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_SR")) +  
                                     unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "OM_RM")))/fish.range)
# "OPTLEN_10_16", "OPTLEN_12_18"        
        Tot_GP_PR_TR[cntr,80] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                                  !!sym(paste("TMn.",kk,sep=""))>10 & !!sym(paste("TMn.",kk,sep=""))<16), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,81] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                                  !!sym(paste("TMn.",kk,sep=""))>12 & !!sym(paste("TMn.",kk,sep=""))<18), "LENGTH")))/TOT_LEN[ii,jj]
# "OPTLEN_10_16_OMSR", "OPTLEN_12_18_OMSR"        
        Tot_GP_PR_TR[cntr,82] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                                  !!sym(paste("TMn.",kk,sep=""))>10 & !!sym(paste("TMn.",kk,sep=""))<16 & OM_SR>0), 
                                                                  "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,83] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                                  !!sym(paste("TMn.",kk,sep=""))>12 & !!sym(paste("TMn.",kk,sep=""))<18 & OM_SR>0), 
                                                                  "LENGTH")))/TOT_LEN[ii,jj]
# "OPTLEN_10_16_OMSRM", "OPTLEN_12_18_OMSRM"
        Tot_GP_PR_TR[cntr,84] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & !!sym(paste("TMn.",kk,sep=""))>10 & 
                                                                  !!sym(paste("TMn.",kk,sep=""))<16 & sum(OM_RM, OM_SR)>0), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,85] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & !!sym(paste("TMn.",kk,sep=""))>12 & 
                                                                  !!sym(paste("TMn.",kk,sep=""))<18 & sum(OM_RM, OM_SR)>0), "LENGTH")))/TOT_LEN[ii,jj]

# "OPTLEN_10_16_SO<4", "OPTLEN_12_18_SO<4"        
        Tot_GP_PR_TR[cntr,86] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>10 & !!sym(paste("TMn.",kk,sep=""))<16 &
                                                            ST_ORD<=4), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,87] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>12 & !!sym(paste("TMn.",kk,sep=""))<18 &
                                                            ST_ORD<=4), "LENGTH")))/TOT_LEN[ii,jj]
# "OPTLEN_10_16_OMSR_SO<4", "OPTLEN_12_18_OMSR_SO<4"        
        Tot_GP_PR_TR[cntr,88] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>10 & !!sym(paste("TMn.",kk,sep=""))<16 & OM_SR>0 &
                                                            ST_ORD<=4), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,89] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>12 & !!sym(paste("TMn.",kk,sep=""))<18 & OM_SR>0 &
                                                            ST_ORD<=4), "LENGTH")))/TOT_LEN[ii,jj]
 # "OPTLEN_10_16_OMSRM_SO<4", "OPTLEN_12_18_OMSRM_SO<4"
        Tot_GP_PR_TR[cntr,90] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & !!sym(paste("TMn.",kk,sep=""))>10 & 
                                                            !!sym(paste("TMn.",kk,sep=""))<16 & sum(OM_RM, OM_SR)>0 &
                                                            ST_ORD<=4), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,91] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & !!sym(paste("TMn.",kk,sep=""))>12 & 
                                                            !!sym(paste("TMn.",kk,sep=""))<18 & sum(OM_RM, OM_SR)>0 &
                                                            ST_ORD<=4), "LENGTH")))/TOT_LEN[ii,jj]

# "OPTLEN_10_16_SO>4", "OPTLEN_12_18_SO>4"        
        Tot_GP_PR_TR[cntr,92] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>10 & !!sym(paste("TMn.",kk,sep=""))<16 &
                                                            ST_ORD>4), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,93] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>12 & !!sym(paste("TMn.",kk,sep=""))<18 &
                                                            ST_ORD>4), "LENGTH")))/TOT_LEN[ii,jj]
# "OPTLEN_10_16_OMSR_SO>4", "OPTLEN_12_18_OMSR_SO>4"        
        Tot_GP_PR_TR[cntr,94] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>10 & !!sym(paste("TMn.",kk,sep=""))<16 & OM_SR>0 &
                                                            ST_ORD>4), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,95] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                            !!sym(paste("TMn.",kk,sep=""))>12 & !!sym(paste("TMn.",kk,sep=""))<18 & OM_SR>0 &
                                                            ST_ORD>4), "LENGTH")))/TOT_LEN[ii,jj]
# "OPTLEN_10_16_OMSRM_SO>4", "OPTLEN_12_18_OMSRM_SO>4"
        Tot_GP_PR_TR[cntr,96] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & !!sym(paste("TMn.",kk,sep=""))>10 & 
                                                            !!sym(paste("TMn.",kk,sep=""))<16 & sum(OM_RM, OM_SR)>0 &
                                                            ST_ORD>4), "LENGTH")))/TOT_LEN[ii,jj]
        Tot_GP_PR_TR[cntr,97] <- sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & !!sym(paste("TMn.",kk,sep=""))>12 & 
                                                            !!sym(paste("TMn.",kk,sep=""))<18 & sum(OM_RM, OM_SR)>0 &
                                                            ST_ORD>4), "LENGTH")))/TOT_LEN[ii,jj]
        
      }
    }
  }
}

saveRDS(Tot_GP_PR_TR, file.path(path, "summary.basin.year.week.GPPR.RData"))

#Tot_GP_PR_TR <- readRDS(file.path(path, "summary.basin.year.week.GPPR.RData"))



Tot_GP_PR_TR_df <- as.data.frame(Tot_GP_PR_TR)

#New code from CJ to make additional vectors of fish production taking area into account, 
#This code is outside of the loops above so that we don't have to run code for 24h to make a new tot_gp_pr_tr object

Tot_GP_PR_TR_temp <- matrix(, nrow=length(basin)*length(year)*length(weeks), ncol=4)

cntr = 0

for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    if(dim(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]))[1]>0){
      print(c(basin[ii], year[jj]))
      for (kk in 1:length(weeks)){
        #print(kk)
        cntr = cntr+1
        
        Tot_GP_PR_TR_temp[cntr,1] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP040.",kk,sep=""))) * 
                                           unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")) * 
                                           (0.54 * exp(0.824 * unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "ST_ORD"))))) 
        Tot_GP_PR_TR_temp[cntr,2] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("OMP100.",kk,sep=""))) * 
                                           unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")) *  
                                           (0.54 * exp(0.824 * unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "ST_ORD")))))
        Tot_GP_PR_TR_temp[cntr,3] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP040.",kk,sep=""))) * 
                                           unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")) *  
                                           (0.54 * exp(0.824 * unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "ST_ORD")))))
        Tot_GP_PR_TR_temp[cntr,4] <- sum(unlist(select(na.omit(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj])), paste("SCP100.",kk,sep=""))) * 
                                           unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "LENGTH")) *
                                           (0.54 * exp(0.824 * unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]), "ST_ORD")))))
        
      }
    }
  }
}




Tot_GP_PR_TR_temp_df <- as.data.frame(Tot_GP_PR_TR_temp)

Tot_GP_PR_TR_df<-cbind(Tot_GP_PR_TR_df,Tot_GP_PR_TR_temp_df)

names(Tot_GP_PR_TR_df)[292:295]<-c("OMP040area","OMP100area","SCP040area","SCP100area")




#This is code to compute data plotted in Figure 4 in manuscript, panel b

#Calculate contribution of seasonally warm habtiat to annual growth potential across riverscape, for 14 basins of Interior Columbia River 
#Years of data vary by river basin, so some basin-year combos are missing
res1<-as.data.frame(matrix(nrow=length(basin)*length(year), ncol=4))

for (ii in 1:length(basin)){
  #jd is just a temporary object to subset data into, not john day
  jd<-na.omit(all.dat[all.dat$BASIN==unique(all.dat$BASIN)[ii],] )#grab basin ii
 
   for (jj in 1:length(year)){
    
jd15<-jd[jd$YEAR==unique(all.dat$YEAR)[jj],] #grab year jj

if(dim(jd15)[1]<1){  #if there are no data, enter NA and go to next iteration
  res1[jj + 6*(ii-1),1] <- basin[ii]
  res1[jj + 6*(ii-1),2]<-year[jj]
  res1[jj + 6*(ii-1),3]<-NA
  res1[jj + 6*(ii-1),4]<-NA
}else{
  
jd15_warm<-jd15[jd15$MAX_T_17==1,]   #grab the warm reaches max_T_thresholdtemp

if(dim(jd15_warm)[1]<1){  #if there are no warm reaches, enter 0.00099 and go to next iteration, use .00099 because functionally zero but allows us to test code
  res1[jj + 6*(ii-1),1]<-basin[ii]
  res1[jj + 6*(ii-1),2]<-year[jj]
  res1[jj + 6*(ii-1),3]<-0.00099
  res1[jj + 6*(ii-1),4]<-0.00099
}else{
  
jd15_warm_growth<-jd15_warm[,188:233] * jd15_warm$LENGTH #now pull out all columns (i.e., the weeks/reaches) that contain growth potential, multiply by length since reaches are not equal length
jd15_warm_growth[jd15_warm_growth<0]<-0  #subset for positive values only

jd15_growth<-jd15[,188:233] *jd15$LENGTH  
jd15_growth[jd15_growth<0]<-0 #subset positive, i.e., ignore negative winter growth potential

propgrowthwarm<-sum(jd15_warm_growth)/sum(jd15_growth)
proplengthwarm<-sum(jd15_warm$LENGTH)/sum(jd15$LENGTH)

res1[jj + 6*(ii-1),1] <- basin[ii]
res1[jj + 6*(ii-1),2]<-year[jj]
res1[jj + 6*(ii-1),3]<-propgrowthwarm #what fraction of fish growth potential * area comes from warm areas
res1[jj + 6*(ii-1),4]<-proplengthwarm #what fraction of network length exceeds temp threshold
}}
}
}
res<-na.omit(res1) #omit basin-year combos that lack data
res1ord<-with(res,reorder(V1,-V3,median)) #sort

dev.off()
#alternative boxplot version of Fig 4b
boxplot(res$V3~res1ord,col=cubicl(14),ylab="Production potential from warm habitat (g*m2/y)", xlab= "Basin")
#Fig 4b panel, reset directory to save plot
postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/scatter.ps",width = 3, height = 3)
ggplot(data=res,aes(x=reorder(V1,-V3),y=V3,color=V2)) +geom_jitter(width=.2,height=0,size=2)+
  theme_classic()+labs(x = "Basin", y="Growth from warm habitat")+labs(colour ="Year")+ theme(axis.text.x = element_text(angle = 90, size=8))
dev.off()
#summary stat
x1<-aggregate(list(gprop=res$V3),by= list(basin=res$V1),mean)
x2<-aggregate(list(gpropsd=res$V3),by= list(basin=res$V1),sd)
mean(x1$gprop)
sd(x1$gprop)

#count years of data 
ag.1<-aggregate(all.dat$OMP040.1,by=list(all.dat$BASIN,all.dat$YEAR),length)
ag.3<-aggregate(ag.1$Group.2,by=list(ag.1$Group.1),length)
mean(ag.3$x)
boxplot((res$V3/res$V4)~res$V1,col=cubicl(14),ylab="GR_warm : Len_warm") #boxplots of prowarmgrowth / propwarmlength -- so standardized by quantity


#Code to plot area of optimal habitat over time
#objects renamed with lower-case a before them

      
aTot_GP_PR_TR_colnames <- c("WEEK", "YEAR", "BASIN", 
                            "len1", "len2", "len3", "len4", 
                            "len5", "len6", "len7", "len8")
aTot_GP_PR_TR <- matrix(, nrow=length(basin)*length(year)*length(weeks), ncol=length(aTot_GP_PR_TR_colnames))



colnames(aTot_GP_PR_TR) <- aTot_GP_PR_TR_colnames

cntr = 0

# Set the bounds of temperature considered optimal here
lo<- 10
hi<- 16

for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    if(dim(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]))[1]>0){
      print(c(basin[ii], year[jj]))
      for (kk in 1:length(weeks)){
        #print(kk)
        cntr = cntr+1
        # "WEEK", "YEAR", "BASIN"        
        aTot_GP_PR_TR[cntr,1] <- kk
        aTot_GP_PR_TR[cntr,2] <- jj
        aTot_GP_PR_TR[cntr,3] <- ii

 # For each stream order, calculate sum of length that is within temperature bounds, then multiply by estimated width
        
aTot_GP_PR_TR[cntr,4] <-ordw$width[1] * sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==1), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,5] <-ordw$width[2] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==2), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,6] <-ordw$width[3] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==3), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,7] <-ordw$width[4] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                  !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==4), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,8] <-ordw$width[5] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==5), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,9] <-ordw$width[6] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==6), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,10] <-ordw$width[7] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==7), "LENGTH")))#/TOT_LEN[ii,jj]
aTot_GP_PR_TR[cntr,11] <-ordw$width[8] *  sum(unlist(select(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj] & 
                                                    !!sym(paste("TMn.",kk,sep=""))>lo & !!sym(paste("TMn.",kk,sep=""))<hi &
                                                    ST_ORD==8), "LENGTH")))#/TOT_LEN[ii,jj]
      }
    }
  }
}

adf <- na.omit(as.data.frame(aTot_GP_PR_TR))
adf[,12]<-apply(adf[,4:6],1,sum) #grouping low versus high order
adf[,13]<-apply(adf[,7:10],1,sum)
adf[,14]<-apply(adf[,4:11],1,sum)
names(adf)[12]<-"lenlo"
names(adf)[13]<-"lenhi"
names(adf)[14]<-"lenall"


###JA: Code to plot area of thermally optimal habitat by stream order, over the year

nbasins<-length(unique(adf$BASIN))

blues<-c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
thick<-seq(1,5,length=7)
#yy<-5
for(bb in 1:nbasins){
  for(yy in 1:3){
  temp<-na.omit(adf[adf$BASIN==bb &adf$YEAR==yy,])
  if(dim(temp)[1]==0){
    plot(c(1,1), main="",type="n",axes=F)} else{
  plot(temp$WEEK,temp$len1/max(temp$len1),type="n",main=paste(basin[bb],year[yy],sep=" "))
  lines(temp$WEEK,temp$len1/max(temp$len1),col=blues[3],lwd=thick[1])
  lines(temp$WEEK,temp$len2/max(temp$len2),col=blues[4],lwd=thick[2])
  lines(temp$WEEK,temp$len3/max(temp$len3),col=blues[5],lwd=thick[3])
  lines(temp$WEEK,temp$len4/max(temp$len4),col=blues[6],lwd=thick[4])
  lines(temp$WEEK,temp$len5/max(temp$len5),col=blues[7],lwd=thick[5])
  lines(temp$WEEK,temp$len6/max(temp$len6),col=blues[8],lwd=thick[6])
  lines(temp$WEEK,temp$len7/max(temp$len7),col=blues[9],lwd=thick[7])
    }
  }
}


# In function form for making single panels for figure, we assembled panels outside of R because their arrangement was complex

plot.expcont.order<-function(basin, year){
    temp<-na.omit(adf[adf$BASIN==basin &adf$YEAR==year,])
    if(dim(temp)[1]==0){
      plot(c(1,1), main="",type="n",axes=F)} else{
        plot(temp$WEEK,temp$len1/max(temp$len1),type="n",main="",xlab="Week", ylab="Proportion",las=1)
        lines(temp$WEEK,temp$len1/max(temp$len1),col=blues[3],lwd=thick[1])
        lines(temp$WEEK,temp$len2/max(temp$len2),col=blues[4],lwd=thick[2])
        lines(temp$WEEK,temp$len3/max(temp$len3),col=blues[5],lwd=thick[3])
        lines(temp$WEEK,temp$len4/max(temp$len4),col=blues[6],lwd=thick[4])
        lines(temp$WEEK,temp$len5/max(temp$len5),col=blues[7],lwd=thick[5])
        lines(temp$WEEK,temp$len6/max(temp$len6),col=blues[8],lwd=thick[6])
        lines(temp$WEEK,temp$len7/max(temp$len7),col=blues[9],lwd=thick[7])
      }

}

plot.expcont.order(12,6)
##not standardized by max
#yy<-5
for(bb in 1:nbasins){
  for(yy in 1:3){
    temp<-na.omit(adf[adf$BASIN==bb &adf$YEAR==yy,])
    if(dim(temp)[1]==0){
      plot(c(1,1), main="",type="n",axes=F)} else{
        plot(temp$WEEK,temp$len7,type="n",main=paste(basin[bb],year[yy],sep=" "),ylim=c(0,max(temp)))
        lines(temp$WEEK,temp$len1,col=blues[3],lwd=thick[1])
        lines(temp$WEEK,temp$len2,col=blues[4],lwd=thick[2])
        lines(temp$WEEK,temp$len3,col=blues[5],lwd=thick[3])
        lines(temp$WEEK,temp$len4,col=blues[6],lwd=thick[4])
        lines(temp$WEEK,temp$len5,col=blues[7],lwd=thick[5])
        lines(temp$WEEK,temp$len6,col=blues[8],lwd=thick[6])
        lines(temp$WEEK,temp$len7,col=blues[9],lwd=thick[7])
      }
  }
}

#Figure 4 panel C
#plot to compare Expansion-contraction in 2021 vs. 2015
plot.gline<-function(basin1,year1,color.1){
  temp<-na.omit(adf[adf$BASIN==basin1 &adf$YEAR==year1,])
  lines(temp$WEEK,temp$lenall/max(temp$lenall),col=color.1,lwd=thick[4])
}
dev.off()
postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/12v15.ps",width = 6, height = 6)
#Mar 21 = week 11
#June21 = week 23
#Sep21 =  week 34
#Dec21 =  week 45

par(mfrow = c(2,1),las=1,
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
cl1<-polychrome(36)[1]
cl2<-polychrome(36)[3]
cl3<-polychrome(36)[6]

x.labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan")

plot(adf$WEEK,adf$lenall/max(adf$lenall),type="n",ylab="",xaxt="n",bty="n",yaxt="n")
     axis(1,at=c(0,seq(4,46,by=8)),labels= rep("",7))
     axis(2,at=seq(0,1,.5))
plot.gline(1,2,cl1); text(8,.9,"ASO",col=cl1)
plot.gline(2,2,cl2); text(8,.8,"ENT",col=cl2)
plot.gline(7,2,cl3); text(8,.7,"MET",col=cl3)
plot(adf$WEEK,adf$lenall/max(adf$lenall),type="n",ylab="",xaxt="n",bty="n",yaxt="n")
 axis(1,at=c(0,seq(4,46,by=8)), labels=c("",x.labels[seq(2,12,2)]))
 axis(2,at=seq(0,1,.5))
plot.gline(1,5,cl1)
plot.gline(2,5,cl2)
plot.gline(7,5,cl3)
text(8,.9,"ASO",col=cl1)
text(8,.8,"ENT",col=cl2)
text(8,.7,"MET",col=cl3)
dev.off()
## For just low order vs. high 
par(mfrow=c(length(basin),length(year)))
par(mar=c(1,1,1,1))
blues<-c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
thick<-seq(1,5,length=7)
#yy<-5
for(bb in 1:nbasins){
  for(yy in 1:length(year)){
    temp<-na.omit(adf[adf$BASIN==bb &adf$YEAR==yy,])
    if(dim(temp)[1]==0){
      plot(c(1,1), main="",type="n",axes=F)} else{
    plot(temp$WEEK,temp$lenlo/max(temp$lenlo),type="n",main=paste(basin[bb],year[yy],sep=" "),axes=F)
    lines(temp$WEEK,temp$lenlo/max(temp$lenlo),col=blues[5],lwd=thick[4])
    lines(temp$WEEK,temp$lenhi/max(temp$lenhi),col=blues[9],lwd=thick[4])}
  
  }
}

##For all habitat
par(mfrow=c(length(basin),length(year)))
par(mar=c(1,1,1,1))
for(bb in 1:14){
  for(yy in 1:length(year)){
    temp<-na.omit(adf[adf$BASIN==bb &adf$YEAR==yy,])
    if(dim(temp)[1]==0){
      plot(c(1,1), main="",type="n",axes=F)} else{
        plot(temp$WEEK,temp$lenall,type="n",main=paste(basin[bb],year[yy],sep=" "),axes=F)
        lines(temp$WEEK,temp$lenall,col=blues[5],lwd=thick[4])
        }
    
  }
}



#plot annual variation in subset of basins
## For just low order vs. high 

blues<-c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
thick<-seq(1,5,length=7)
#yy<-5
plot.expcont.binary<-function(bsn,yr){
    temp<-na.omit(adf[adf$BASIN==bsn &adf$YEAR==yr,])
    if(dim(temp)[1]==0){
      plot(c(1,1), main="",type="n",axes=F)} else{
        plot(temp$WEEK,temp$lenlo/max(temp$lenlo),type="n",main=" ",las=1,xlab="Week",ylab="Proportion of habitat")
        lines(temp$WEEK,temp$lenlo/max(temp$lenlo),col=blues[5],lwd=thick[4])
        lines(temp$WEEK,temp$lenhi/max(temp$lenhi),col=blues[9],lwd=thick[4])}
    
  }
plot.expcont.binary(3,5)



plot.expcont.area<-function(basinx){
  temp<-na.omit(adf[adf$BASIN==basinx,])
plot(temp$WEEK,temp$lenall/max(temp$lenall),type="n",main="",las=1)
  for(yy in 1:length(year)){  
  temp<-na.omit(adf[adf$BASIN==basinx &adf$YEAR==yy,])
  if(dim(temp)[1]>0){
        lines(temp$WEEK,temp$lenall/max(temp$lenall),col=cubicl(6)[yy],lwd=thick[4])
      }
    
  }
}

plot.expcont.area(3)



#plot by week - Optimal Length (OPTLEN_10_16, OPTLEN_12_18 )
dev.off()
linetypes <- c(1,3)
par(cex = 0.6)
par(mar = c(0,0,0,0), oma = c(6,6,1,1))
cntr=0
for (ii in 14:14){
  #for (ii in 4:4){
  #for (jj in 1:length(year)){
  for (jj in 4:4){
    
    
        plot(unlist(select(filter(aTot_GP_PR_TR_df, BASIN==ii, YEAR==jj), WEEK), use.names=FALSE),  
             col=cntr, lty=linetypes[1], lwd=4, axes=TRUE)
        basinvec <- ii
      }
   
  }




#coarsening temporal resolution -- aggregating over time

#Summarize Tot_GP_PR_TR over years 

SumY_Tot_GP_PR_TR_colnames <- c("BASIN", "WEEK",
                                "OPTLEN_10_16", "OPTLEN_12_18",
                                "OPTLEN_10_16_OMSR", "OPTLEN_12_18_OMSR",
                                "OPTLEN_10_16_OMSRM", "OPTLEN_12_18_OMSRM")

SumY_Tot_GP_PR_TR <- matrix(, nrow=length(basin)*length(weeks), ncol=length(SumY_Tot_GP_PR_TR_colnames))

colnames(SumY_Tot_GP_PR_TR) <- SumY_Tot_GP_PR_TR_colnames

cntr = 0

for (ii in 1:length(basin)){
  for (kk in 1:length(weeks)){
    
    cntr = cntr + 1
    
    SumY_Tot_GP_PR_TR[cntr,1] <- ii
    SumY_Tot_GP_PR_TR[cntr,2] <- kk
    SumY_Tot_GP_PR_TR[cntr,3] <- mean(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$WEEK==kk), 'OPTLEN_10_16']))
    SumY_Tot_GP_PR_TR[cntr,4] <- mean(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$WEEK==kk), 'OPTLEN_12_18']))
    SumY_Tot_GP_PR_TR[cntr,5] <- mean(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$WEEK==kk), 'OPTLEN_10_16_OMSR']))
  

  SumY_Tot_GP_PR_TR[cntr,6] <- mean(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$WEEK==kk), 'OPTLEN_12_18_OMSR']))
    SumY_Tot_GP_PR_TR[cntr,7] <- mean(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$WEEK==kk), 'OPTLEN_10_16_OMSRM']))
    SumY_Tot_GP_PR_TR[cntr,8] <- mean(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$WEEK==kk), 'OPTLEN_12_18_OMSRM']))
  }
}

saveRDS(SumY_Tot_GP_PR_TR, file.path(path, "summary.basin.week.GPPR.RData"))

#SumY_Tot_GP_PR_TR <- readRDS(file.path(path, "summary.basin.week.GPPR.RData"))

SumY_Tot_GP_PR_TR_df <- as.data.frame(SumY_Tot_GP_PR_TR)

#
#


#Summarize Tot_GP_PR_TR over Weeks 

basin_year <- 0
for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    basin_year = basin_year+ifelse((dim(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]))[1]>0), 1, 0)
  }
}

SumW_Tot_GP_PR_TR_colnames <- c("BASIN", "YEAR",
                                "OM040GPSU13", "OM100GPSU13", "SC040GPSU13", "SC100GPSU13", 
                                "OM040GPSU17", "OM100GPSU17", "SC040GPSU17", "SC100GPSU17",
                                "OM040GPSU13_OMSR", "OM100GPSU13_OMSR", "SC040GPSU13_OMSR", "SC100GPSU13_OMSR", 
                                "OM040GPSU17_OMSR", "OM100GPSU17_OMSR", "SC040GPSU17_OMSR", "SC100GPSU17_OMSR",
                                "OM040GPSU13_OMSRM", "OM100GPSU13_OMSRM", "SC040GPSU13_OMSRM", "SC100GPSU13_OMSRM", 
                                "OM040GPSU17_OMSRM", "OM100GPSU17_OMSRM", "SC040GPSU17_OMSRM", "SC100GPSU17_OMSRM")
 
SumW_Tot_GP_PR_TR <- matrix(, nrow=basin_year, ncol=length(SumW_Tot_GP_PR_TR_colnames))

colnames(SumW_Tot_GP_PR_TR) <- SumW_Tot_GP_PR_TR_colnames

cntr = 0

for (ii in 1:length(basin)){
  for (jj in 1:length(year)){
    if(dim(filter(all.dat, BASIN==basin[ii] & YEAR==year[jj]))[1]>0){
      
      cntr = cntr + 1
# "BASIN", "YEAR",    
      SumW_Tot_GP_PR_TR[cntr,1] <- ii
      SumW_Tot_GP_PR_TR[cntr,2] <- jj
# "OM040GPSU13", "OM100GPSU13", "SC040GPSU13", "SC100GPSU13",    
      SumW_Tot_GP_PR_TR[cntr,3] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GPSU13']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GP']))
      SumW_Tot_GP_PR_TR[cntr,4] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GPSU13']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GP']))
      SumW_Tot_GP_PR_TR[cntr,5] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GPSU13']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GP']))
      SumW_Tot_GP_PR_TR[cntr,6] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GPSU13']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GP']))
# "OM040GPSU17", "OM100GPSU17", "SC040GPSU17", "SC100GPSU17",    
      SumW_Tot_GP_PR_TR[cntr,7] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GPSU17']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GP']))
      SumW_Tot_GP_PR_TR[cntr,8] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GPSU17']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GP']))
      SumW_Tot_GP_PR_TR[cntr,9] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GPSU17']))/
                                   sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GP']))
      SumW_Tot_GP_PR_TR[cntr,10] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GPSU17']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GP']))
#"OM040GPSU13_OMSR", "OM100GPSU13_OMSR", "SC040GPSU13_OMSR", "SC100GPSU13_OMSR"    
      SumW_Tot_GP_PR_TR[cntr,11] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GPSU13_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GP_OMSR']))
      SumW_Tot_GP_PR_TR[cntr,12] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GPSU13_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GP_OMSR']))
      SumW_Tot_GP_PR_TR[cntr,13] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GPSU13_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GP_OMSR']))
      SumW_Tot_GP_PR_TR[cntr,14] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GPSU13_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GP_OMSR']))
# "OM040GPSU17_OMSR", "OM100GPSU17_OMSR", "SC040GPSU17_OMSR", "SC100GPSU17_OMSR",   
      SumW_Tot_GP_PR_TR[cntr,15] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GPSU17_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GP_OMSR']))
      SumW_Tot_GP_PR_TR[cntr,16] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GPSU17_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GP_OMSR']))
      SumW_Tot_GP_PR_TR[cntr,17] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GPSU17_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GP_OMSR']))
      SumW_Tot_GP_PR_TR[cntr,18] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GPSU17_OMSR']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GP_OMSR']))
# "OM040GPSU13_OMSRM", "OM100GPSU13_OMSRM", "SC040GPSU13_OMSRM", "SC100GPSU13_OMSRM"  
      SumW_Tot_GP_PR_TR[cntr,19] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GPSU13_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GP_OMSRM']))
      SumW_Tot_GP_PR_TR[cntr,20] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GPSU13_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GP_OMSRM']))
      SumW_Tot_GP_PR_TR[cntr,21] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GPSU13_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GP_OMSRM']))
      SumW_Tot_GP_PR_TR[cntr,22] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GPSU13_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GP_OMSRM']))
# "OM040GPSU17_OMSRM", "OM100GPSU17_OMSRM", "SC040GPSU17_OMSRM", "SC100GPSU17_OMSRM"  
      SumW_Tot_GP_PR_TR[cntr,23] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GPSU17_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM040GP_OMSRM']))
      SumW_Tot_GP_PR_TR[cntr,24] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GPSU17_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'OM100GP_OMSRM']))
      SumW_Tot_GP_PR_TR[cntr,25] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GPSU17_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC040GP_OMSRM']))
      SumW_Tot_GP_PR_TR[cntr,26] <- sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GPSU17_OMSRM']))/
                                    sum(na.omit(Tot_GP_PR_TR_df[(Tot_GP_PR_TR_df$BASIN==ii & Tot_GP_PR_TR_df$YEAR==jj), 'SC100GP_OMSRM']))
    }
  }  
}

saveRDS(SumW_Tot_GP_PR_TR, file.path(path, "summary.basin.GPPR.RData"))

#SumW_Tot_GP_PR_TR <- readRDS(file.path(path, "summary.basin.GPPR.RData"))

SumW_Tot_GP_PR_TR_df <- as.data.frame(SumW_Tot_GP_PR_TR)


#



#JBA: Code to explore annual positive growth across landscape



junk2 <- array(0, dim=dim(all.dat)[1])
for (ii in 1:dim(all.dat)[1]){
  for (jj in weeks){
    if ((all.dat[ii, paste("SCP100.",jj,sep="")]<0) &&
        !is.na(all.dat[ii, paste("SCP100.",jj,sep="")])){
      junk2[ii] = junk2[ii] + all.dat[ii, paste("SCP100.",jj,sep="")]
    }
  }
}


all.dat$SCP100_neg <- junk2  


  junk1 <- array(0, dim=dim(all.dat)[1])
  for (ii in 1:dim(all.dat)[1]){
    for (jj in weeks){
      if ((all.dat[ii, paste("SCP100.",jj,sep="")]>0) &&
          !is.na(all.dat[ii, paste("SCP100.",jj,sep="")])){
        junk1[ii] = junk1[ii] + all.dat[ii, paste("SCP100.",jj,sep="")]
      }
    }
  }
  all.dat$SCP100_pos <- junk1 
  #Should be able to add junk to all.dat if you want as a new column
  hist(na.omit(all.dat$OMP100_pos[all.dat$OMP100_pos<500]), breaks = 20)

  junk2 <- array(0, dim=dim(all.dat)[1])
  for (ii in 1:dim(all.dat)[1]){
   
    junk2[ii]<-max(all.dat[ii, paste("TMn.",weeks,sep="")])
      
    }

  all.dat$Tmax <- junk2 

#Figure 4 panel D 
 #plot Cumulative annual G
 postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/TGlines.ps",width = 4, height = 4)
 
 all.dat.d<-na.omit(all.dat[all.dat$Tmax<30 ,])
 plot(all.dat.d$Tmax,all.dat.d$SCP100_pos,type="n",xlim=c(7,27),ylim=c(0,70),las=1,cex=.5, bty="n",
      xlab="Max. temperature", ylab="Cum. pos. growth")
 for (ii in 1:14) {
   
   temp<-na.omit(all.dat[all.dat$BASIN==basin[ii]  &  all.dat$Tmax<30,])
   if(dim(temp)[1]>1){
     xp<-sort(temp$Tmax)
     yp<-temp$SCP100_pos[order(temp$Tmax)]
     loesl<-loess(yp~xp)
     #plot(temp$Tmax,jitter(temp$SCP100_durg50,3), xlim=c(5,27),ylim=c(0,35),col=c(cubicl(14)[as.numeric(temp$BASIN)],alpha=.5),
     #xlab="Max annual temperature", ylab="Weeks > 0.5 Gmax")
     lines(x=xp,y=predict(loesl),lwd=3,col= "gray")#cubicl(14)[as.numeric(temp$BASIN)])
     #points(8,seq(36,69,length=15)[ii],col=cubicl(14)[as.numeric(temp$BASIN)], pch=19,cex=1)
     #text(9.5,seq(36,69,length=15)[ii],col=cubicl(14)[as.numeric(temp$BASIN)],labels=basin[ii],cex=1)
     #lines(smooth.spline(temp$Tmax,temp$SCP100_durg50),col=cubicl(14)[ii],lwd=3)
     #points(temp$Tmax,temp$SCP100_durg50,cex=1,col=cubicl(14)[ii])
   }
 }
 
 xp<-sort(all.dat.d$Tmax)
 yp<-all.dat.d$SCP100_pos[order(all.dat.d$Tmax)]
 loesl3<-loess(yp~xp)
 #plot(all.dat.d$Tmax,jitter(all.dat.d$SCP100_durg50,3), xlim=c(5,27),ylim=c(0,35),col=c(cubicl(14)[as.numeric(all.dat.d$BASIN)],alpha=.5),
 #  xlab="Max annual temperature", ylab="Weeks > 0.5 Gmax")
 lines(x=xp,y=predict(loesl3),lwd=6)
 #points(8,seq(36,69,length=15)[15],pch=19,cex=1)
 #text(9.5,seq(36,69,length=15)[15],cex=1,labels="ALL")
 
 dev.off()
 
library(scales)
library(plotrix)
##PLOT MAPS OF OCT - JUL GROWTH POTENTIAL for Figure 4 Panel a

    gdiffbins<-seq(-2,3.5,by=.5)
    cols2<-parula(length(gdiffbins))
    greenpurp1<-c('#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837')
  
    #Function to generate maps in fig 4a
    map.gdiff<-function(basinm,yearm){
    par(mar=c(.1,.1,.1,.1))
    plot(unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), MIDPNT_X), use.names=FALSE), 
         unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), MIDPNT_Y), use.names=FALSE), 
         axes=FALSE, type="n", ylab="", xlab="")
    
    #box(lty=1)
    gdif<-unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), SCP100.35), use.names=FALSE) - 
                   unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), SCP100.26), use.names=FALSE)
                 
    segments(unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), STARTPNT_X), use.names=FALSE), 
             unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), STARTPNT_Y), use.names=FALSE),
             unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), ENDPNT_X), use.names=FALSE), 
             unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), ENDPNT_Y), use.names=FALSE),
            lwd = unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), ST_ORD),use.names=FALSE)^1.5/2.5,
             col = ifelse(gdif <gdiffbins[1], cols2[1], 
                          ifelse(gdif <gdiffbins[2], cols2[1], 
                                 ifelse(gdif <gdiffbins[3], cols2[1], 
                                        ifelse(gdif <gdiffbins[4], cols2[1], 
                                               ifelse(gdif <gdiffbins[5], cols2[6], 
                                                      ifelse(gdif <gdiffbins[6], cols2[6], 
                                                             ifelse(gdif <gdiffbins[7], cols2[11], 
                                                                    ifelse(gdif <gdiffbins[8], cols2[11], 
                                                                           ifelse(gdif <gdiffbins[9], cols2[11], 
                                                                                  ifelse(gdif <gdiffbins[10], cols2[11], 
                                                                                         ifelse(gdif <gdiffbins[11], cols2[11],cols2[11]))))))))))), 
             ylab="", xlab="")
    
    #mtext(basin[ii], side=3, line=-1, outer=FALSE, adj=1, cex=0.8)
    }
  #Function to make inlaid histogramps for fig 4a
    hist.gdiff<-function(basinm,yearm){
      par(mar=c(3,2,1,1))
      breakseq<-seq(-2, 3.5,by=.5)
     
        if(dim(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]))[1]>0){ # if there are data to plot (have to use year[yearm] here)
          
          gp1<- unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), SCP100.35), use.names=FALSE)  #gp october
          gp2<-unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), SCP100.26), use.names=FALSE) #gp july
          gp3<-gp1-gp2
          gp3 <- na.locf(gp3)
          #breakseq<-seq(-2, 5,by=.5)
          orders<-unlist(select(filter(all.dat, BASIN==basin[basinm] & YEAR==year[yearm]), ST_ORD),use.names=FALSE)
          orders[orders>7]<-7
          widths1<-ordw$width[orders]
          #whyearm<-density(x=gp3,weights = widths1)
          #weighted.hist(x=gp3,w=widths1, plot=T,breaks=breakseq,yaxt="n",ylab="",col=cols2[c(1,1,1,6,6,11,11,11,11,11,11,11,11)])
          hist(gp3,breaks=breakseq,yaxt="n",ylab="",col=cols2[c(1,1,1,6,6,11,11,11,11,11,11,11,11)], main="")
          #return(median(gp3))
        }
    }
    pal.bands(parula(12))
    par(mfrow=c(1,1))
    
    postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/map.pane1.ps",width = 3, height = 3)
    map.gdiff(12,6)
    dev.off()
     postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/map.pane2.ps",width = 3, height = 3)
    map.gdiff(11,6)
    dev.off()
    postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/map.pane3.ps",width = 3, height = 3)
    map.gdiff(14,6)
     dev.off()
     
     postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/hist.pane1.ps",width = 3, height = 3)
     hist.gdiff(12,6)
     dev.off()
     
postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/hist.pane2.ps",width = 3, height = 3)
    hist.gdiff(11,6)
    dev.off()
    
    postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/hist.pane3.ps",width = 3, height = 3)
    hist.gdiff(14,6)
    dev.off()

# Code to plot 3 constrasting lines of Growth potential over time, for inlaid graph in Fig 4a  
#plot 3 contrastring GR's
grl.dat<-all.dat[all.dat$BASIN=="UGR" & all.dat$YEAR==2016,188:233]  #188


# 25 is cold  #461 bimodal
postscript(file="/Users/jonathanarmstrong/Dropbox/OngoingResearch/GrowthRegimes/LST_FIg2Components/GRlines.ps",width = 3, height = 3)
xvals<-1:46
yvals<-grl.dat[461,]
pdat<-cbind(xvals,as.numeric(yvals))
fml<-loess(pdat[,2]~pdat[,1],span=.35)
plot(pdat[,1],pdat[,2],type="n", xlab="",ylab="",cex.lab=.5,bty="n")
lines(pdat[,1],predict(fml),type="l",col=cols2[11],lwd=3)

yvals<-grl.dat[621,]
pdat<-cbind(xvals,as.numeric(yvals))
fml<-loess(pdat[,2]~pdat[,1],span=.35)
lines(pdat[,1],predict(fml),type="l",col=cols2[1],lwd=3)
which.max(grl.dat[,28])

yvals<-grl.dat[48,]
pdat<-cbind(xvals,as.numeric(yvals))
fml<-loess(pdat[,2]~pdat[,1],span=.35)
lines(pdat[,1],predict(fml),type="l",col=cols2[6],lwd=3)
which.max(grl.dat[,28])
dev.off()

#Figure 2 temperature data are extracted using the code below, 
  
#Get growth and temp lines for JDR -- these are for illustrative purposes, so year doesn't really matter, we use 2014
jdr.dat<-all.dat[all.dat$BASIN=="JDR" & all.dat$YEAR==2014,188:233]  
jdr.dat.t<-all.dat[all.dat$BASIN=="JDR" & all.dat$YEAR==2014,4:49]  

#find contrasting temperature patterns based on July temp
which.max(jdr.dat.t[,28])
which.min(jdr.dat.t[,28])

plot(as.numeric(jdr.dat[790,1:45]),type="l",ylim=c(-1,2.5))
lines(as.numeric(jdr.dat[2434,1:45]),lty=2)
lines(as.numeric(jdr.dat[1,1:45]),lty=3)

rb.dat<-as.data.frame(cbind(as.numeric(jdr.dat.t[790,1:45]), as.numeric(jdr.dat.t[2434,1:45]),as.numeric(jdr.dat.t[1,1:45])))
names(rb.dat)<-c("twarm","tcold","tcool")
rb.dat.g<-as.data.frame(cbind(as.numeric(jdr.dat[790,1:45]), as.numeric(jdr.dat[2434,1:45]),as.numeric(jdr.dat[1,1:45])))
names(rb.dat.g)<-c("gwarm","gcold","gcool")
rb.dat.all<-cbind(rb.dat,rb.dat.g)
write.csv(rb.dat.all,file="JDR_temp_growth2014.csv")



    
#Figure 3 maps of growth potential
#function form for generating panels, which are assembled outside of R

greenpurp<-parula(7)
plot.gmap<-function(pbasin,pyear,pweek){
  

      plot(unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), MIDPNT_X), use.names=FALSE), 
           unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), MIDPNT_Y), use.names=FALSE), 
           
           axes=FALSE, type="n", ylab="", xlab="") 
      #box(lty=1)     
      
      
      segments(unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), STARTPNT_X), use.names=FALSE), 
               unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), STARTPNT_Y), use.names=FALSE),
               unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), ENDPNT_X), use.names=FALSE), 
               unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), ENDPNT_Y), use.names=FALSE),
               
               
               
               col = ifelse(unlist(select(filter(all.dat, BASIN==basin[pbasin], YEAR==year[pyear]), 
                                          paste("SCP100.",pweek,sep="")), use.names=FALSE) < 0, greenpurp[1] , 
                            ifelse(unlist(select(filter(all.dat, BASIN==basin[pbasin], YEAR==year[pyear]), 
                                                 paste("SCP100.",pweek,sep="")), use.names=FALSE) < .4, greenpurp[2] ,
                                   ifelse(unlist(select(filter(all.dat, BASIN==basin[pbasin], YEAR==year[pyear]), 
                                                        paste("SCP100.",pweek,sep="")), use.names=FALSE) < .8, greenpurp[3] ,
                                          ifelse(unlist(select(filter(all.dat, BASIN==basin[pbasin], YEAR==year[pyear]), 
                                                               paste("SCP100.",pweek,sep="")), use.names=FALSE) < 1.2, greenpurp[4] , 
                                                 ifelse(unlist(select(filter(all.dat, BASIN==basin[pbasin], YEAR==year[pyear]), 
                                                                      paste("SCP100.",pweek,sep="")), use.names=FALSE) < 1.6, greenpurp[5] , 
                                                        ifelse(unlist(select(filter(all.dat, BASIN==basin[pbasin], YEAR==year[pyear]), 
                                                                             paste("SCP100.",pweek,sep="")), use.names=FALSE) < 2, greenpurp[6] , 
                                                               greenpurp[7])))))),
               
               lwd = unlist(select(filter(all.dat, BASIN==basin[pbasin] & YEAR==year[pyear]), ST_ORD),use.names=FALSE)^1.7/5,
               ylab="", xlab="")
      
      
}

plot.expcont.order<-function(basin, year){
  temp<-na.omit(adf[adf$BASIN==basin &adf$YEAR==year,])
  if(dim(temp)[1]==0){
    plot(c(1,1), main="",type="n",axes=F)} else{
      plot(temp$WEEK,temp$len1/max(temp$len1),type="n",main="",xlab="Week", ylab="Proportion",las=1)
      lines(temp$WEEK,temp$len1/max(temp$len1),col=blues[3],lwd=thick[1])
      lines(temp$WEEK,temp$len2/max(temp$len2),col=blues[4],lwd=thick[2])
      lines(temp$WEEK,temp$len3/max(temp$len3),col=blues[5],lwd=thick[3])
      lines(temp$WEEK,temp$len4/max(temp$len4),col=blues[6],lwd=thick[4])
      lines(temp$WEEK,temp$len5/max(temp$len5),col=blues[7],lwd=thick[5])
      lines(temp$WEEK,temp$len6/max(temp$len6),col=blues[8],lwd=thick[6])
      lines(temp$WEEK,temp$len7/max(temp$len7),col=blues[9],lwd=thick[7])
    }
  
}

plot.expcont.area<-function(basinx){
  temp<-na.omit(adf[adf$BASIN==basinx,])
  plot(temp$WEEK,temp$lenall/max(temp$lenall),type="n",main="", xlab= "week", las=1,bty="n",ylab="")
  for(yy in 1:length(year)){  
    temp<-na.omit(adf[adf$BASIN==basinx &adf$YEAR==yy,])
    if(dim(temp)[1]>0){
      lines(temp$WEEK,temp$lenall/max(temp$lenall),col=cubicl(6)[yy],lwd=thick[4])
    }
    
  }
}
#Note both .area and .order and .binary all show area of habitat standardized, not length
blues<-c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
thick<-seq(1,5,length=7)
#yy<-5
plot.expcont.binary<-function(bsn,yr){
  temp<-na.omit(adf[adf$BASIN==bsn &adf$YEAR==yr,])
  if(dim(temp)[1]==0){
    plot(c(1,1), main="",type="n",axes=F)} else{
      plot(temp$WEEK,temp$lenlo/max(temp$lenlo),type="n",main=" ",las=1,xlab="",ylab="",xaxt="n",bty="n")
      lines(temp$WEEK,temp$lenlo/max(temp$lenlo),col=blues[5],lwd=thick[4])
      lines(temp$WEEK,temp$lenhi/max(temp$lenhi),col=blues[9],lwd=thick[4])}
  
}




basin.p<-3
year.p<-6
winter.p<-2
spring.p<-16 
summer.p<-28
fall.p<-39


par(mfrow=c(2,3))

layout(mat = matrix(c(1,2,3,4,5,6), 
                    nrow = 2, 
                    ncol = 3),
       heights = c(1, 1),    # Heights of the two rows
       widths = c(1, 1,1))     # Widths of the two columns
par(mar=c(.5,.5,.5,.5))
plot.gmap(basin.p,year.p,winter.p)
plot.gmap(basin.p,year.p,summer.p)
plot.gmap(basin.p,year.p,spring.p)
plot.gmap(basin.p,year.p,fall.p)
par(mar=c(5,4,2,2))
plot.expcont.binary(basin.p,year.p)
plot.expcont.area(basin.p)
#par(new=F)
par(mar=c(5,2,5,2))
pal.bands(parula,n=7,show.names=F)
#labels = c("< 0","0-0.4","0.4-0.8","0.8-0.12","0.12-0.16","0.16-0.2",">0.2")

#Code for interpreting weeks
jdate<-0:364
nweek<-rep(1:46,each=8)
as.Date(jdate, origin = "2016-01-01")

ddd<-as.data.frame(matrix(nrow=365,ncol=3))
names(ddd)<-c("doy","date","nasa.week")
ddd[,1]<-jdate
ddd[,2]<-as.Date(jdate, origin = "2016-01-01")
ddd[,3]<-nweek[1:365]
names(ddd)<-c("doy","date","nasa.week")


#Figure 2 circular plot in center of figure (assembled outside of R)

library(tidyverse)

# Optional code for shading bars by season (not used in manuscript)
winter <- c(1:10, 39:46)
spring <- c(11:20)
summer <- c(21:31)
hot <- c(24:28)
fall <- c(32:38)
areaopt<-na.omit(adf[adf$BASIN==3 &adf$YEAR==6,])

for(i in 1:46){
  ifelse(areaopt$WEEK[i]%in%winter,areaopt$season[i]<-"winter", 
         ifelse(areaopt$WEEK[i]%in%spring,areaopt$season[i]<-"spring", 
                ifelse(areaopt$WEEK[i]%in%summer,areaopt$season[i]<-"summer",areaopt$season[i]<- "fall")))
                       
}

areaopt$lenall
areaopt$lenlo+areaopt$lenhi
plot(areaopt$WEEK,areaopt$lenhi,col="red")
points(areaopt$WEEK,areaopt$lenlo,col="blue")
#reorganize data for stacked chart
sbch<-as.data.frame(cbind(c(areaopt$len1+areaopt$len2+areaopt$len3+areaopt$len4+areaopt$len5,
                                      areaopt$len6+areaopt$len7+areaopt$len8),
                          c(rep("lo",46),rep("hi",46)), 
                          c(rep(1:46,2))), stringsAsFactors=F)

names(sbch)<-c("area","order","week")
sbch$area<-as.numeric(sbch$area)
sbch$week<-as.numeric(sbch$week)
sbch$season<-c(rep("winter",10),rep("spring",11),rep("summer",11),rep("fall",11),rep("winter",3))
#Stacked circular bar chart
thick.l<-1.5
p <- ggplot(sbch, aes(x=week, y=area, fill = order)) +      
  
  
  geom_bar(stat="identity") +  scale_fill_grey(start=0.5, end=0.2) + geom_text(aes(label=week))+
  geom_vline(xintercept = 10.5,lwd=thick.l)+geom_vline(xintercept = 21.5,lwd=thick.l)+geom_vline(xintercept = 33.5,lwd=thick.l)+geom_vline(xintercept = 44.5,lwd=thick.l)+
ylim(-5e7,170000000) +
  
  
  theme_minimal()+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")    
  ) +
  
  coord_polar(start = 10) #rotates plot
p

#Normal circular bar chart
# Make the plot
p2 <- ggplot(areaopt, aes(x=as.factor(WEEK), y=lenall, fill = as.factor(season))) +      
  geom_bar(stat="identity") + scale_fill_grey()+
 
  ylim(-5e7,170000000) +
  
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")   
  ) +
  
 
  coord_polar(start = 0)
p2

