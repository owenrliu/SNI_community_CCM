### Understory algae and pavement from benthic cover data
library(dplyr)
W_D <- getwd()

#fix dates in the datasets to be more readable by R
betterDates <-function(dat) {
  dates<-as.Date(as.character(dat$Date),format="%m/%d/%y")
  return(dates)
}

# benthic cover data
#**************************************
benthcover.dat <- read.csv(paste0(W_D,"/data/Benthic cover summary data.csv"))
benthcover.dat$Date <- betterDates(benthcover.dat)

# All species, all surveys
spp.key <- read.csv(file=paste(W_D,'/data/Table4_Species_sampled.csv',sep=''),stringsAsFactors = F)
spp.key.bcov <- spp.key %>% filter(DataSet=="Benthic cover")


# quick function to turn a "genus species" into an abbreviated gen.spe identifier
abbr.species.names <- function(x) {
  temp <- strsplit(x," ")[[1]]
  g <- substr(temp[1],1,3)
  spe <- substr(temp[2],1,3)
  paste(g,spe,sep=".")
}

spp.key.bcov$name <- sapply(spp.key.bcov$SpeciesName,abbr.species.names)

# Same for benthic cover, which has slightly different original data structure
benthcovdat <- benthcover.dat %>%
  rename(site=Station,period=Period,date=Date,dens=CoverMean) %>%
  left_join(select(spp.key.bcov,SpeciesCode,name),by="SpeciesCode") %>%
  
  # remove unneeded columns
  select(-SpeciesCode,-date,-Replicates,-CoverSE) %>%
  rename(spp=name) %>%
  arrange(spp,period)

#*************************************
# There's one more problem, which is that there are some periods which were not sampled. It's better to make those values NAs in the data than just skipping over time periods.  This set of code expands the timeseries for each site/species to include NA density values for missing survey times.

periods.all <- data_frame(period=1:63)


# Filling in cover data
sites <- sort(unique(benthcovdat$site))
spps <- sort(unique(benthcovdat$spp))
benthcovdat.full <- data_frame()
for(i in 1:length(sites)) {
  for(j in 1:length(spps)) {
    temp <- benthcovdat %>% filter(site==sites[i],spp==spps[j])
    temp <- full_join(periods.all,temp,by="period") %>% mutate(site=sites[i],spp=spps[j])
    benthcovdat.full <- bind_rows(benthcovdat.full,temp)
  }
}
rm(temp)

#******************************
# As a last data setup step, we combine the two datasets into one, "longform" data
benthcovdat.full$site <- as.character(benthcovdat.full$site)

##*** Choose and Aggregate Species ***********
# We choose commonly occuring species, grouping them by encrusting coralline algae,
# calcareous red algae (Calliarthron spp., Bossiella spp. and 
# Corallina officinalis), foliose and filamentous red algae
# (Rhodymenia californica, Crytopleura spp., Plocamium cartilagineum, Pterosiphonia spp.,
# Nienburgia andersoniana, Laurencia pacifica), and the seasonal brown seaweeds
# Dictyota binghamiae and Desmarestia ligulata

encrust <- "Enc.cor"
calred <- c("Cal.spp","Bos.spp","Cor.off")
folred <- c("Rho.cal","Cry.spp","Plo.car","Pte.spp","Nie.and","Lau.pac")
brown <- c("Dic.bin","Des.lig")
species <- c(encrust,calred,folred,brown)

categorize_algaes <- function(x) {
  sapply(x,function(y) {
    if(y %in% encrust) return("encrust")
    if(y %in% calred) return("calred")
    if(y %in% folred) return("folred")
    if(y %in% brown) return("brown")
    else return(NA)
  })
}


# normalize the time series
#***********************************
benthcovdat.sel <- benthcovdat.full %>% filter(spp %in% species) %>%
  mutate(group=categorize_algaes(spp)) %>%
  group_by(site,period,group) %>%
  summarise(grp.dens=sum(dens,na.rm=F)) %>%
  rename(spp=group,dens=grp.dens)

dat.norm <- benthcovdat.sel %>%
  group_by(spp,site) %>%
  mutate(norm=(dens-mean(dens,na.rm=T))/sd(dens,na.rm=T)) %>%
  ungroup() %>%
  unique()

write.csv(dat.norm,file=paste0(W_D,"/data/benthcov_trim_normalized",Sys.Date(),".csv"),row.names = F)

## test
library(rEDM)
block <- dat.norm %>%
  select(-dens) %>%
  spread(key=spp,value=norm)

species <- c("folred","calred","encrust","brown")
par(mfrow=c(2,2))

block_segs_sites <- block %>%  
  mutate(ind = row_number()) %>% 
  group_by(site) %>%
  summarise(first=first(ind),last=last(ind))

for(i in 1:length(species)) {
  spp <- species[i]
  dat <- block %>% select(matches(spp)) %>% as.data.frame()
  out <- simplex(as.numeric(dat[,1]),lib=c(127,189),E=2:10,silent=T)
  plot(out$E,out$rho,type="l",main=fullkey$long[match(spp,fullkey$short)],
       xlab="Embedding Dimension (E)",ylab=expression(paste("Skill, ",rho)))
}

navfac.cov <- block %>% filter(site==1)
dutch.cov <- block %>% filter(site %in% c(4,5))
daytona.cov <- block %>% filter(site==6)
westend.cov <- block %>% filter(site %in% c(2,3))
sandy.cov <- block %>% filter(site==7)

plot(sandy.cov$period,sandy.cov$brown, main="brown")
plot(sandy.cov$period,sandy.cov$encrust,main="encrust")
plot(sandy.cov$period,sandy.cov$folred,main="folred")
plot(sandy.cov$period,sandy.cov$calred,main="calred")

plot(navfac.cov$period,navfac.cov$brown, main="brown")
plot(navfac.cov$period,navfac.cov$encrust,main="encrust")
plot(navfac.cov$period,navfac.cov$folred,main="folred")
plot(navfac.cov$period,navfac.cov$calred,main="calred")

plot(westend.cov$period,westend.cov$brown, main="brown")
plot(westend.cov$period,westend.cov$encrust,main="encrust")
plot(westend.cov$period,westend.cov$folred,main="folred")
plot(westend.cov$period,westend.cov$calred,main="calred")

plot(dutch.cov$period,dutch.cov$brown, main="brown")
plot(dutch.cov$period,dutch.cov$encrust,main="encrust")
plot(dutch.cov$period,dutch.cov$folred,main="folred")
plot(dutch.cov$period,dutch.cov$calred,main="calred")

plot(daytona.cov$period,daytona.cov$brown, main="brown")
plot(daytona.cov$period,daytona.cov$encrust,main="encrust")
plot(daytona.cov$period,daytona.cov$folred,main="folred")
plot(daytona.cov$period,daytona.cov$calred,main="calred")
