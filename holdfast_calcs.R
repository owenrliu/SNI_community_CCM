## Holdfast Diameter data
library(dplyr)
library(ggplot2)
library(tidyr)

W_D <- getwd()
holdfasts <- read.csv(paste0(W_D,"/data/Giant kelp size frequency.csv"),stringsAsFactors = F)
sum(!is.na(holdfasts$StipeCount))
sum(!is.na(holdfasts$HoldFast))

par(mfrow=c(1,1))
plot(holdfasts$HoldFast,holdfasts$StipeCount,xlab="Holdfast Diameter",ylab="Stipe Count",main="")

holdfasts <- holdfasts %>% mutate(season=ifelse(Period %% 2 == 0, "Spring","Fall")) %>%
  unite(site,Station,Swath)

ggplot(holdfasts,aes(x=HoldFast,y=StipeCount,color=season)) +
  geom_point()

# Mac density data
mac_density <- read.csv(paste0(W_D,"/data/benth_data_trim_normalized2016-09-09.csv"),stringsAsFactors = F) %>%
  filter(spp=="mac")

# Join by period/site
holdfasts <- holdfasts %>% left_join(mac_density,by=c("Period"="period","site"="site")) %>%
  filter(StipeCount != 0)

# Stipe density, spring vs. fall
ggplot(holdfasts,aes(x=season,y=StipeCount,fill=season))+
  geom_boxplot()

ggplot(holdfasts,aes(x=season,y=dens,fill=season))+
  geom_boxplot()

# Relationship between density and stipes and holdfasts
ggplot(holdfasts,aes(x=norm,y=StipeCount,color=season)) +
  geom_point()+
  xlab("Density")+
  ylab("Stipe Count")+
  geom_smooth(formula=log(y)~x)

ggplot(holdfasts,aes(x=dens,y=HoldFast,color=season)) +
  geom_point()+
  xlab("Density")+
  ylab("Holdfast Diameter")+
  geom_smooth(method="lm")
