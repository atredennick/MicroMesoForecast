######################################################
#### Plot yearly intercepts and slopes to look for
#### influential years by vital rate
####

library(ggplot2); library(plyr); library(reshape2)

####
#### Bring in data -----------------------
####
years <- c(1933:1945)
grow_coef <- read.csv("./growth/growthQuants.csv")
surv_coef <- read.csv("./survival/survivalQuants.csv")
rec_coef <- read.csv("./recruitment/recruitmentQuants.csv")

pdf(file = "yearly_coefs.pdf", onefile = TRUE, width = 8, height = 3)

grow_coef <- grow_coef[-grep("betaSpp", x = grow_coef$X),]
beta_id <- grep("beta", x = grow_coef$X)
beta_grow <- grow_coef[beta_id,]
beta_grow$year <- rep(years, each=4)
ggplot(beta_grow, aes(x=as.character(year), y=X50., color=species, group=species))+
  geom_line()+
  geom_point()+
  ggtitle("growth / size slope")
int_id <- grep("intYr", x = grow_coef$X)
int_grow <- grow_coef[int_id,]
int_grow$year <- rep(years, each=4)
ggplot(int_grow, aes(x=as.character(year), y=X50., color=species, group=species))+
  geom_line()+
  geom_point()+
  ggtitle("growth / intercept")

surv_coef <- surv_coef[-grep("betaSpp", x = surv_coef$X),]
beta_id <- grep("beta", x = surv_coef$X)
beta_surv <- surv_coef[beta_id,]
beta_surv$year <- rep(years, each=4)
ggplot(beta_surv, aes(x=as.character(year), y=X50., color=species, group=species))+
  geom_line()+
  geom_point()+
  ggtitle("survival / size slope")
int_id <- grep("intYr", x = surv_coef$X)
int_surv <- surv_coef[int_id,]
int_surv$year <- rep(years, each=4)
ggplot(int_surv, aes(x=as.character(year), y=X50., color=species, group=species))+
  geom_line()+
  geom_point()+
  ggtitle("survival / intercept")

int_id <- grep("intYr", x = rec_coef$X)
int_rec <- rec_coef[int_id,]
int_rec$year <- rep(years, each=4)
ggplot(int_rec, aes(x=as.character(year), y=X50., color=species, group=species))+
  geom_line()+
  geom_point()+
  ggtitle("recruitment / intercept")

dev.off()
