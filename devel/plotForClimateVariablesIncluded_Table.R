testD <- data.frame(x = rep(seq(1,10),3),
                    y = rep(c(1:3), each=10),
                    z = rbinom(30,1,0.5))

library(ggplot2)

ggplot(testD, hjust = 0.5, vjust = 0.5)+
  geom_raster(aes(x=y, y=x, fill=as.character(z)), color="white")+
  scale_fill_manual(breaks = levels(factor(testD$z)),
                    values = c("grey75", "dodgerblue"))+
  geom_hline(aes(yintercept=x+0.5), color="white")+
  geom_vline(aes(xintercept=y+0.5), color="white")+
  coord_equal()