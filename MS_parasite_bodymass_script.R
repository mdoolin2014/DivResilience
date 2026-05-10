#Parasite and host body mass data analysis


##### Parasitism #####

library(ggplot2)
library(ggpubr)
library(ggbreak)
library(ggbeeswarm)
library(gridExtra)
library(reshape2)
library(AICcmodavg)
library(MASS)


#Finding mean and sd of each group
mean(data142hi$total_worms)
sd(data142hi$total_worms)

library(MASS)
mod <- glm.nb(total_worms ~ microb, data=infecteddat)
summary(mod)


#Plotting worm burden
plot1 <- ggplot(data=infecteddat, aes(microb, total_worms)) +
  geom_boxplot(aes(colour=microb), outlier.shape=NA) + 
  geom_beeswarm(aes(shape=microb, color=microb), cex=4, size=3, alpha=0.9) +
  #geom_point(aes(shape=sex, fill=microb), position=position_jitter(0.1), size=10, alpha=0.9) +
  scale_colour_manual(values=c("#C66EF5","#FDD964")) +
  #scale_fill_manual(values=c("#C66EF5","#FDD964")) +
  scale_shape_manual(values=c(16, 15)) + #With real shapes for sex. 
  ggtitle("All worm counts, D14") +
  theme_bw() + 
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=11),
        legend.title=element_text(size=13, face="bold"), 
        legend.text=element_text(size=10), 
        legend.position="top") + 
  labs(title="Worm burden",  
       colour= "Microbiome", 
       shape="Microbiome",
       subtitle="", 
       x="Microbiome group",
       y="Total worm burden") 
#stat_compare_means(method="anova") 
#  stat_compare_means(comparisons=my_comparisons)

##geom_text(mapping = aes(label = mouse), size = 4, vjust = 2)
plot2 <- plot1 + 
  ylim(c(0,130)) +
  scale_y_break(c(40, 120), ticklabels=c(120, 130))

plot2



##### Body mass #####

mem.bod <- lme(mass ~ microb + day_num + treatment, random=~1|mouse, data=dat)
mem.bod
anova(mem.bod)






