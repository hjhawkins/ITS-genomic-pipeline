##########################################################################
################### ALPHA DIVERSITY (within sample) ######################
##########################################################################

#you can prune OTUs that are not present in any of the samples - BUT DON'T TRIM MORE THAN THAT! many richness estimates are modeled on singletons and doubletons in the abundance data. You need to leave them in the dataset if you want a meaningful estimate.

### Plot with the plot_richness function from the phyloseq package

p <- plot_richness(phyf2,x = "FH",measures=c("Shannon"), title = "Alpha Diveristy - Shannon") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw() 
p

p <- plot_richness(phyf2,x = "FH",measures=c("Observed"), title = "Alpha Diveristy - Observed") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("Chao1"), title = "Alpha Diveristy - Chao1") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("ACE"), title = "Alpha Diveristy - ACE") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("Simpson"), title = "Alpha Diveristy - Simpson") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p
p <- plot_richness(phyf2,x = "FH",measures=c("InvSimpson"), title = "Alpha Diveristy - InvSimpson") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()
p <- plot_richness(phyf2,x = "FH",measures=c("Fisher"), title = "Alpha Diveristy - Fisher") + facet_grid(zone ~ site) + geom_boxplot() + theme_bw()


pdf(paste0(outDir,"/alpha_diversity_by_treatment.pdf"))
p
dev.off()


# Compute the Shannon diversity associated with each sample and join it with sample annotation.

alpha_Shannon <- estimate_richness(phyf2, split = TRUE, measure = "Shannon")
alpha_Shannon$Sample.ID <- rownames(alpha_Shannon) %>% as.factor()
ps_samp <- Env %>% left_join(alpha_Shannon, by = "Sample.ID") 
#Rem : probl = some samples were removed from sample table

#other way:
est <- estimate_richness(phyf2, split = TRUE, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")) # gives "shannon" for all samples.
shann <- cbind(est,sample_data(phyf2)[,c("fire", "fire.freq", "herbiv", "site", "zone", "site_FH_zone", "site_ttm_zone")]) #add a column "fire " and "herbiv"

write.csv(shann, file = "shannon.csv")
shann <- read.csv("shannon.csv")

alpha_diversity <- estimate_richness(silk_working_control, measure = c("Shannon", "Observed"))
alpha_diversity
H <- alpha_diversity$Shannon
S1 <- alpha_diversity$Observed
S <- log(S1)
evenness <- H/S
evenness
alpha_diversity$Evenness = evenness
alpha_diversity


# plot the calculated shannon

alpha_boxplot <- ggplot(ps_samp, aes(FH, Shannon)) + geom_boxplot()+facet_grid(site ~ Zone)
alpha_boxplot

alpha_boxplot_savanna <- ggplot(subset(ps_samp, site == "N" | site == "L"), aes(FH, Shannon)) +geom_boxplot()+facet_grid(site ~ Zone)
alpha_boxplot_savanna
alpha_boxplot_grassland <- ggplot(subset(ps_samp, site == "U" | site == "B"), aes(FH, Shannon)) +geom_boxplot()+facet_grid(. ~ site)
alpha_boxplot_grassland


# make a table with the average alpha diversity per site, zone and treatment

diversity_means <- ps_samp %>%
  group_by(site, Zone, Site.ttm.zone) %>%
  summarise(mean_div = mean(Shannon)) %>%
  arrange(mean_div)


### test the normality : do an histogram of shannon measures + Kolmogorov-Smirnov test + Shapiro-Wilk test. rem: need to it per group we test, not for the whole dataset.

h<- ggplot(shann, aes(x = Shannon)) + geom_histogram() + theme_bw()
h

ks.test(shann$Shannon, "pnorm", mean = mean(shann$Shannon), sd=sd(shann$Shannon)) #The null hypothesis of the K-S test is that the distribution is normal.Therefore, if p-value of the test is >0.05, we do not reject the null hypothesis and conclude that the distribution in question is not statistically different from a normal distribution.
#D = 0.088895, p-value = 0.09619
#alternative hypothesis: two-sided

shapiro.test(shann$Shannon) # tests the null hypothesis is that the population is normally distributed.have greater power when compared to the K-S test.
#W = 0.92997, p-value = 5.6e-08 => a lot smaller than 0.05 => significantly differ from a normal distribution.


###	Kruskal-Wallis rank sum test : tot test the impact of fire and herbivores.

t_fire <- kruskal.test(shann$Shannon~shann$fire)
t_fire # Kruskal-Wallis chi-squared = 0.10579, df = 1, p-value = 0.745

t_herbiv <- kruskal.test(shann$Shannon~shann$herbiv)
t_herbiv # Kruskal-Wallis chi-squared = 0.076529, df = 1, p-value = 0.7821

### ANOVA
d <- data.frame(lev=shann$site_FH_zone, y=shann$Shannon)

aov <- aov(shann$Shannon~shann$site_FH_zone, data=shann)
summary(aov)
#                       Df Sum Sq Mean Sq F value     Pr(>F)    
#shann$site_FH_zone     22  19.62  0.8920   5.773 6.45e-12 ***
#Residuals             169  26.11  0.1545                     
#That tells us that there is a difference, but does not tell us which means are different. A Tukey's Honest Significant Difference (HSD) test can do pairwise comparisons of the means to find this out. We will use the HSD.test function from the agricolae package since it provides grouping codes that are useful for graphing.https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
#rem:try this for the letters!

#other way: to try do a post hoc
TUKEY <- TukeyHSD(aov, ordered = FALSE, conf.level = 0.95) # add "shann$site_FH_zone" ?
print(TUKEY)
generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}
# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY, "shann$site_FH_zone")
# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)
# Draw the basic boxplot
a <- boxplot(shann$Shannon~shann$site_FH_zone , ylim=c(min(data$value) , 1.1*max(data$value)) , col=my_colors[as.numeric(LABELS[,1])] , ylab="value" , main="")
# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )
#Add the labels
text( c(1:nlevels(data$treatment)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=my_colors[as.numeric(LABELS[,1])] )



### linear modeling on alpha

alpha_div_model <- lm(Shannon ~ fire + herbiv, data = shann)
summary(alpha_div_model)
plot(alpha_div_model)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  4.37681    0.13945  31.387   <2e-16 ***
#fireN       -0.03652    0.17241  -0.212    0.832    
#herbivN      0.14346    0.17001   0.844    0.400  
#Residual standard error: 0.4915 on 189 degrees of freedom
#Multiple R-squared:  0.001791,	Adjusted R-squared:  -0.008772 
#F-statistic: 0.1695 on 2 and 189 DF,  p-value: 0.8442

Pval <- anova(alpha_div_model) [] #creates a nice table

qqnorm(rstudent(alpha_div_model)); abline(0,1) 

#We can test whether the residuals are normally distributed (then we don't need the data to be normally distributed): Kolmogorov-Smirnov test. If p value >0.05 => we are OK

ks.test(rstudent(alpha_div_model), pnorm, mean = mean(rstudent(alpha_div_model)), sd=sd((rstudent(alpha_div_model))))
#data:  rstudent(alpha_div_model)
#D = 0.084009, p-value = 0.133
#alternative hypothesis: two-sided    

alpha_div_model2 <- lm(Shannon ~ fire * herbiv, data = temp)
summary(alpha_div_model2)
plot(alpha_div_model2)

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    5.39133    0.06515  82.755   <2e-16 ***
#fireN          0.02101    0.10301   0.204    0.839    
#herbivN        0.04210    0.09213   0.457    0.648    
#fireN:herbivN -0.12259    0.14458  -0.848    0.398    
#Residual standard error: 0.4919 on 188 degrees of freedom
#Multiple R-squared:  0.005593,	Adjusted R-squared:  -0.01027 
#F-statistic: 0.3525 on 3 and 188 DF,  p-value: 0.7874

#REM : U1, U22 and U2 are outliers (but within the cook's distance...)

Pval <- anova(alpha_div_model2) [] #creates a nice table

qqnorm(rstudent(alpha_div_model2)); abline(0,1) 

#We can test whether the residuals are normally distributed (then we don't need the data to be normally distributed): Kolmogorov-Smirnov test. If p value >0.05 => we are OK

ks.test(rstudent(alpha_div_model2), pnorm, mean = mean(rstudent(alpha_div_model2)), sd=sd((rstudent(alpha_div_model2))))
#data:  rstudent(alpha_div_model2)
#D = 0.093381, p-value = 0.07027


shannon_U <- read.csv("shannon_U.csv")
head(shannon_U)
alpha_div_U <- lm(Shannon ~ fire.freq * herbiv, data = shannon_U)
summary(alpha_div_U)
plot(alpha_div_U)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)         5.278607   0.227896  23.162   <2e-16 ***
#fire.freqC          0.076155   0.322293   0.236   0.8149    
#fire.freqT         -0.001613   0.338024  -0.005   0.9962    
#herbivN            -0.027993   0.322293  -0.087   0.9314    
#fire.freqC:herbivN -1.179046   0.455791  -2.587   0.0152 *  
#fire.freqT:herbivN -0.010332   0.478038  -0.022   0.9829    
#Residual standard error: 0.5582 on 28 degrees of freedom
#Multiple R-squared:  0.4231,	Adjusted R-squared:  0.3201 
#F-statistic: 4.107 on 5 and 28 DF,  p-value: 0.006367

alpha_div_U2 <- lm(Shannon ~ fire.freq + herbiv, data = shannon_U)
summary(alpha_div_U2)
plot(alpha_div_U2)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.488194   0.207009  26.512   <2e-16 ***
#fire.freqC  -0.513368   0.251689  -2.040   0.0503 .  
#fire.freqT  -0.006779   0.263974  -0.026   0.9797    
#herbivN     -0.447166   0.211461  -2.115   0.0429 *  
#Residual standard error: 0.6165 on 30 degrees of freedom
#Multiple R-squared:  0.2461,	Adjusted R-squared:  0.1707 
#F-statistic: 3.264 on 3 and 30 DF,  p-value: 0.03495

shannon_B <- read.csv("shannon_B.csv")
head(shannon_B)
alpha_div_B <- lm(Shannon ~ fire.freq * herbiv, data = shannon_B)
summary(alpha_div_B)
plot(alpha_div_B)
#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          5.5052     0.1897  29.017   <2e-16 ***
#fire.freqB          -0.1380     0.2683  -0.514    0.612    
#fire.freqC          -0.5206     0.3118  -1.670    0.110    
#herbivN             -0.1557     0.2814  -0.553    0.586    
#fire.freqB:herbivN   0.3768     0.3980   0.947    0.355    
#fire.freqC:herbivN       NA         NA      NA       NA    
#Residual standard error: 0.4647 on 21 degrees of freedom
#Multiple R-squared:  0.2504,	Adjusted R-squared:  0.1077 
#F-statistic: 1.754 on 4 and 21 DF,  p-value: 0.1759

alpha_div_B2 <- lm(Shannon ~ fire.freq + herbiv, data = shannon_B)
summary(alpha_div_B2)
plot(alpha_div_B2)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  5.41955    0.16638  32.573   <2e-16 ***
#fire.freqB   0.03324    0.19769   0.168   0.8680    
#fire.freqC  -0.62332    0.29156  -2.138   0.0439 *  
#herbivN      0.03268    0.19852   0.165   0.8708    
#Residual standard error: 0.4636 on 22 degrees of freedom
#Multiple R-squared:  0.2185,	Adjusted R-squared:  0.1119 
#F-statistic:  2.05 on 3 and 22 DF,  p-value: 0.1362
