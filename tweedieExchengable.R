library(dplyr)
library(haven)
library(stats)
library(statmod)
library(tweedie)
library(labelled)
library(foreign)
library(purrr)
library(stringr)
## library(Rmisc) for CI calculation of the mean, unload to prevent conflict of dplyr and plyr
## library(ggplot2) For plots, unload to prevent conflict of dplyr and plyr
df1 <- read.csv("C:/redd/paper3/tweedie_kheus_data/Rcodes2013/kheus2013_data_questionnaire/data/gok1.csv")
head(df1)
View(df1)
attach(df1)
#dog <- "The quick brown dog"    
#str_to_title(dog)
#df1$county = str_to_title(County)
#################---------------------------------------------------for bubble summary data
########-----example starts here
id<-c(1,1,1,1,2,2,2,2)
x<-c(1,0,0,1,0,1,1,0)
djk<-data.frame(id,x)
head(djk)

a1<-djk %>%
  group_by(id) %>%
  summarise(sum_of_1 = sum(x!=0),
            pct= round((sum_of_1/n())*100))
a1
##########----example ends here
#percentages of admissions with reasons
###---------------------------------------------------------Code line courtesy of NOAH mutahi .PhD starts here
a2<-df1 %>%
  group_by(county) %>%
  summarise(meanDist=mean(Distance),
            sum_of_1 = sum(Malaria_fever!=0),pct1= round((sum_of_1/n())*100),
            sum_of_2 = sum(reason_staffQualified!=0),pct2= round((sum_of_2/n())*100),
            sum_of_3 = sum(reason_Less_costly!=0),pct3= round((sum_of_2/n())*100))
head(a2)
write.csv(a2, file="C:/redd/paper3/maps_paper_write/outPlots/staffQLessCost.txt")
###---------------------------------------------------------Code line courtesy of NOAH mutahi .PhD Ends here

##----------------------------------------------------------mean distance and 95%ci
##for counties    load the Rmisc package for CI calculation
a3<-df1 %>%
  group_by(county) %>%
  dplyr::summarise(mean_dist = mean(Distance), 
                   Lower_CI = CI(Distance)[3], 
                   Upper_CI = CI(Distance)[1],
                   min_dist=min(Distance),
                   max_dist=max(Distance))
a3
write.csv(a3, file="C:/redd/paper3/maps_paper_write/outPlots/meanCI.txt")
#For kenya
a4<-df1 %>%
  group_by(fac_owned_govt) %>%
  dplyr::summarise(mean_Dist = mean(newDistBelowOneKm), 
                   Upper_CI = CI(newDistBelowOneKm)[1], 
                   Lower_CI = CI(newDistBelowOneKm)[3])
a4

power=tweedie.profile(Distance~w_index+age+residence+gender+Population.1000.+
                        PovertyIndex+Malaria_fever+reason_staffQualified+reason_Less_costly+reason_medicineavailable,
                      p.vec=seq(1.3,1.6,length=10),
                      do.plot=TRUE,do.ci=TRUE, method="interpolation", data = df1)
p=power$p.max
glmmodel<-glm(Distance~w_index+age+residence+gender+Population.1000.+
                PovertyIndex+Malaria_fever+reason_staffQualified+reason_Less_costly+reason_medicineavailable,
              family=tweedie(var.power=p,link.power=0),x=TRUE, data=df1)
fits<-glmmodel$fitted.values
beta=glmmodel$coefficients
phi = power$phi.max
n=length(Distance)
r=glmmodel$rank
dev=sum(tweedie.dev(Distance,fits,p))
devold=100*dev
epsilon = 1e-8

while (abs(dev - devold)/(0.1 + abs(dev)) > epsilon) {
  p.residuals=(Distance-fits)/sqrt(fits^p)
  p.residuals
  new.phi<-sum(p.residuals^2)/(n-r)
  new.phi
  
  df1$resd2<-p.residuals
  
  res_total_n_for_alph<-df1 %>%
    group_by(county) %>%
    summarise(resd = sum(combn(resd2, 2, FUN = prod)),n = n())%>%  mutate(sampl_n = n*(n-1)/2)
  res_total_n_for_alph 
  
  
  sum_res_alpha<-sum(res_total_n_for_alph$resd)
  sum_res_alpha
  sum_n_alpha<-sum(res_total_n_for_alph$sampl_n)
  sum_n_alpha
  alpha<-(1/new.phi)*(1/(sum_n_alpha-r))*sum_res_alpha
  alpha
  
  i=0
  R <- matrix(NA, n,n)
  for (i in 1:n){
    for (j in 1:n){
      if (i==j){
        R[i,j] <- 1.
      } else {
        R[i,j] <- alpha
      }
    }
  }
  A<-diag(fits)^(p/2)
  V = (A %% R  %% A)
  xmat=glmmodel$x
  D = matrix(nrow=n,ncol=r)
  for(i in (1:r)){
    D[,i]=fits*xmat[,i]
  }
  C=t(D) %% solve(V) %% D
  B=t(D) %% solve(V) %% (Distance-fits)
  beta1=beta + (solve(C) %*% B)
  fits<-exp(t(beta) %*% t(xmat))
  fits<-as.vector(fits)
  devold<-dev
  dev<-sum(tweedie.dev(Distance,fits,p))
}

quasi<-sum((Distance*fits^(1-p)/(1-p))-((fits^(2-p))/(2-p)))
qicu<-(-2*quasi)+(2*r)

marginal=(1/n)*sum(Distance) 
top=sum((Distance-fits)^2)
bottom=sum((Distance-marginal)^2) 
R2 = 1 - (top / bottom)
qicu
R2
alpha