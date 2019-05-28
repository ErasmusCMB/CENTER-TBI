library(caret)
library(rms)
library(data.table)
library(PCAmixdata)
library(bgravesteijn)
library(ggplot2)
library(MASS)
library(lubridate)
library(mice)
library(tableone)
library(miceadds)
library(knitr)
library(nnet)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(readr)

dt      <- data.frame(read_csv("Data/Data_6-12-2018.csv"))
imaging <- data.frame(read_csv("Data/Imaging_18-03-2019.csv"))

###########initial data analysis#################
factorvars <- colnames(dt)[c(5:8,19:23)]
for(i in factorvars){
  dt[,i] <- factor(dt[,i])
}

levels(dt$Subject.GOSE6monthEndpointDerived)        <- 2:8
levels(dt$InjuryHx.EDComplEventHypoxia)             <- c(0,1,1,NA)
levels(dt$InjuryHx.EDComplEventHypotension)         <- c(0,1,1,NA)
levels(dt$InjuryHx.EDICPMonitoring)                 <- c(0,1,NA) 
levels(dt$InjuryHx.EmergSurgInterventionsExtraCran) <- c(0,1,NA)
levels(dt$InjuryHx.EmergSurgInterventionsIntraCran) <- c(0,1,NA)
levels(dt$InjuryHx.InjCause)                        <- c("RTA", "Fall", "Other", 
                                                         "Violence/suicide", "Violence/suicide",
                                                         "Violence/suicide", "Other", NA)


dt$InjuryHx.BestOfExternaISS <- NULL
dt$mei                       <- apply(dt[,c("InjuryHx.PelvicGirdleAIS",
                                            "InjuryHx.ThoracicSpineAIS",
                                            "InjuryHx.AbdomenPelvicContentsAIS",
                                            "InjuryHx.UpperExtremitiesAIS",
                                            "InjuryHx.ExternaAIS",
                                            "InjuryHx.FaceAIS",
                                            "InjuryHx.ThoraxChestAIS",
                                            "InjuryHx.LumbarSpineAIS")],
                                      MARGIN = 1, FUN = max, na.rm=TRUE)
dt$mei                       <- as.numeric(dt$mei>3)
dt$mei                       <- factor(dt$mei)

imaging    <- data.table(imaging)
imaging.cn <- colnames(imaging)
imaging    <- imaging[,.(first(Imaging.SubduralHematomaSubacuteChronic),
                         first(Imaging.TraumaticSubarachnoidHemorrhage),
                         first(Imaging.EpiduralHematoma),
                         first(Imaging.SubduralHematomaAcute),
                         first(Imaging.SkullFracture),
                         first(Imaging.SubduralCollectionMixedDensity),
                         first(Imaging.CisternalCompression),
                         first(Imaging.MidlineShift),
                         first(Imaging.MassLesion),
                         first(Imaging.IntraventricularHemorrhage),
                         first(Imaging.TAI),
                         first(Imaging.Contusion)),
                      by=gupi]
colnames(imaging) <- imaging.cn
imaging <- data.frame(imaging)

for(i in colnames(imaging)[-1]){
  imaging[,i] <- factor(imaging[,i])
}

for(i in colnames(imaging)){
  if(length(levels(imaging[,i]))==3){
    levels(imaging[,i]) <- c(0,1,1)
  }
  if(length(levels(imaging[,i]))==2){
    levels(imaging[,i]) <- c(0,1)
  }
}

dt <- merge(dt, imaging, by="gupi")

##add los and til
til <- data.table(read.csv("Data/til_11-1-2019.csv"))
til <- til[,median(DailyTIL.TotalTIL, na.rm = TRUE), by=gupi]
til$V1 <- ifelse(is.na(til$V1),0,til$V1)
colnames(til) <- c("gupi", "mediantil")
dt <- merge(dt, til, by="gupi", all.x=TRUE)

los <- read.csv("Data/disch_dates_11-1-2019.csv")
str <- read.csv("Data/stratum_11-1-2019.csv")
los <- merge(los, str, by="gupi")
los$los <- difftime(time1 = ymd(los$Hospital.HospDischDate), 
                    time2 = ymd("1970-1-1"), units = "day")
los$los[los$los<0|los$los>100] <- NA
los$loicus <- difftime(time1 = ymd(los$Hospital.ICUDischDate), 
                       time2 = ymd("1970-1-1"), units = "day")
los$loicus[los$loicus<0|los$loicus>100] <- NA
los$loicus[los$Subject.PatientType!=3&is.na(los$loicus)] <- 0 
los <- los[,c("gupi", "los", "loicus")]
los$loicus <- as.numeric(los$loicus)
los$los    <- as.numeric(los$los)
dt <- merge(dt, los, by="gupi", all.x=TRUE)

dt <- dt[,c("gupi","Subject.Age",
            "InjuryHx.GCSMotorBaselineDerived","InjuryHx.GCSScoreBaselineDerived",
            "InjuryHx.PupilsBaselineDerived","InjuryHx.EDComplEventHypoxia", 
            "InjuryHx.EDComplEventHypotension","Subject.GOSE6monthEndpointDerived",
            "InjuryHx.EDICPMonitoring","InjuryHx.EmergSurgInterventionsIntraCran",
            "InjuryHx.EmergSurgInterventionsExtraCran","Subject.Sex",
            "InjuryHx.InjCause", "mei",
            "Imaging.SubduralHematomaSubacuteChronic","Imaging.TraumaticSubarachnoidHemorrhage",
            "Imaging.EpiduralHematoma","Imaging.SubduralHematomaAcute",
            "Imaging.SkullFracture","Imaging.SubduralCollectionMixedDensity",  
            "Imaging.CisternalCompression","Imaging.MidlineShift",
            "Imaging.MassLesion","Imaging.IntraventricularHemorrhage", 
            "Imaging.TAI", "Imaging.Contusion",
            "los", "loicus", "mediantil")]

## Descriptives#####
vars    <- colnames(dt)[-1][c(1,11,12,2,3,4,5,6,13,14:25,8:10,7, 26:28)]
contvar <- colnames(dt)[c(2,3,4, 27:29)]
catvar  <- vars[-which(vars%in%contvar)]

tbl1   <- CreateTableOne(vars = vars, data = dt, factorVars = catvar)
tbl1.p <- print(tbl1, nonnorm=contvar, missing=TRUE, printToggle=FALSE)
kable(tbl1.p)
write.csv2(tbl1.p, file = "Output/tbl1.csv")

## impute data ###########
dti <- mice(dt, maxit=0)
meth <- dti$method
pred <- dti$predictorMatrix

pred[, c("gupi")]        <- 0
meth[which(meth=="pmm")] <- "midastouch"

dti <- mice(dt, method=meth, predictorMatrix = pred, m=1, seed = 1234, printFlag = FALSE)
plot(dti)
densityplot(dti)
save(dti, file = "Data/imputed.RData")

#########PCA on imaging data###########
load("Data/imputed.RData")
dti <- mice::complete(dti,1)
dti.pca <- dti[,which(substr(colnames(dt),start = 1, stop=7)=="Imaging")]

res.pcamix <- PCAmix(X.quali = dti.pca, 
                     rename.level=TRUE, graph= FALSE)

#cumulative explained variance
cum.var <- data.frame(dim=rep(0:length(res.pcamix$eig[,"Cumulative"]),length(dti.pca)), 
                      cum.var=c(0,res.pcamix$eig[,"Cumulative"]),
                      imp=rep(0:length(res.pcamix$eig[,"Cumulative"]),
                              each=length(dti.pca)))
cum.var <- data.table(cum.var)
cum.var <- cum.var[,.(mean(cum.var), min(cum.var), max(cum.var)), by=dim]


plot(x=cum.var$dim, y=cum.var$V2, 
     xlab="Principal Component", 
     ylab="Cumulative Proportion of Variance Explained", 
     type="l", 
     ylim=c(0,100),
     xlim=c(0,12))

#loadings per component
plot(res.pcamix, choice ="levels", xlim=c(-1.5,2.5))
plot(res.pcamix, choice = "sqload", coloring.var=TRUE, leg=TRUE, posleg="topright")

# add first four principal components to dataset
pred <- predict(res.pcamix, 
                X.quali=dti.pca,
                rename.level=TRUE)[,1:4]
dti <- cbind(dti,pred)

#####Clustering#### 
library(cluster)

#define right variables
cluster <- dti[,-c(1,3,8,9:12,15:24)]

set.seed(1234)
## calculate gower's distance
gower_dist <- daisy(cluster,
                    metric="gower",
                    type=list(logratio = 3))
## identify 4 clusters based on gower's distance
pam_fit <- pam(x = gower_dist,
               diss=TRUE,
               k=4)

# add clusters to dataframe
dt$cluster  <- pam_fit$clustering
dti$cluster <- pam_fit$clustering

###### describe clusters
tbl_clust <- CreateTableOne(vars = vars, 
                            data = data.frame(dt),
                            strata = "cluster", 
                            factorVars = catvar)
tbl_clust.p <- print(tbl_clust, 
                     nonnorm=contvar, 
                     missing=TRUE, 
                     printToggle=FALSE)
kable(tbl_clust.p)

ggplot(dt, aes(x=Subject.Age))+geom_histogram()+facet_wrap(~cluster)


############### most defining characteristics

#add clusters to data
cluster_fit <- cbind(cluster, fit=pam_fit$clustering)

set.seed(123)

# vector to store results in
res_p_r2 <- rep(0, 11)

# model with no variables
nullmod <- multinom(fit~1, data=cluster_fit)

# full model, with all variables
fullmod <- multinom(fit~., data=cluster_fit)

#calculate r2 of full model 
r2_calc <- function(data, model, nullmodel){
  CoxSnell <- 1 - exp(-(2/nrow(data) * (logLik(model) -  logLik(nullmodel))))
  r2_full <- CoxSnell/(1 - exp((2 * nrow(data)^(-1)) * logLik(nullmodel)))
  return(r2_full)
}

## full r2
r2_full <- r2_calc(data = cluster_fit,model = fullmod, nullmodel = nullmod)

for (i in 1:11){# for each variable
  fit <- multinom(fit ~ ., data=cluster_fit[,-i]) #fit the model without the variable
  #compare the new r2 with the full r2
  res_p_r2[i] <- r2_full - r2_calc(data = cluster_fit,
                                   model = fit,
                                   nullmodel = nullmod)
}


## plot the partial r2's
res_pr2 <-sprintf("%.2f",res_p_r2)
names(res_pr2) <- colnames(cluster_fit)[1:11]
plot.df <- data.frame(part_r2 = res_p_r2, var = names(res_pr2))
plot.df$var <- factor(plot.df$var, levels=plot.df[order(plot.df$part_r2),]$var)
levels(plot.df$var) <- c("Hypotension", "Pupils", "Dim4", "Hypoxia", "Dim3", "Dim2", "Age", "Dim1", "GCS", "Major extracranial injury", "Injury cause")
ggplot(plot.df, aes(x=var, y=part_r2))+geom_bar(stat="identity")+
  coord_flip()+labs(x="Variable", y="Partial Nagelkerke R2")+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=13))
p3


#############Stability######
n.try <- 4 #number of iterations
#to save the results (the cluster per iteration)
cluster.fit <- expand.grid(id=1:nrow(dti),try=1:n.try, fit=NA) 

# cluster n.try times
set.seed(1234)
for(j in 1:n.try){
  index <- sample(x = 1:nrow(cluster), size = nrow(cluster), replace = TRUE)
  gower_dist <- daisy(cluster[index,],
                      metric="gower",
                      type=list(logratio = 3))
  pam_fit <- pam(x = gower_dist,
                 diss=TRUE,
                 k=4)
  cluster.fit[cluster.fit$try==j,][index,]$fit <- pam_fit$clustering
}

cluster.fit$fit <- factor(cluster.fit$fit)


### plot the flow over the clusters
alluv4   <- cluster.fit
alluv4_gr = alluv4  %>%
  spread( key = try, value = fit ) %>%
  dplyr::select( - id ) %>%
  group_by_all() %>%
  count() %>%
  ungroup() %>%
  mutate( alluvium = row_number() ) %>%
  rename( weight = n ) %>%
  gather( key = 'try', value = 'stratum', - weight, -alluvium ) %>%
  mutate( try = forcats::as_factor(try) )
ggplot(alluv4_gr[!is.na(alluv4_gr$stratum),],
             aes(x = try, stratum = stratum, 
                 alluvium = alluvium, y = weight, 
                 fill = stratum, label = stratum)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom")+
  theme_bw()




##measure of stability:
##of those from the same cluster, 
##who stay together in the same cluster, 
##how much stay together in one group 
##in the next try
n.try <- 1000
n.cluster <- 4

#to save the results in
cluster_id <- vector("list", length = n.try) # the cluster ids for each patient
p.stability <- rep(0, n.try)                 # the stability metric per iteration

set.seed(1234)
for(j in 1:n.try){
  # cluster with a random sample (with replacement)
  index <- sample(x = 1:nrow(cluster), size = nrow(cluster), replace = TRUE)
  gower_dist <- daisy(cluster[index,],
                      metric="gower",
                      type=list(logratio = 3))
  pam_fit <- pam(x = gower_dist,
                 diss=TRUE,
                 k=n.cluster)
  
  #save id and cluster
  cluster_id[[j]] <- data.frame(id=index, cluster=pam_fit$clustering)
  
  if(j>1){ # only if you have more clustered more than once (otherwise you cannot compare)
    first_try <- cluster_id[[(j-1)]]   # the cluster ids and clusters in the previous iteration
    second_try <- cluster_id[[j]]      # the cluster ids and clusters in this iteration
    which_id_present <- second_try$id  # the ids in the this iteration
    first_try_also_second_try <- first_try[first_try$id%in%second_try$id,] # which ids and their corresponding clusters are in both iteration
    stability_p_in_loop <- rep(0, n.cluster) # to save the results in
    for(k in 1:n.cluster){ # for each cluster
      first_grp1 <- first_try_also_second_try[first_try_also_second_try$cluster==k,
                                              "id"] # the id's of cluster k in the previous iteration
      # the proportion of patients in the current iteration within the same cluster as those of in the previous iteration
      stability_p_in_loop[k] <- max(prop.table(table(second_try[second_try$id%in%first_grp1,
                                                                "cluster"])))0
    }
    p.stability[j] <- mean(stability_p_in_loop) #the average over the four clusters
  }
  
}

p.stability <- p.stability[-1] # because the first iteration cannot be compared tot the 0st iteration

quantile(p.stability, probs=c(0.5,0.025,0.975))


#########compare prognosis of the groups
fit_mort <- glm(Subject.GOSE6monthEndpointDerived==2~
                  Subject.Age+
                  InjuryHx.GCSMotorBaselineDerived+
                  InjuryHx.PupilsBaselineDerived+
                  InjuryHx.EDComplEventHypotension+
                  InjuryHx.EDComplEventHypoxia+
                  Subject.Sex+
                  Imaging.TraumaticSubarachnoidHemorrhage+
                  Imaging.EpiduralHematoma,
                family="binomial", data=dti)
fit_GOSE <- polr(Subject.GOSE6monthEndpointDerived~
                   Subject.Age+
                   InjuryHx.GCSMotorBaselineDerived+
                   InjuryHx.PupilsBaselineDerived+
                   InjuryHx.EDComplEventHypotension+
                   InjuryHx.EDComplEventHypoxia+
                   Subject.Sex+
                   Imaging.TraumaticSubarachnoidHemorrhage+
                   Imaging.EpiduralHematoma,
                 data=dti)
dti$pred_mort <- predict(fit_mort, type = "response")

gose_df <- data.frame(prop.table(table(dti$cluster, dti$Subject.GOSE6monthEndpointDerived), 1))
colnames(gose_df) <- c("cluster", "GOSE", "Percentage")


gose_plot <- ggplot(gose_df, aes(x=cluster, fill=GOSE, y=Percentage))+
  geom_bar(stat = "identity")+scale_fill_brewer(palette="RdYlGn")+
  theme(axis.title.x=element_text(size=24),
        axis.text.x=element_text(size=22),
        axis.title.y=element_text(size=24),
        axis.text.y=element_text(size=22),
        legend.title=element_text(size=24),
        legend.text=element_text(size=22),
        legend.key.height = unit(1, units = "cm"))
gose_plot



ggplot(dti, aes(x=pred_mort))+geom_histogram()+facet_wrap(~cluster)

##calibration plots
source("C://Users/276018/Dropbox/onderzoek/How to/Val.prob.ci.2/auc.nonpara.mw.R")
source("C://Users/276018/Dropbox/onderzoek/How to/Val.prob.ci.2/ci.auc.R")
source("C://Users/276018/Dropbox/onderzoek/How to/Val.prob.ci.2/val.prob.ci.2.R")
par(mfrow=c(2,2))
for(i in 1:4){
  val.prob.ci.2(p = dti$pred_mort[dti$cluster==i], 
                y = dti$Subject.GOSE6monthEndpointDerived[dti$cluster==i]==2,
                logit = "p", main=paste("Cluster", i))
}

lrm(Subject.GOSE6monthEndpointDerived~InjuryHx.GCSScoreBaselineDerived, data=dti)
lrm(Subject.GOSE6monthEndpointDerived~mei+
      InjuryHx.GCSScoreBaselineDerived+
      InjuryHx.InjCause, data=dti)
lrm(Subject.GOSE6monthEndpointDerived~mei+
      InjuryHx.GCSScoreBaselineDerived+
      InjuryHx.InjCause+
      dim1+
      Subject.Age+
      dim2, data=dti)
lrm(Subject.GOSE6monthEndpointDerived~InjuryHx.PupilsBaselineDerived+
      InjuryHx.GCSScoreBaselineDerived+
      Subject.Age, data=dti)
lrm(Subject.GOSE6monthEndpointDerived~InjuryHx.PupilsBaselineDerived+
      InjuryHx.GCSScoreBaselineDerived+
      Subject.Age+
      dim1, data=dti)



##In  the nine proposed classes
class <- rep(0, nrow(df))
class[df$RTA==1&
df$mei==1&
dt$InjuryHx.GCSScoreBaselineDerived>12]<-1
class[df$RTA==1&
df$mei==1&
(dt$InjuryHx.GCSScoreBaselineDerived<=12&
dt$InjuryHx.GCSScoreBaselineDerived>8)]<-2
class[df$RTA==1&
df$mei==1&
dt$InjuryHx.GCSScoreBaselineDerived<=8]<-3
class[df$RTA==1&
df$mei==0&
dt$InjuryHx.GCSScoreBaselineDerived>12]<-4
class[df$RTA==1&
df$mei==0&
(dt$InjuryHx.GCSScoreBaselineDerived<=12&
dt$InjuryHx.GCSScoreBaselineDerived>8)]<-5
class[df$RTA==1&
df$mei==0&
dt$InjuryHx.GCSScoreBaselineDerived<=8]<-6
class[df$RTA==0&
dt$InjuryHx.GCSScoreBaselineDerived>12]<-7
class[df$RTA==0&
(dt$InjuryHx.GCSScoreBaselineDerived<=12&
dt$InjuryHx.GCSScoreBaselineDerived>8)]<-8
class[df$RTA==0&
dt$InjuryHx.GCSScoreBaselineDerived<=8]<-9

class[class==0] <- NA
dt$class <- factor(class)

df.class_gose <- data.frame(prop.table(table(dt$class, dt$Subject.GOSE6monthEndpointDerived), 1))
colnames(df.class_gose) <- c("class", "GOSE", "freq")
df.class_gose$class<-factor(df.class_gose$class, levels=9:1)
gose_classes <-ggplot(df.class_gose, aes(x=class, y=freq*100, fill=GOSE))+
  geom_bar(stat="identity")+
  scale_x_discrete("N", labels=table(df$class)[9:1])+
  coord_flip()+
  ylab("%")+
  scale_fill_brewer(palette="RdYlGn")+
  theme(axis.title.x=element_text(size=24),
        axis.text.x=element_text(size=22),
        axis.title.y=element_text(size=24),
        axis.text.y=element_text(size=22),
        legend.title=element_text(size=24),
        legend.text=element_text(size=22),
        legend.key.height = unit(1, units = "cm"))

gose_classes


####ct scans by proposoed class
dt$sah <- ifelse(dt$Imaging.TraumaticSubarachnoidHemorrhage==1, 1, 0)
dt$dai <- ifelse(dt$Imaging.TAI==1, 1, 0)
dt$cont <- ifelse(dt$Imaging.Contusion==1, 1, 0)
dt$hem <- ifelse(dt$Imaging.EpiduralHematoma==1|
dt$Imaging.EpiduralHematoma==1|
dt$Imaging.SubduralHematomaAcute==1|
dt$Imaging.SubduralCollectionMixedDensity==1|
dt$Imaging.IntraventricularHemorrhage==1, 1, 0)
df.ct_class <- data.table(dt)[,.(prop.table(table(sah))[2],
prop.table(table(dai))[2],
prop.table(table(cont))[2],
prop.table(table(hem))[2]),
by=class]
colnames(df.ct_class) <- c("class", "sah", "dai", "cont", "hen")
df.ct_class.l <- expand.grid(class=df.ct_class$class,
ct=c("sah", "dai", "cont", "hen"),
p=0)
for(i in c("sah", "dai", "cont", "hen")){
df.ct_class.l[df.ct_class.l$ct==i,]$p<-data.frame(df.ct_class)[,i]
}
df.ct_class.l <- df.ct_class.l[!is.na(df.ct_class$class),]

####se
nctcls <- data.table(dt)[,.(table(sah)[2],
table(dai)[2],
table(cont)[2],
table(hem)[2]),
by=class]
colnames(nctcls)<-c("class",  levels(df.ct_class.l$ct))
df.ct_class.l$n <- 0
df.ct_class.l$class <-factor(df.ct_class.l$class)
for(i in 1:length(levels(df.ct_class.l$class))){
df.ct_class.l$n[df.ct_class.l$class==i] <- unlist(c(nctcls[nctcls$class==i,]))[-1]
}
df.ct_class.l$se <- sqrt((df.ct_class.l$p*(1-df.ct_class.l$p))/df.ct_class.l$n)

df.ct_class.l$lo <- df.ct_class.l$p-1.96*df.ct_class.l$se
df.ct_class.l$hi <- df.ct_class.l$p+1.96*df.ct_class.l$se

new.lvl <- levels(df.ct_class.l$ct)[order(aggregate(p~ct,data =  df.ct_class.l, FUN=mean)$p, decreasing = TRUE)]
df.ct_class.l$ct <- factor(df.ct_class.l$ct, levels=new.lvl)
df.ct_class.l[df.ct_class.l<0]<-0
plot.ct.class <-ggplot(df.ct_class.l[!is.na(df.ct_class.l$class),], 
aes(x=ct, y=p*100, ymin=lo*100, ymax=hi*100))+
geom_linerange(cex=2)+
geom_point(cex=5)+
scale_x_discrete("CT finding", labels=c("SAH", "DAI", "Contustion", "Hematoma"))+
scale_y_continuous(limits=c(0,100))+
facet_wrap(~class)+
ylab("%")+
theme(axis.title.x=element_text(size=24),
axis.text.x=element_text(size=22, angle=90),
axis.title.y=element_text(size=24),
axis.text.y=element_text(size=22),
legend.title=element_text(size=24))
plot.ct.class
pdf("Output/plot.ct.class.pdf")
plot.ct.class
dev.off()

####by fitted class

df.ct_class <- data.table(dt)[,.(prop.table(table(Imaging.SubduralHematomaSubacuteChronic))[2],
prop.table(table(Imaging.TraumaticSubarachnoidHemorrhage))[2],
prop.table(table(Imaging.EpiduralHematoma))[2],
prop.table(table(Imaging.SubduralCollectionMixedDensity))[2],
prop.table(table(Imaging.CisternalCompression))[2],
prop.table(table(Imaging.MidlineShift))[2],
prop.table(table(Imaging.MassLesion))[2],
prop.table(table(Imaging.IntraventricularHemorrhage))[2]),
by=cluster]
colnames(df.ct_class) <- c("cluster", "SDH.c", "t.SAH", "EDH", "SDCMD", "Cistern", 
"Shift", "Mass", "IV.H")
df.ct_class.l <- expand.grid(cluster=df.ct_class$cluster,
ct=c("SDH.c", "t.SAH", "EDH", "SDCMD", "Cistern", 
"Shift", "Mass", "IV.H"),
p=0)
for(i in c("SDH.c", "t.SAH", "EDH", "SDCMD", "Cistern", 
"Shift", "Mass", "IV.H")){
df.ct_class.l[df.ct_class.l$ct==i,]$p<-data.frame(df.ct_class)[,i]
}



##standard error
nctcls <- data.table(dt)[,.(table(sah)[2],
table(dai)[2],
table(cont)[2],
table(hem)[2]),
by=cluster]
colnames(nctcls)<-c("cluster",  levels(df.ct_class.l$ct))
df.ct_class.l$n <- 0
for(i in 1:length(unique(df.ct_class.l$cluster))){
df.ct_class.l$n[df.ct_class.l$cluster==i] <- unlist(c(nctcls[nctcls$cluster==i,]))[-1]
}
df.ct_class.l$se <- sqrt((df.ct_class.l$p*(1-df.ct_class.l$p))/df.ct_class.l$n)

df.ct_class.l$lo <- df.ct_class.l$p-1.96*df.ct_class.l$se
df.ct_class.l$hi <- df.ct_class.l$p+1.96*df.ct_class.l$se

new.lvl <- levels(df.ct_class.l$ct)[order(aggregate(p~ct,data =  df.ct_class.l, FUN=mean)$p, decreasing = TRUE)]
df.ct_class.l$ct <- factor(df.ct_class.l$ct, levels=new.lvl)

df.ct_class.l[df.ct_class.l<0]<-0
plot.ct.class <-ggplot(df.ct_class.l[df.ct_class.l$ct!="SDH.c"&
!is.na(df.ct_class.l$cluster),], 
aes(x=ct, y=p*100, ymin=lo*100, ymax=hi*100))+
geom_linerange(cex=2)+
geom_point(cex=5)+
scale_x_discrete("CT finding", labels=c("t-SAH", "Mass", "Cistern", "IV-H",
"Shift","EDH", "SDCMD"))+
facet_wrap(~cluster)+
scale_y_continuous(limits=c(0,100))+
ylab("%")+
theme(axis.title.x=element_text(size=24),
axis.text.x=element_text(size=22, angle=90),
axis.title.y=element_text(size=24),
axis.text.y=element_text(size=22),
legend.title=element_text(size=24))
plot.ct.class
pdf("Output/plot.ct.cluster.pdf")
plot.ct.class
dev.off()


fit_sah<- glm(sah~class, family=binomial,data=dt)
summary(fit_sah)
fit_sah<- glm(dai~class, family=binomial,data=dt)
summary(fit_sah)
fit_sah<- glm(cont~class, family=binomial,data=dt)
summary(fit_sah)
fit_sah<- glm(hem~class, family=binomial,data=dt)
summary(fit_sah)
