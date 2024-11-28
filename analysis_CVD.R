
# Load data and packages -----------------------------------------------------

#Load packages
if (!require(pacman)) install.packages('pacman')
library(pacman)
p_load(tidyverse,survival,ggplot2,flexsurv,survminer,mice,meta,
       pec,riskRegression,gridExtra)


longDF <- complete(imp_tot, "long", include=TRUE) #Get long format of imputed data

n_imp <- 10 #Number of imputations
n_pt <- length(longDF[longDF$.imp==1,]$ID) #Number of individual patients

#Get data in list format
listDF <- lapply(1:n_imp,function(x){
  dat <- longDF %>% filter(.imp==x)
})

#Weight data to account for competing risks
full_weight <- lapply(1:n_imp,function(x){
  crprep(Tstop=listDF[[x]]$TTE_CVDDeath, status=listDF[[x]]$status, 
         id=listDF[[x]]$.id,
         keep=names(listDF[[x]]),
         data=listDF[[x]]) %>% select(-c(4,7))
})

rm(mi.res)

Studies <- c("HRS","CoLaus","HealthABC","SHARE")


# Fit cause-specific CIF model ----------------------------------------

fit.model <- list()

for(i in 1:n_imp){ #Fit model in each imputation dataset
  
  fit.model[[i]] <- flexsurvspline(Surv(Tstart,Tstop,status==1)~
                                     cAge+Gender+cBMI+CurDrink+Smoker+
                                     TreatHyperten_adj+Choltreat+
                                     DiabetesDurCat+AntiDiab_adj+
                                     cHbA1c+cTOTC+cHDLC,
                                   k=0, #Flexible parametric model with 0 knots is equivalent to a Weibull model
                                   scale="hazard",
                                   data=full_weight[[i]],
                                   weights=weight.cens)
}

# Pool coefficients -------------------------------------------------------

##Get coefficients
sum_imp <- lapply(fit.model,function(x){x$res})

Q <- sapply(sum_imp,function(x){x[,1]})
u <- sapply(sum_imp,function(x){x[,4]})

pool_imp <- lapply(1:nrow(Q),function(x){
  pool.scalar(Q=Q[x,],U=u[x,],n=n_pt,k=1)
})

coef_imp <- t(sapply(pool_imp,function(x){c(x$qbar,x$ubar)}))
colnames(coef_imp) <- c("estimate","se")
rownames(coef_imp) <- rownames(Q)

coef_imp <- as.data.frame(coef_imp)

coef_imp$lower <- coef_imp$estimate -1.96*coef_imp$se
coef_imp$upper <- coef_imp$estimate +1.96*coef_imp$se


# Assess apparent performance ---------------------------------------------

#Function that calculates performance measures
calculate_performance <- function(model = NULL, time = NULL, 
                                  status = NULL, data.used = NULL,
                                  time.horizon = NULL, primary.event = NULL,
                                  incl.calibration=TRUE,
                                  incl.oe=TRUE,
                                  incl.Cindex.W=TRUE,
                                  predicted.Score = NULL){
  
  pred <- predictRisk(
    object = model,
    cause = primary.event,
    newdata = data.used,
    times = time.horizon
  )
  
  if(incl.Cindex.W==TRUE){ #C-index (Wolber's)
    cindx.W <- unlist(cindex(object=list(pred),
                             data=data.used,
                             formula=Surv(TTE_CVDDeath,status)~1,
                             eval.times=time.horizon)$AppCindex)
    
    cindx.W <- as.numeric(cindx.W[1])
    
    ##Use Hanley's variance of the c-index as described in  doi: https://doi.org/10.1136/bmj.i6460
    N <- length(data.used$ID)
    s <- sum(data.used$CVD_event)
    t <- N-s
    varlogitc <- (1+(N/2-1)*(1-cindx.W)/(2-cindx.W)+(N/2-1)*cindx.W/(1+cindx.W)) / (cindx.W*(1-cindx.W)*s*t)
    cindx.W.logitse <- as.numeric(sqrt(varlogitc))
    
  } else {
    cindx.W <- NULL
    cindx.W.logitse <- NULL
  }
  
  if(!is.null(predicted.Score)){
    pred <- predicted.Score
  }
  
  if(incl.calibration==TRUE){
    
    score_vdata <- Score(
      list("sdh_validation" = pred),
      formula = Hist(TTE_CVDDeath, status) ~ 1, 
      cens.model = "km", 
      data = data.used, 
      conf.int = TRUE, 
      times = time.horizon,
      metrics = c("auc", "brier"),
      summary = c("ipa"), 
      cause = primary.event,
      plots = "calibration" )
    
    #Calibration
    pseudos <- data.frame(score_vdata$Calibration$plotframe)
    pseudos <- pseudos[order(pseudos$risk), ]
    pseudos$cll_pred <- log(-log(1 - pseudos$risk))  
    
    # Fit model for calibration intercept
    fit_cal_int <- geese(
      pseudovalue ~ offset(cll_pred), 
      data = pseudos,
      id = ID, 
      scale.fix = TRUE, 
      family = gaussian,
      mean.link = "cloglog",
      corstr = "independence", 
      jack = TRUE  )
    
    cal.int <- summary(fit_cal_int)$mean$estimate
    cal.int.se <- summary(fit_cal_int)$mean$san.se
    
    # Fit model for calibration slope
    fit_cal_slope <- geese(
      pseudovalue ~ offset(cll_pred) + cll_pred,
      data = pseudos,
      id = ID, 
      scale.fix = TRUE, 
      family = gaussian,
      mean.link = "cloglog",
      corstr = "independence", 
      jack = TRUE  )
    
    cal.slope <- 1 + summary(fit_cal_slope)$mean["cll_pred",]$estimate
    cal.slope.se <- summary(fit_cal_slope)$mean["cll_pred",]$san.se
    
    #Calibration plot
    smooth_pseudos <- predict(
      stats::loess(pseudovalue ~ risk, data = pseudos, degree = 1, span = 0.66), 
      se = TRUE)
    
  } else {
    pseudos <- NULL
    smooth_pseudos <- NULL
    cal.int <- NULL
    cal.int.se <- NULL
    cal.slope <- NULL
    cal.slope.se <- NULL
  }
  
  if(incl.oe==TRUE){
    # Calculate O/E
    obj <- summary(survfit(Surv(TTE_CVDDeath, event) ~ 1, 
                           data = data.used), 
                   times = time.horizon)
    aj <- list("obs" = obj$pstate[, primary.event + 1], "se" = obj$std.err[, primary.event + 1])
    OE <- aj$obs / mean(pred)
    OE.se <- aj$se / aj$obs
  } else {
    OE <- NULL
    OE.se <- NULL
  }
  
  returnVec <- list()
  returnVec[[1]] <- c("OE" = OE,
                      "OE.se" = OE.se,
                      "cal.int" = cal.int,
                      "cal.int.se"= cal.int.se,
                      "cal.slope" = cal.slope,
                      "cal.slope.se"= cal.slope.se,
                      "cindex.W" = cindx.W,
                      "cindex.W.logitse" = cindx.W.logitse
  )
  returnVec[[2]] <- pseudos
  returnVec[[3]] <- smooth_pseudos
  
  
  return(returnVec)
}


#Calculate 5-year performance

time.horizon <- 5*12
apparent.competing.list <- list()

for(i in 1:n_imp){
  apparent.competing.list[[i]]<- calculate_performance(model = fit.model[[i]],
                                                       time = listDF[[i]]$TTE_CVDDeath,
                                                       status = listDF[[i]]$status,
                                                       data.used = listDF[[i]],
                                                       time.horizon = time.horizon,
                                                       primary.event = 1)
}

apparent.indicators <- sapply(apparent.competing.list,`[[`,1)

source("functions_pool.R") #Load helper functions to pool performance estimates

app.oe <- pool_oe(est=apparent.indicators[1,],
                  se=apparent.indicators[2,],
                  n=n_pt)

app.cal.int <- pool_cal(est=apparent.indicators[3,],
                        se=apparent.indicators[4,],
                        n=n_pt)

app.cal.slope <- pool_cal(est=apparent.indicators[5,],
                          se=apparent.indicators[6,],
                          n=n_pt)

app.cindx.W <- pool_cindx(cindx=apparent.indicators["cindex.W",],
                          se=apparent.indicators["cindex.W.logitse",],
                          n=n_pt,
                          method="W")

apparent <- rbind("OE ratio" = app.oe,
                  "Calibration intercept" = app.cal.int,
                  "Calibration slope" = app.cal.slope,
                  "Wolber's C-index" = app.cindx.W)


# Calibration plot at 5-years

##5-year CVD risk prediction
pred.5y <- list()
for(i in 1:n_imp){
  pred.5y[[i]] <- predictRisk(
    object = fit.model[[i]],
    cause = 1,
    newdata = listDF[[i]],
    times = time.horizon
  )
}

## Combine values over imputation datasets
pseudos <- lapply(apparent.competing.list,`[[`,2)
smooth_pseudos <- lapply(apparent.competing.list,`[[`,3)

ps.risk <- Reduce(`+`, lapply(pseudos, `[[`, "risk"))/n_imp
ps.val <- Reduce(`+`, lapply(pseudos, `[[`, "pseudovalue"))/n_imp
ps.fit <- Reduce(`+`, lapply(smooth_pseudos, `[[`, "fit"))/n_imp
ps.df <- Reduce(`+`, lapply(smooth_pseudos, `[[`, "df"))/n_imp
ps.se <- Reduce(`+`, lapply(smooth_pseudos, `[[`, "se.fit"))/n_imp
predicted.competing.1 <- lapply(pred.5y,function(x)log(x/(1-x)))
predicted.competing.1 <- exp(Reduce(`+`, predicted.competing.1)/n_imp)/
  (1+exp(Reduce(`+`, predicted.competing.1)/n_imp))

centile_LP <- cut(predicted.competing.1 ,
                  breaks=quantile(predicted.competing.1 , prob = seq(0,1,0.1)),
                  labels=c(1:10),include.lowest=TRUE)
obj <- summary(survfit(Surv(TTE_CVDDeath, event) ~ centile_LP, 
                       data = listDF[[1]]), 
               times = 5*12)
aj <- data.frame(
  "centile_LP" = 1:10,
  "obs" = obj$pstate[, 2], 
  "se" = obj$std.err[, 2])

pts <- data.frame(pred=predicted.competing.1,centile_LP=as.integer(centile_LP))
pts <- pts %>%
  group_by(centile_LP) %>%
  summarise(mean.pred=mean(pred),
            sd.pred=sd(pred))
pts <- pts %>% left_join(aj)

spike_bounds <- c(-0.075, 0)
bin_breaks <- seq(0, 0.6, length.out = 100 + 1)
freqs <- table(cut(predicted.competing.1, breaks = bin_breaks))
bins <- bin_breaks[-1]
freqs_valid <- freqs[freqs > 0]
freqs_rescaled <- spike_bounds[1] + (spike_bounds[2] - spike_bounds[1]) * 
  (freqs_valid - min(freqs_valid)) / (max(freqs_valid) - min(freqs_valid))

## Build plot
plot(
  x = ps.risk, 
  y = ps.val,
  xlim = c(0, 1), 
  ylim = c(spike_bounds[1], 1),
  yaxt = "n",
  frame.plot = FALSE,
  xlab = "Estimated risks",
  ylab = "Observed outcome proportions", 
  type = "n"
)
axis(2, seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
polygon(
  x = c(ps.risk, rev(ps.risk)),
  y = c(
    pmax(ps.fit - qt(0.975, ps.df) * ps.se, 0),
    rev(ps.fit + qt(0.975, ps.df) * ps.se)
  ),
  border = FALSE,
  col = "lightgray"
)
abline(a = 0, b = 1, col = "gray")
lines(x = ps.risk, y = ps.fit, lwd = 2)
segments(
  x0 = bins[freqs > 0], 
  y0 = spike_bounds[1], 
  x1 = bins[freqs > 0], 
  y1 = freqs_rescaled,
       lwd=1.2
)
points(pts$mean.pred, pts$obs,pch=20)
arrows(x0=pts$mean.pred, y0=pts$obs-1.96*pts$se, x1=pts$mean.pred, y1=pts$obs+1.96*pts$se, code=3, angle=90, length=0.1, lwd=1)
arrows(x0=pts$mean.pred-pts$sd.pred, y0=pts$obs, x1=pts$mean.pred+pts$sd.pred, y1=pts$obs, code=3, angle=90, length=0.1, lwd=1)

                                
#Kaplan-Meier curves of risk groups

ci_fit <- list()
pred.df <- list()
pred.df_sum <- list()
for (i in 1:n_imp){
  listDF[[i]] <- cbind(listDF[[i]],pred.5y=pred.5y[[i]])
  listDF[[i]] <- listDF[[i]] %>%
    mutate(
      #5-year risk groups set to half of 10-year rsik groups as described in ESC diabetes guidelines
      # (doi:https://doi.org/10.1093/eurheartj/ehad192)
      Group = ifelse(pred.5y<0.025,"LoW",
                     ifelse(pred.5y<0.05,"Moderate",
                            ifelse(pred.5y<0.1,"High",
                                   ifelse(pred.5y>=0.1,"Very high",NA))))
      )
  
  ci_fit[[i]] <- cuminc(ftime = listDF[[i]]$TTE_CVDDeath/12, #Months to Years
                        fstatus=listDF[[i]]$status,
                        group=listDF[[i]]$Group,
                        cencode=0)
  ci_fit[[i]] <- ci_fit[[i]][-c(4:6)]
  
  pred.df[[i]] <- predictRisk(
    object = fit.model[[i]],
    cause = 1,
    newdata = listDF[[i]],
    times = seq(0,200)
  )
  
  pred.df[[i]] <- cbind(as.data.frame(pred.df[[i]]),Group=listDF[[i]]$Group)
  
  pred.df_sum[[i]] <- pred.df[[i]] %>% 
    group_by(Group) %>%
    summarise(across(everything(), mean))
  
}

comp.se <- ggcompetingrisks(ci_fit[[1]],xlab="Months",multiple_panels=F,conf.int=T)

se.vhigh <- filter(comp.se[["data"]],name=="Very high 1")$std
se.high <- filter(comp.se[["data"]],name=="High 1")$std
se.mod <- filter(comp.se[["data"]],name=="Moderate 1")$std
#There was no patient with low 5-year risk of CVD

ci_fit[[1]][["High 1"]]$ymin <- ci_fit[[1]][["High 1"]]$est + 1.96*se.high
ci_fit[[1]][["High 1"]]$ymax <- ci_fit[[1]][["High 1"]]$est - 1.96*se.high

ci_fit[[1]][["Moderate 1"]]$ymin <- ci_fit[[1]][["Moderate 1"]]$est + 1.96*se.mod
ci_fit[[1]][["Moderate 1"]]$ymax <- ci_fit[[1]][["Moderate 1"]]$est - 1.96*se.mod

ci_fit[[1]][["Low 1"]]$ymin <- ci_fit[[1]][["Low 1"]]$est + 1.96*se.low
ci_fit[[1]][["Low 1"]]$ymax <- ci_fit[[1]][["Low 1"]]$est - 1.96*se.low

## Build plot
p.km <- plot(ci_fit[[1]],
             xlab="Years",
             color=c("darkviolet","darkgoldenrod3","deepskyblue3"),
             lty=c(2,2),
             wh=c(0,3),
             xlim=c(0,10))

polygon(c(ci_fit[[1]][["Very high 1"]]$time, rev(ci_fit[[1]][["Very high 1"]]$time)), 
        c(ci_fit[[1]][["Very high 1"]]$ymin, rev(ci_fit[[1]][["Very high 1"]]$ymax)),
        col = adjustcolor("deepskyblue3", 0.1),
        border = 0)
polygon(c(ci_fit[[1]][["High 1"]]$time, rev(ci_fit[[1]][["High 1"]]$time)), 
        c(ci_fit[[1]][["High 1"]]$ymin, rev(ci_fit[[1]][["High 1"]]$ymax)),
        col = adjustcolor("darkviolet", 0.1),
        border = 0)
polygon(c(ci_fit[[1]][["Moderate 1"]]$time, rev(ci_fit[[1]][["Moderate 1"]]$time)), 
        c(ci_fit[[1]][["Moderate 1"]]$ymin, rev(ci_fit[[1]][["Moderate 1"]]$ymax)),
        col = adjustcolor("darkgoldenrod3", 0.1),
        border = 0)

lines(seq(0,200)/12, pred.df_sum[[1]][3,-1], col = "deepskyblue3", type = "l", lty = 1, lwd=2)
lines(seq(0,200)/12, pred.df_sum[[1]][1,-1], col = "darkviolet", type = "l", lty = 1, lwd=2)
lines(seq(0,200)/12, pred.df_sum[[1]][2,-1], col = "darkgoldenrod3", type = "l", lty = 1, lwd=2)

legend("topleft", 
       legend = c(
         "Observed",
         "Predicted"
       ), 
       col = c(rep("black",2)), 
       lty = c(2,1),
       bty = "n", 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.05))

legend("topleft", 
       legend = c(
         "Low (<2.5%)",
         "Moderate (2.5% to <5%)",
         "High (5% to <10%)",
         "Very high (>=10%)"
       ), 
       col = c("chartreuse4","darkgoldenrod3","darkviolet","deepskyblue3"), 
       lty = c(1,1,1,1),
       bty = "n", 
       cex = 1, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.3, 0.05))


# Apparent performance by study ------------------------------------------
N.study <- length(Studies)

apparent.study.list <- list("Indicators"=list(),
                            "Pseudo values"=list(),
                            "Smoothed pseudo values"=list())

#Loop through each study and calculate performance
for (i in 1:N.study){ 
  
  apparent.study.list.each <- list()
  
  for (j in 1:n_imp){
    studDF <- listDF[[j]] %>% filter(Study==Studies[i])
    n_pt_stud <- length(studDF$ID)
    
    apparent.study.list.each[[j]]<- calculate_performance(model = fit.model[[j]],
                                                          time = studDF$TTE_CVDDeath,
                                                          status = studDF$status,
                                                          data.used = studDF,
                                                          time.horizon = time.horizon,
                                                          primary.event = 1)
  }
  
  app.study.indicators <- sapply(apparent.study.list.each,`[[`,1)
  
  app.stud.oe <- pool_oe(est=app.study.indicators[1,],
                         se=app.study.indicators[2,],
                         n=n_pt_stud,
                         transform=F)
  
  app.stud.cal.int <- pool_cal(est=app.study.indicators[3,],
                               se=app.study.indicators[4,],
                               n=n_pt_stud)
  
  app.stud.cal.slope <- pool_cal(est=app.study.indicators[5,],
                                 se=app.study.indicators[6,],
                                 n=n_pt_stud)

  
  app.stud.cindx.W <- pool_cindx(cindx=app.study.indicators["cindex.W",],
                                 se=app.study.indicators["cindex.W.logitse",],
                                 n=n_pt,
                                 method="W",
                                 transform=F)
  
  
  apparent.study.list[[1]][[i]] <- c("OE ratio" = app.stud.oe,
                                     "Calibration intercept" = app.stud.cal.int,
                                     "Calibration slope" = app.stud.cal.slope,
                                     "Wolber's C-index" = app.stud.cindx.W)
  
  apparent.study.list[[2]][[i]] <- lapply(apparent.study.list.each,`[[`,2)
  apparent.study.list[[3]][[i]] <- lapply(apparent.study.list.each,`[[`,3)
  
}

app.study.indic <- as.matrix(bind_rows(apparent.study.list[[1]]))

#Meta-analyse results
oe.meta <- metagen(app.study.indic[,"OE ratio.Estimate"],
                   app.study.indic[,"OE ratio.logSE"],
                   studlab=Studies)

cal.int.meta <- metagen(app.study.indic[,"Calibration intercept.Estimate"],
                        app.study.indic[,"Calibration intercept.SE"],
                        studlab=Studies)

cal.slope.meta <- metagen(app.study.indic[,"Calibration slope.Estimate"],
                          app.study.indic[,"Calibration slope.SE"],
                          studlab=Studies)

cindx.W.meta <- metagen(app.study.indic[,"Wolber's C-index.Estimate"],
                        app.study.indic[,"Wolber's C-index.logitSE"],
                        studlab=Studies)

apparent.meta <- rbind(app.study.indic,
                       c(oe.meta$TE.random,oe.meta$seTE.random,oe.meta$lower.random,oe.meta$upper.random,
                         cal.int.meta$TE.random,cal.int.meta$seTE.random,cal.int.meta$lower.random,cal.int.meta$upper.random,
                         cal.slope.meta$TE.random,cal.slope.meta$seTE.random,cal.slope.meta$lower.random,cal.slope.meta$upper.random,
                         cindx.W.meta$TE.random,cindx.W.meta$seTE.random,cindx.W.meta$lower.random,cindx.W.meta$upper.random
                       )
)
apparent.meta[,1] <- exp(apparent.meta[,1])
apparent.meta[,3] <- exp(apparent.meta[,3])
apparent.meta[,4] <- exp(apparent.meta[,4])
apparent.meta[,13] <- exp(apparent.meta[,13])/(1+exp(apparent.meta[,13]))
apparent.meta[,15] <- exp(apparent.meta[,15])/(1+exp(apparent.meta[,15]))
apparent.meta[,16] <- exp(apparent.meta[,16])/(1+exp(apparent.meta[,16]))

heterogeneity.meta <- data.frame(
  OE = c(oe.meta$I2,oe.meta$tau2,oe.meta$pval.Q),
  Cal.Intercept = c(cal.int.meta$I2,cal.int.meta$tau2,cal.int.meta$pval.Q),
  Cal.Slope = c(cal.slope.meta$I2,cal.slope.meta$tau2,cal.slope.meta$pval.Q),
  CIndex.W = c(cindx.W.meta$I2,cindx.W.meta$tau2,cindx.W.meta$pval.Q),
)
rownames(heterogeneity.meta) <- c("I2","Tau2","Q P-val")


# Apparent performance by subgroup ----------------------------------------

subgroups <- c(
  AgeL = expression(cAge>=5),
  AgeH = expression(cAge<5),
  Men = expression(Gender==0),
  Women = expression(Gender==1)
)

N.subgrp <- length(subgroups)

apparent.subgrp.list <- list("Indicators"=list(),
                             "Pseudo values"=list(),
                             "Smoothed pseudo values"=list())

#Loop throigh each subgroup and calculate performance
for (i in 1:N.subgrp){
  apparent.subgrp.list.each <- list()
  for (j in 1:n_imp){
    subgrpDF <- listDF[[j]] %>% filter(eval(subgroups[i]))
    n_pt_subgrp <- length(subgrpDF$ID)
    
    apparent.subgrp.list.each[[j]]<- calculate_performance(model = fit.model[[j]],
                                                           time = subgrpDF$TTE_CVDDeath,
                                                           status = subgrpDF$status,
                                                           data.used = subgrpDF,
                                                           time.horizon = time.horizon,
                                                           primary.event = 1)
    
    
  }
  
  app.subgrp.indicators <- sapply(apparent.subgrp.list.each,`[[`,1)
  
  app.subgrp.oe <- pool_oe(est=app.subgrp.indicators[1,],
                           se=app.subgrp.indicators[2,],
                           n=n_pt_subgrp)
  
  app.subgrp.cal.int <- pool_cal(est=app.subgrp.indicators[3,],
                                 se=app.subgrp.indicators[4,],
                                 n=n_pt_subgrp)
  
  app.subgrp.cal.slope <- pool_cal(est=app.subgrp.indicators[5,],
                                   se=app.subgrp.indicators[6,],
                                   n=n_pt_subgrp)
  
  app.subgrp.cindx <- pool_cindx(cindx=app.subgrp.indicators[7,],
                                 se=app.subgrp.indicators[8,],
                                 n=n_pt_subgrp)
  
  app.subgrp.cindx.W <- pool_cindx(cindx=app.subgrp.indicators["cindex.W",],
                                   se=app.subgrp.indicators["cindex.W.logitse",],
                                   n=n_pt,
                                   method="W")
  
  apparent.subgrp.list[[1]][[i]] <- c("OE ratio" = app.subgrp.oe,
                                      "Calibration intercept" = app.subgrp.cal.int,
                                      "Calibration slope" = app.subgrp.cal.slope,
                                      "Wolber's C-index" = app.subgrp.cindx.W)
  
  apparent.subgrp.list[[2]][[i]] <- lapply(apparent.subgrp.list.each,`[[`,2)
  apparent.subgrp.list[[3]][[i]] <- lapply(apparent.subgrp.list.each,`[[`,3)
}

res.Subgroup <- bind_rows(apparent.subgrp.list[[1]])
rownames(res.Subgroup) <- names(subgroups)


# Internal validation via bootstrapping -------------------

set.seed(202405)

n.bootstrap <- 500
boot.eachfold <- matrix(NA, n.bootstrap, 8*1)

start_time <- Sys.time()
for (i in 1:1){ #Only bootstrap in 1 imputation dataset due to model runtime
  
  for(j in 1:n.bootstrap){
    
    boot.sample <- sample(n_pt, replace = T) 
    
    # create bootstrap sample
    imp.boot <- as.data.frame(lapply(listDF[[i]], function(x){x[boot.sample]})) 
    
    # weight bootstrap sample
    f1 <- crprep(Tstop=imp.boot$TTE_CVDDeath, status=imp.boot$status, 
                 keep=names(imp.boot),
                 # shorten=FALSE,
                 data=imp.boot) %>% select(-c(4))
    
    
    # fit model
    fit.competing.boot <- flexsurvspline(Surv(Tstart,Tstop,status==1)~
                                           cAge+Gender+cBMI+CurDrink+Smoker+
                                           TreatHyperten_adj+Choltreat+
                                           DiabetesDurCat+AntiDiab_adj+
                                           cHbA1c+cTOTC+cHDLC,
                                         k=0,
                                         scale="hazard",
                                         data=f1,
                                         weights=weight.cens)
    
    #calculate performance in bootstrap data
    competing.boot <- calculate_performance(model = fit.competing.boot,
                                            time = imp.boot$TTE_CVDDeath, 
                                            status = imp.boot$status,
                                            data.used = imp.boot,
                                            time.horizon = time.horizon,
                                            primary.event = 1)
    
    
    #calculate performance in original data
    competing.test <- calculate_performance(model = fit.competing.boot,
                                            time = listDF[[i]]$TTE_CVDDeath, 
                                            status = listDF[[i]]$status,
                                            data.used = listDF[[i]],
                                            time.horizon = time.horizon,
                                            primary.event = 1)
    
    #Subtract performance in original data from bootstrap data to calculate optimism
    boot.eachfold[j,] <- competing.boot[[1]] - competing.test[[1]] 
    
    print(paste0("iteration done: ", j))
  }
}
end_time <- Sys.time()
end_time - start_time

mean.optimism.competing <- apply(boot.eachfold, 2, mean)


# Internal validation via 10-fold cross-validation ------------------------

set.seed(202405)

n.folds <- 10
cv.eachfold <- list()
cv.competing <- matrix(NA, n_imp, 8)

for (i in 1:n_imp){ #Loop through each imputation datasets
  
  #Create folds
  folds <- createFolds(listDF[[i]]$ID, k = 10, list = TRUE, returnTrain = FALSE)
  
  for(j in 1:n.folds){ #Loop through each fold
    
    dev.folds <- unlist(folds[-j])
    test.fold <- unlist(folds[j])
    
    # create sample
    dev.dat <- listDF[[i]][dev.folds,]
    
    # weight sample
    f1 <- crprep(Tstop=dev.dat$TTE_CVDDeath, status=dev.dat$status, 
                 keep=names(dev.dat),
                 data=dev.dat) %>% select(-c(4))
    
    fit.competing.cv <- flexsurvspline(Surv(Tstart,Tstop,status==1)~
                                         cAge+Gender+cBMI+CurDrink+Smoker+
                                         TreatHyperten_adj+Choltreat+
                                         DiabetesDurCat+AntiDiab_adj+
                                         cHbA1c+cTOTC+cHDLC,
                                       k=0,
                                       scale="hazard",
                                       data=f1,
                                       weights=weight.cens)
    
    # predict in test data
    test.dat <- listDF[[i]][test.fold,]
    
    cv.eachfold[[j]] <- calculate_performance(model = fit.competing.cv,
                                              time = test.dat$TTE_CVDDeath, 
                                              status = test.dat$status,
                                              data.used = test.dat,
                                              time.horizon = time.horizon,
                                              primary.event = 1)
    rm(dev.dat,test.dat,f1,fit.competing.cv)
    
  }
  cv.indicators <- sapply(cv.eachfold,`[[`,1)
  cv.competing[i,] <- apply(cv.indicators, 1, mean)
  
  print(paste0("imputation done: ", i))
}

cv.oe <- pool_oe(est=cv.competing[,1],
                 se=cv.competing[,2],
                 n=n_pt)

cv.cal.int <- pool_cal(est=cv.competing[,3],
                       se=cv.competing[,4],
                       n=n_pt)

cv.cal.slope <- pool_cal(est=cv.competing[,5],
                         se=cv.competing[,6],
                         n=n_pt)

cv.cindx.W <- pool_cindx(cindx=cv.competing[,7],
                         se=cv.competing[,8],
                         n=n_pt,
                         method="W")

cv.res <- rbind("OE ratio" = cv.oe,
                "Calibration intercept" = cv.cal.int,
                "Calibration slope" = cv.cal.slope,
                "Wolber's C-index" = cv.cindx.W)


# Internal-external CV -------------------

clusters <- list("HRS","HealthABC","SHARE")
N.clust <- length(clusters)
data.in <- data.leftout <- list()

`%out%` <- Negate(`%in%`)

leftout.performance.competing <- list("Indicators"=list(),
                                      "Pseudo values"=list(),
                                      "Smoothed pseudo values"=list())

for (i in 1:N.clust){ #Loop through each cluster
  
  leftout.performance.competing.each <- list()
  for(j in 1:n_imp){ #Loop through each imputation dataset
    data.in <- listDF[[j]][listDF[[j]]$Study %out% clusters[[i]],]
    data.leftout <- listDF[[j]][listDF[[j]]$Study %in% clusters[[i]],]
    n_pt_leftout <- length(data.leftout$ID)
    
    # weight sample
    f1 <- crprep(Tstop=data.in$TTE_CVDDeath, status=data.in$status, 
                 keep=names(data.in),
                 data=data.in) %>% select(-c(4))
  
    fit.competing.cv <- flexsurvspline(Surv(Tstart,Tstop,status==1)~
                                         cAge+Gender+cBMI+CurDrink+Smoker+
                                         TreatHyperten_adj+Choltreat+
                                         DiabetesDurCat+AntiDiab_adj+
                                         cHbA1c+cTOTC+cHDLC,
                                       k=0,
                                       scale="hazard",
                                       data=f1,
                                       weights=weight.cens)
    
    leftout.performance.competing.each[[j]] <- calculate_performance(
      model = fit.competing.cv,
      time = data.leftout$TTE_CVDDeath, 
      status = data.leftout$status,
      data.used = data.leftout,
      time.horizon = time.horizon,
      primary.event = 1
    )
    
    rm(fit.competing.cv)
  }
  
  leftout.indicators <- sapply(leftout.performance.competing.each,`[[`,1)
  
  leftout.oe <- pool_oe(est=leftout.indicators[1,],
                        se=leftout.indicators[2,],
                        n=n_pt_leftout)
  
  leftout.cal.int <- pool_cal(est=leftout.indicators[3,],
                              se=leftout.indicators[4,],
                              n=n_pt_leftout)
  
  leftout.cal.slope <- pool_cal(est=leftout.indicators[5,],
                                se=leftout.indicators[6,],
                                n=n_pt_leftout)
  
  leftout.cindx.W <- pool_cindx(cindx=leftout.indicators[7,],
                                se=leftout.indicators[8,],
                                n=n_pt,
                                method="W")
  
  leftout.performance.competing[[1]][[i]] <- c("OE ratio" = leftout.oe,
                                               "Calibration intercept" = leftout.cal.int,
                                               "Calibration slope" = leftout.cal.slope,
                                               "Wolber's C-index" = leftout.cindx.W)
  leftout.performance.competing[[2]][[i]] <- sapply(leftout.performance.competing.each,`[[`,2)
  leftout.performance.competing[[3]][[i]] <- sapply(leftout.performance.competing.each,`[[`,3)
}

# performance per cluster
leftout.performance.competing[[1]]

res.IECV <- bind_rows(leftout.performance.competing[[1]])
rownames(res.IECV) <- c("HRS","HealthABC","SHARE")

# Meta-analyse results
meta.oe <- metagen(TE=log(res.IECV$OE.ratio.Estimate),
                   seTE=res.IECV$OE.ratio.logSE,
                   studlab=res.IECV$NA.,
                   title="Observed-to-expected ratio",
                   n.e = c(2860,386,3461),
                   func.backtransf=exp)

forest(meta.oe,
       common=F,
       random=T,
       text.random = "Summary estimate",
       rightcols = c("effect","ci"),
       rightlabs = c("Estimate","95% CI"),
       leftcols = c("studlab","n.e"),
       leftlabs = c("Leftout study","N"),
       header.line=T,
       smlab="Observed-to-expected ratio",
       pooled.totals = F
)
plot.oe <- grid.grab()


meta.calib <- metagen(TE=res.IECV$Calibration.intercept.Estimate,
                      seTE=res.IECV$Calibration.intercept.SE,
                      studlab=res.IECV$NA.,
                      title="Calibration intercept",
                      n.e = c(2860,386,3461))

forest(meta.calib,
       common=F,
       random=T,
       text.random = "Summary estimate",
       rightcols = c("effect","ci"),
       rightlabs = c("Estimate","95% CI"),
       leftcols = c("studlab","n.e"),
       leftlabs = c("Leftout study","N"),
       header.line=T,
       smlab="Calibration intercept",
       pooled.totals = F
)
plot.calib.intcpt <- grid.grab()

meta.calib.slope <- metagen(TE=res.IECV$Calibration.slope.Estimate,
                            seTE=res.IECV$Calibration.slope.SE,
                            studlab=res.IECV$NA.,
                            title="Calibration slope",
                            n.e = c(2860,386,3461))

forest(meta.calib.slope,
       common=F,
       random=T,
       text.random = "Summary estimate",
       rightcols = c("effect","ci"),
       rightlabs = c("Estimate","95% CI"),
       leftcols = c("studlab","n.e"),
       leftlabs = c("Leftout study","N"),
       header.line=T,
       smlab="Calibration slope",
       ref=1,
       pooled.totals = F
)
plot.calib.slope <- grid.grab()

func.antilogit <- function(x){exp(x)/(1+exp(x))}

meta.Wolberc <- metagen(TE=car::logit(res.IECV$Wolber.s.C.index.Estimate),
                        seTE=res.IECV$Wolber.s.C.index.logitSE,
                        studlab=res.IECV$NA.,
                        title="C-index",
                        n.e = c(2860,386,3461),
                        func.backtransf = func.antilogit
)

forest(meta.Wolberc,
       common=F,
       random=T,
       text.random = "Summary estimate",
       rightcols = c("effect","ci"),
       rightlabs = c("Estimate","95% CI"),
       leftcols = c("studlab","n.e"),
       leftlabs = c("Leftout study","N"),
       xlim=c(0,1),
       ref=0.5,
       header.line=T,
       smlab="C-index",
       pooled.totals = F
)
plot.cindex <- grid.grab()

grid.newpage()
grid.arrange(plot.oe,plot.calib.intcpt,plot.calib.slope,plot.cindex, ncol=1)


# Comparison with SCORE2-Diabetes and PREVENT ASCVD ----------------------------------------------

# Coefficients SCORE2-Diabetes (first is for men, second for women)
coef.age <- c(0.5368,0.6624)
coef.smoking <- c(0.4774,0.6139)
coef.sbp <- c(0.1322,0.1421)
coef.diabetes <- c(0.6457,0.8096)
coef.totc <- c(0.1102,0.1127)
coef.hdlc <- c(-0.1087,-0.1568)
coef.smokage <- c(-0.0672,-0.1122)
coef.sbpage <- c(-0.0268,-0.0167)
coef.diabage <- c(-0.0983,-0.1272)
coef.totcage <- c(-0.0181,-0.0200)
coef.hdlcage <- c(0.0095,0.0186)
coef.agediag <- c(-0.0998,-0.118)
coef.hba1c <- c(0.0955,0.1173)
coef.egfr <- c(-0.0591,-0.0640)
coef.egfrsq <- c(0.0058,0.0062)
coef.hba1cage <- c(-0.0134,-0.0196)
coef.egfrage <- c(0.0115,0.0169)

# Calculate linear predictor
lp.score2diabetes.men <- expression(coef.age[1]*(cAge+70-60)/5 + 
                                      coef.smoking[1]*(as.numeric(Smoker)-1) + 
                                      coef.sbp[1]*(cSBP+130/20-120/20) + 
                                      coef.diabetes[1]*1 + 
                                      coef.totc[1]*(cTOTC+5-6)/1 + 
                                      coef.hdlc[1]*(cHDLC/0.5) + 
                                      coef.smokage[1]*(cAge+70-60)/5*(as.numeric(Smoker)-1) + 
                                      coef.sbpage[1]*(cAge+70-60)/5*(cSBP+130/20-120/20) + 
                                      coef.diabage[1]*(cAge+70-60)/5*1 +
                                      coef.totcage[1]*(cAge+70-60)/5*(cTOTC+5-6)/1 + 
                                      coef.hdlcage[1]*(cAge+70-60)/5*(cHDLC/0.5) + 
                                      coef.agediag[1]*(cDiabetesAge/5) + 
                                      coef.hba1c[1]*((10.93*(cHbA1c+6.5)-23.5)-31)/9.34 + 
                                      coef.egfr[1]*clneGFR +
                                      coef.egfrsq[1]*clneGFR*clneGFR +
                                      coef.hba1cage[1]*(cAge+70-60)/5*((10.93*(cHbA1c+6.5)-23.5)-31)/9.34 + 
                                      coef.egfrage[1]*(cAge+70-60)/5*clneGFR)

lp.score2diabetes.women <- expression(coef.age[2]*(cAge+70-60)/5 + 
                                        coef.smoking[2]*(as.numeric(Smoker)-1) + 
                                        coef.sbp[2]*(cSBP+130/20-120/20) + 
                                        coef.diabetes[2]*1 + 
                                        coef.totc[2]*(cTOTC+5-6)/1 + 
                                        coef.hdlc[2]*(cHDLC/0.5) + 
                                        coef.smokage[2]*(cAge+70-60)/5*(as.numeric(Smoker)-1) + 
                                        coef.sbpage[2]*(cAge+70-60)/5*(cSBP+130/20-120/20) + 
                                        coef.diabage[2]*(cAge+70-60)/5*1 +
                                        coef.totcage[2]*(cAge+70-60)/5*(cTOTC+5-6)/1 + 
                                        coef.hdlcage[2]*(cAge+70-60)/5*(cHDLC/0.5) + 
                                        coef.agediag[2]*(cDiabetesAge/5) + 
                                        coef.hba1c[2]*((10.93*(cHbA1c+6.5)-23.5)-31)/9.34 + 
                                        coef.egfr[2]*clneGFR +
                                        coef.egfrsq[2]*clneGFR*clneGFR +
                                        coef.hba1cage[2]*(cAge+70-60)/5*((10.93*(cHbA1c+6.5)-23.5)-31)/9.34 + 
                                        coef.egfrage[2]*(cAge+70-60)/5*clneGFR)


#Coefficients for PREVENT ASCVD (first is for men, second for women)
coef.const.p <- c(-3.500655,-3.819975)
coef.age.p <- c(0.7099847,0.719883)
coef.nonhdlc.p <- c(0.1658663,0.1176967)
coef.hdlc.p <- c(-0.1144285,-0.151185)
coef.sbp.p <- c(-0.2837212,-0.0835358)
coef.sbp2.p <- c(0.3239977,0.3592852)
coef.diabetes.p <- c(0.7189597,0.8348585)
coef.smoking.p <- c(0.3956973,0.4831078)
coef.egfr.p <- c(0.3690075,0.4864619)
coef.egfr2.p <- c(0.0203619,0.0397779)
coef.antihyp.p <- c(0.2036522,0.2265309)
coef.statin.p <- c(-0.0865581,-0.0592374)
coef.antihyp.sbp.p <- c(-0.0322916,-0.0395762)
coef.statin.nonhdlc.p <- c(0.114563,0.0844423)
coef.age.nonhdlc.p <- c(-0.0300005,-0.0567839)
coef.age.hdlc.p <- c(0.0232747,0.0325692)
coef.age.sbp.p <- c(-0.0927024,-0.1035985)
coef.age.diabetes.p <- c(-0.2018525,-0.2417542)
coef.age.smoker.p <- c(-0.0970527,-0.0791142)
coef.age.egfr.p <- c(-0.1217081,-0.1671492)

# Calculate linear predictor
lp.prevent.men <- expression(coef.const.p[1] +
                               coef.age.p[1]*(cAge+70-55)/10 + 
                               coef.nonhdlc.p[1]*((cTOTC+5)-(cHDLC+1.3)-3.5) +
                               coef.hdlc.p[1]*(cHDLC)/0.3 +
                               coef.sbp.p[1]*(min(cSBP*20+130,110)-110)/20 +                                
                               coef.sbp2.p[1]*(max(cSBP*20+130,110)-130)/20 +     
                               coef.diabetes.p[1]*1 +
                               coef.smoking.p[1]*(as.numeric(Smoker)-1) + 
                               coef.egfr.p[1]*(min(exp(clneGFR*0.15+4.5),60)-60)/-15 +
                               coef.egfr2.p[1]*(max(exp(clneGFR*0.15+4.5),60)-90)/-15 +
                               coef.antihyp.p[1]*(as.numeric(TreatHyperten_adj)-1) +
                               coef.statin.p[1]*(as.numeric(Choltreat)-1) +
                               coef.antihyp.sbp.p[1]*(as.numeric(TreatHyperten_adj)-1)*(max(cSBP*20+130,110)-130)/20 +
                               coef.statin.nonhdlc.p[1]*(as.numeric(Choltreat)-1)*((cTOTC+5)-(cHDLC+1.3)-3.5) +
                               coef.age.nonhdlc.p[1]*(cAge+70-55)/10*((cTOTC+5)-(cHDLC+1.3)-3.5) +
                               coef.age.hdlc.p[1]*(cAge+70-55)/10*(cHDLC)/0.3 + 
                               coef.age.sbp.p[1]*(cAge+70-55)/10*(max(cSBP*20+130,110)-130)/20 +
                               coef.age.diabetes.p[1]*(cAge+70-55)/10*1 +
                               coef.age.smoker.p[1]*(cAge+70-55)/10*(as.numeric(Smoker)-1) + 
                               coef.age.egfr.p[1]*(cAge+70-55)/10*(min(exp(clneGFR*0.15+4.5),60)-60)/-15)

lp.prevent.women <- expression(coef.const.p[2] +
                                 coef.age.p[2]*(cAge+70-55)/10 + 
                                 coef.nonhdlc.p[2]*((cTOTC+5)-(cHDLC+1.3)-3.5) +
                                 coef.hdlc.p[2]*(cHDLC)/0.3 +
                                 coef.sbp.p[2]*(min(cSBP*20+130,110)-110)/20 +                                
                                 coef.sbp2.p[2]*(max(cSBP*20+130,110)-130)/20 +     
                                 coef.diabetes.p[2]*1 +
                                 coef.smoking.p[2]*(as.numeric(Smoker)-1) + 
                                 coef.egfr.p[2]*(min(exp(clneGFR*0.15+4.5),60)-60)/-15 +
                                 coef.egfr2.p[2]*(max(exp(clneGFR*0.15+4.5),60)-90)/-15 +
                                 coef.antihyp.p[2]*(as.numeric(TreatHyperten_adj)-1) +
                                 coef.statin.p[2]*(as.numeric(Choltreat)-1) +
                                 coef.antihyp.sbp.p[2]*(as.numeric(TreatHyperten_adj)-1)*(max(cSBP*20+130,110)-130)/20 +
                                 coef.statin.nonhdlc.p[2]*(as.numeric(Choltreat)-1)*((cTOTC+5)-(cHDLC+1.3)-3.5) +
                                 coef.age.nonhdlc.p[2]*(cAge+70-55)/10*((cTOTC+5)-(cHDLC+1.3)-3.5) +
                                 coef.age.hdlc.p[2]*(cAge+70-55)/10*(cHDLC)/0.3 + 
                                 coef.age.sbp.p[2]*(cAge+70-55)/10*(max(cSBP*20+130,110)-130)/20 +
                                 coef.age.diabetes.p[2]*(cAge+70-55)/10*1 +
                                 coef.age.smoker.p[2]*(cAge+70-55)/10*(as.numeric(Smoker)-1) + 
                                 coef.age.egfr.p[2]*(cAge+70-55)/10*(min(exp(clneGFR*0.15+4.5),60)-60)/-15)


#Calculate risks using SCORE2-Diabetes and PREVENT ASCVD and add to dataset
for (i in 1:n_imp){
  listDF[[i]] <- listDF[[i]] %>% mutate(
    lp_score2diabetes = ifelse(Gender==1,eval(lp.score2diabetes.women),eval(lp.score2diabetes.men)),
    risk_uncal = 1-ifelse(Gender==1,0.9776,0.9605)^exp(lp_score2diabetes),
    scale1 = ifelse(Region_mod==1,
                    ifelse(Gender==1,-0.3143,-0.1565),
                    ifelse(Region_high==1,
                           ifelse(Gender==1,0.5710,0.3207),
                           ifelse(Gender==1,-0.7380,-0.5699))),
    scale2 = ifelse(Region_mod==1,
                    ifelse(Gender==1,0.7701,0.8009),
                    ifelse(Region_high==1,
                           ifelse(Gender==1,0.9369,0.9360),
                           ifelse(Gender==1,0.7019,0.7476))),
    risk_cal = 1-exp(-exp(scale1+scale2*log(-log(1-risk_uncal)))), #Calibrated risk according to SCORE2-Diabetes
    
    lp_prevent = ifelse(Gender==1,eval(lp.prevent.women),eval(lp.prevent.men)),
    risk_prevent = exp(lp_prevent)/(1+exp(lp_prevent)), #Calibrated risk according to PREVENT ASCVD
    obsv_10y = ifelse(TTE_CVDDeath>10*12,0,CVD_event) 
  )
}

#Conduct 10-fold cross-validation to compare model performances

set.seed(202405)

n.folds <- 10
cv.eachfold <- list()
cv.eachfold.Score <- list()
cv.eachfold.Prevent <- matrix(NA, n_imp, 8)
cv.competing <- matrix(NA, n_imp, 8)
cv.competing.Score <- matrix(NA, n_imp, 8)
cv.competing.Prevent <- matrix(NA, n_imp, 8)

for (i in 1:n_imp){
  
  folds <- createFolds(filter(listDF[[i]])$ID, k = 10, list = TRUE, returnTrain = FALSE)
  
  for(j in 1:n.folds){
    
    dev.folds <- unlist(folds[-j])
    test.fold <- unlist(folds[j])
    
    # create sample
    dev.dat <- listDF[[i]][dev.folds,]
    
    # weight sample
    f1 <- crprep(Tstop=dev.dat$TTE_CVDDeath, status=dev.dat$status,
                 keep=names(dev.dat),
                 data=dev.dat) %>% select(-c(4))
    
    fit.competing.cv <- flexsurvspline(Surv(Tstart,Tstop,status==1)~
                                         cAge+Gender+cBMI+CurDrink+Smoker+
                                         TreatHyperten_adj+Choltreat+
                                         DiabetesDurCat+AntiDiab_adj+
                                         cHbA1c+cTOTC+cHDLC,
                                       k=0,
                                       scale="hazard",
                                       data=f1,
                                       weights=weight.cens)
    
    fit.Score <- FGR(Hist(TTE_CVDDeath,status)~lp_score2diabetes,
                     data=dev.dat,
                     cause=1)
    
    # predict in test data
    f2 <- listDF[[i]][test.fold,]
    
    ## Performance CARE-DM
    cv.eachfold[[j]] <- calculate_performance(model = fit.competing.cv,
                                              time = f2$TTE_CVDDeath,
                                              status = f2$status,
                                              data.used = f2,
                                              time.horizon = 10*12,
                                              primary.event = 1)
    
    ## Performance SCORE2-Diabetes
    cv.eachfold.Score[[j]] <- calculate_performance(model = fit.Score,
                                                    time = f2$TTE_CVDDeath,
                                                    status = f2$status,
                                                    data.used = f2,
                                                    time.horizon = 10*12,
                                                    primary.event = 1,
                                                    predicted.Score=f2$risk_cal)
    
    ## Performance PREVENT ASCVD
    
    ### OE ratio
    obj <- summary(survfit(Surv(TTE_CVDDeath, event) ~ 1, data = f2), times = 10*12)
    aj <- list("obs" = obj$pstate[, 2], "se" = obj$std.err[, 2])
    OE <- aj$obs / mean(listDF[[1]]$risk_prevent)
    OE.se <- aj$se / aj$obs
    
    ### Calibation and C-index
    perf.Prevent <- val.prob.ci.2(f2$risk_prevent,f2$obsv_10y)
    
    N <- length(f2$ID)
    s <- sum(f2$CVD_event)
    t <- N-s
    varlogitc <- (1+(N/2-1)*(1-cindx.W)/(2-cindx.W)+(N/2-1)*cindx.W/(1+cindx.W)) / (cindx.W*(1-cindx.W)*s*t)
    cindx.W.logitse <- as.numeric(sqrt(varlogitc))
    
    
    cv.eachfold.Prevent[j,] <- c("OE" = OE,
                                 "OE.se" = OE.se,
                                 "cal.int" = perf.Prevent$Calibration$Intercept[1],
                                 "cal.int.se"= (perf.Prevent$Calibration$Intercept[3]-perf.Prevent$Calibration$Intercept[2])/3.92,
                                 "cal.slope" = perf.Prevent$Calibration$Slope[1],
                                 "cal.slope.se"= (perf.Prevent$Calibration$Slope[3]-perf.Prevent$Calibration$Slope[2])/3.92,
                                 "cindex" = perf.Prevent$Cindex[1],
                                 "cindex.logitse" = cindx.logitse
    )
    
    
  }
  cv.indicators <- sapply(cv.eachfold,`[[`,1)
  cv.competing[i,] <- apply(cv.indicators, 1, mean)
  
  cv.indicators.Score <- sapply(cv.eachfold.Score,`[[`,1)
  cv.competing.Score[i,] <- apply(cv.indicators.Score, 1, mean)
  
  cv.competing.Prevent[i,] <- colMeans(cv.eachfold.Prevent, na.rm = FALSE)
  
  print(paste0("imputation done: ", i))
}

#Summarize results CARE-DM
cv.oe <- pool_oe(est=cv.competing[,1],
                 se=cv.competing[,2],
                 n=n_pt)

cv.cal.int <- pool_cal(est=cv.competing[,3],
                       se=cv.competing[,4],
                       n=n_pt)

cv.cal.slope <- pool_cal(est=cv.competing[,5],
                         se=cv.competing[,6],
                         n=n_pt)

cv.cindx.W <- pool_cindx(cindx=cv.competing[,7],
                         se=cv.competing[,8],
                         n=n_pt)

cv.res <- rbind("OE ratio" = cv.oe,
                "Calibration intercept" = cv.cal.int,
                "Calibration slope" = cv.cal.slope,
                "C-index" = cv.cindx.W)


#Summarize results SCORE2-Diabetes
cv.oe.Score <- pool_oe(est=cv.competing.Score[,1],
                       se=cv.competing.Score[,2],
                       n=n_pt)

cv.cal.int.Score <- pool_cal(est=cv.competing.Score[,3],
                             se=cv.competing.Score[,4],
                             n=n_pt)

cv.cal.slope.Score <- pool_cal(est=cv.competing.Score[,5],
                               se=cv.competing.Score[,6],
                               n=n_pt)

cv.cindx.W.Score <- pool_cindx(cindx=cv.competing.Score[,7],
                               se=cv.competing.Score[,8],
                               n=n_pt,
                               method="W")

cv.res.Score <- rbind("OE ratio" = cv.oe.Score,
                      "Calibration intercept" = cv.cal.int.Score,
                      "Calibration slope" = cv.cal.slope.Score,
                      "C-index" = cv.cindx.W.Score)

#Summarize results PREVENT ASCVD
cv.oe.Prevent <- pool_oe(est=cv.competing.Prevent[,1],
                         se=cv.competing.Prevent[,2],
                         n=n_pt)

cv.cal.int.Prevent <- pool_cal(est=cv.competing.Prevent[,3],
                               se=cv.competing.Prevent[,4],
                               n=n_pt)

cv.cal.slope.Prevent <- pool_cal(est=cv.competing.Prevent[,5],
                                 se=cv.competing.Prevent[,6],
                                 n=n_pt)

cv.cindx.W.Prevent <- pool_cindx(cindx=cv.competing.Prevent[,7],
                                 se=cv.competing.Prevent[,8],
                                 n=n_pt,
                                 method="W")

cv.res.Prevent <- rbind("OE ratio" = cv.oe.Prevent,
                        "Calibration intercept" = cv.cal.int.Prevent,
                        "Calibration slope" = cv.cal.slope.Prevent,
                        "C-index" = cv.cindx.Prevent)

comp.tot <- rbind(New=cv.res,SCORE2Diabetes=cv.res.Score,PREVENT=cv.res.Prevent)


# Decision-analysis curve -------------------------------------------------
source("stdca.R") #Load code from course "Statistical Methods for Risk Prediction and Prognostic Models" 
                  # adapted to multiply imputed dataset

dca_dat <- list()
for (i in 1:n_imp){
  pred.surv <- as.matrix(1-predict(fit.model[[i]],newdata=listDF[[i]],type="survival",times=time.horizon)[,2])
  dca_dat[[i]] <- cbind(listDF[[i]],pred.surv)
}


dca_results <- stdca(data=dca_dat, outcome="status", ttoutcome="TTE_CVDDeath", timepoint=time.horizon, 
                     predictors=c(".pred_survival"), xstop=0.5, loess.span=0.2, smooth=TRUE,
                     cmprsk=T)
dca_results


