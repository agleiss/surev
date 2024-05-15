"f.surev"<-
  function(coxfit)
    # update of f.surev() from "surev" package for R-versions starting from 3.0.0
    # by taking account of new structur of cph() output object in lines 71-76
    # by Hanna Sinkovec (2019)
  {
    #function which calculates  M(t_(j)) for a given time t_(j)
    f.Mt <- function(tempo, tutti.tempi, stima.surv, tempi.evento, Stj, ind.censura, num.sogg)
    {
      Stj1 <- unique(Stj[tempi.evento == tempo])
      primo <- rep(1 - Stj1, num.sogg)
      primo[tutti.tempi <= tempo] <- 0
      secondo <- Stj1 * (1 - ind.censura)
      secondo[tutti.tempi > tempo] <- 0
      terzo <- ind.censura * (((1 - Stj1) * Stj1)/stima.surv + Stj1 * (1 - Stj1/stima.surv))
      terzo[tutti.tempi > tempo] <- 0
      terzo[is.na(terzo)] <- 0
      ris <- primo + secondo + terzo
      return(sum(ris)/num.sogg)
    }
    #function which calculates M(t_(j)|x) for a given time t_(j)
    f.Mt.cox <- function(tempo, tutti.tempi, stima.surv, tempi.evento, Stj0, ind.censura, num.sogg, lin.pred)
    {
      #S(t_(j)|x_i)
      Stj00 <- unique(Stj0[tempi.evento == tempo])
      Stj1 <- Stj00^exp(lin.pred)
      primo <- 1 - Stj1
      primo[tutti.tempi <= tempo] <- 0
      secondo <- Stj1 * (1 - ind.censura)
      secondo[tutti.tempi > tempo] <- 0
      terzo <- ind.censura * (((1 - Stj1) * Stj1)/stima.surv + Stj1 * (1 - Stj1/stima.surv))
      terzo[tutti.tempi > tempo] <- 0
      terzo[is.na(terzo)] <- 0
      ris <- primo + secondo + terzo
      return(sum(ris)/num.sogg)
    }
    #function that evaluates the indexes to be used in the attribution, for each time of interest, 
    #of survival estimates 
    #obtained from Kaplan-Meier and Cox regression 
    f.assegna.surv <- function(tempo, tempi.eventi)
    {
      if(any(tempo == tempi.eventi)) {
        pos <- (c(1:length(tempi.eventi)) * as.numeric(tempo == 
                                                         tempi.eventi))
        pos <- pos[pos != 0]
      }
      else {
        tmp <- (tempo - tempi.eventi)
        if(all(tmp < 0))
          pos <- NA
        else {
          tmp[tmp < 0] <- Inf
          pos <- order(tmp)[1]
        }
      }
      return(pos)
    }
    #beginning of function 
    tsurv <- as.numeric(coxfit$y[, 1])
    surv <- as.numeric(coxfit$y[, 2])
    num.sogg <- length(tsurv)
    km <- survfit(Surv(tsurv, surv) ~ 1)
    tempi.eventi <- km$time[km$n.event != 0]
    pos.surv <- apply(as.matrix(tsurv), 1, f.assegna.surv, tempi.eventi)
    surv.tot.km <- (km$surv[km$n.event != 0])[pos.surv]
    surv.tot.km[is.na(surv.tot.km)] <- 1
    ind.censura <- as.numeric(!as.logical(surv))
    surv.tj <- km$surv[km$n.event != 0]
    Mt <- apply(as.matrix(tempi.eventi), 1, f.Mt, tsurv, surv.tot.km, tempi.eventi, surv.tj, ind.censura, num.sogg)
    numero.eventi <- km$n.event[km$n.event != 0]
    #if(length(coxfit$time) - 1 == length(tempi.eventi)) indici.da.eliminare <- 1 else indici.da.eliminare <- c(1, length(coxfit$time))
    #surv0.tj.cox <- coxfit$surv[ - indici.da.eliminare]
    #surv0.tot.cox <- (coxfit$surv[ - indici.da.eliminare])[pos.surv]
    surv0.tj.cox <- coxfit$surv[which(coxfit$time%in%tempi.eventi)]
    surv0.tot.cox <- (coxfit$surv[which(coxfit$time%in%tempi.eventi)])[pos.surv]
    surv.tot.cox <- surv0.tot.cox^exp(coxfit$linear.predictors)
    Mtx <- apply(as.matrix(tempi.eventi), 1, f.Mt.cox, tsurv, surv.tot.cox, tempi.eventi, surv0.tj.cox, ind.censura, num.sogg, coxfit$linear.predictors)
    Gkm <- survfit(Surv(tsurv, ind.censura) ~ 1)
    tempi.censure <- Gkm$time[Gkm$n.event != 0]	
    #modifica per tenere conto di possibilita' di dati senza censure
    if(!length(tempi.censure)) cens.tot.km <- rep(1, length(tempi.eventi))else {pos.surv.censure <- apply(as.matrix(tempi.eventi), 1, f.assegna.surv, tempi.censure)
    cens.tot.km <- (Gkm$surv[Gkm$n.event != 0])[pos.surv.censure]
    cens.tot.km[tempi.eventi < min(Gkm$time[Gkm$n.event != 0])] <- 1
    }
    pesi <- numero.eventi/cens.tot.km
    peso.tot <- sum(pesi)
    D <- sum(Mt * pesi)/peso.tot
    Dx <- sum(Mtx * pesi)/peso.tot
    V <- (D - Dx)/D
    if(any(Mt == 0)) {
      Vw <- sum((Mt[Mt != 0] - Mtx[Mt != 0])/Mt[Mt != 0] * pesi[Mt != 
                                                                  0])/sum(pesi[Mt != 0])}else Vw <- sum((Mt - Mtx)/Mt * pesi)/peso.tot
    
    return(list(Model = coxfit$call, D=D, Dx=Dx, V=V, Vw=Vw))
  }

#library(survminer)
#data(myeloma)
#coxfit<-cph(Surv(time, event) ~ CCND1 + CRIM1, y=TRUE, surv=TRUE, method="breslow", type="kaplan-meier", data=myeloma)
#Result<-f.surev(coxfit)	                                          
