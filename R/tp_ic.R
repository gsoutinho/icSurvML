#' Estimation of Transition Probabilities in Progressive Illness-Death Models
#'
#' @description Computes transition probability estimators for progressive illness-death
#' models under a landmarking framework, including the Aalen–Johansen estimator
#' (AJ), the landmark estimator based on Kaplan–Meier with midpoint imputation
#' (LM), and the proposed Landmark Turnbull estimator (LMTB) for interval-censored
#' intermediate events.
#'
#' @param data Data frame containing the multi-state data structure. The dataset
#'   should include the variables required for the selected estimator, namely:
#'   \code{time1.mean}, \code{event1}, \code{Stime}, \code{event},
#'   \code{left.rec}, \code{right.rec}, and \code{id}.
#'
#' @param s Landmark time used to condition the transition probabilities.
#'
#' @param method Character string specifying the estimation method.
#'   Available options are:
#'   \describe{
#'     \item{\code{"AJ"}}{Aalen–Johansen estimator.}
#'     \item{\code{"LM"}}{Landmark estimator based on Kaplan–Meier and midpoint imputation.}
#'     \item{\code{"LMTB"}}{Landmark Turnbull estimator for interval-censored intermediate events.}
#'   }
#'
#' @param conf Logical. If \code{TRUE}, confidence intervals are computed whenever
#'   supported by the selected method. Defaults to \code{FALSE}.
#'
#' @return
#' Depending on the selected method:
#' \describe{
#'   \item{\code{"AJ"} and \code{"LM"}}{An object returned by
#'   \code{tprob}, compatible with the corresponding
#'   \code{summary()} and \code{plot()} methods.}
#'
#'   \item{\code{"LMTB"}}{A data frame of class \code{"flexIC_LMAJ"}
#'   containing the estimated transition probabilities over time:
#'   \code{p11(s,t)}, \code{p12(s,t)}, \code{p13(s,t)},
#'   \code{p22(s,t)}, and \code{p23(s,t)}.}
#' }
#'
#' @details
#' The \code{"LMTB"} method implements the subsampling landmark strategy
#' proposed for interval-censored intermediate events. Turnbull intervals are
#' estimated through \code{get_intervals()}, and survival probabilities are
#' obtained using \code{survival_ic()}.

tp_ic<- function(data, s=s, method= "AJ", conf = FALSE){
  
  db<-data
  
  #s=365
  
  conf = FALSE
  
  #method= "AJ"
  
  if(method=="AJ"){
    
    res.aj <- tprob(survIDM(time1.mean, event1, Stime, event) ~ 1, s = s,
                    method = "AJ", conf = FALSE, data = db)
    #res.aj$est
    #res.aj.f<-summary(res.aj, time=seq(s, max(db$Stime), seq=1))
    #res.aj.f2<-res.aj.f$est[2]
    
    return(res.aj)
    
  }
  
  
  if(method=="LM"){
    
    res.lm <- tprob(survIDM(time1.mean, event1, Stime, event) ~ 1, s = s,
                    method = "LM", conf = FALSE, data = db)
    
    return(res.lm)
    
  }

  if(method=="LMTB"){
    
    db$prob.inc_ini<-ifelse(s<=db$left.rec, 1, ifelse(s>=db$right.rec, 0, NA))
    
    db$prob.inc<-ifelse(is.na(db$prob.inc_ini), round((db$right.rec-s)/
                        (db$right.rec-db$left.rec),4), db$prob.inc_ini)
    
    m<-round(sum(db$prob.inc, na.rm = T),0); m
    
    db$prob.inc_f<-db$prob.inc/m
    
    #bootstrap simulation
    
    M<-1 #500   
    
    times <- seq(s, max(db$Stime), seq=1)#max(db$Stime) #c(1,2,3,4,5,6,7,8)
    
    EST.all<-NULL  #p11
    
    EST.TB.p13.all<-NULL  #p13
    
    db$prob.inc_f<-ifelse(is.na(db$prob.inc_f), 0, db$prob.inc_f)
    
    for(i in 1:M){
      
      subsample <- sort(sample(1:max(db$id), size = m, replace = TRUE, 
                               prob = db$prob.inc_f))#; subsample
      
      db.s<-db[subsample,]; db.s
      
      #p11:
      
      #tdata <- TNBintervals(db.s$left.rec, db.s$right.rec, nboot = 5)
      
      tbdata <- get_intervals(db.s$left.rec, db.s$right.rec, group = NULL, 
                                 nboot = 5)
      
      #est <- TNBsurvival(data = tdata, times = times, conf = FALSE, 
      #conf.level = 0.95)
      
      est<-survival_ic(tbdata, times = times, conf = TRUE)#conf.level = 0.95
      
      #est #mesmos resultados
      
      EST.all<-cbind(EST.all, est[,2])
      
      #p13:
      
      km11<- survfit(Surv(Stime, event) ~ 1, data= db.s)
      
      est.TB.p13<-as.numeric(as.data.frame(summary(km11, 
                                                   times =  times)[6])[,1])
      
      EST.TB.p13.all<-cbind(EST.TB.p13.all, est.TB.p13)
      
    }
    
    db$prob2.inc_ini<-ifelse(s<=db$left.rec, 0, ifelse(s>= db$right.rec, 1, NA))
    
    db$prob2.inc_ini<-ifelse(s>=db$Stime, 0, db$prob2.inc_ini )
    
    db$prob2.inc<-ifelse(is.na(db$prob2.inc_ini),round((s-db$left.rec)/
                      (db$right.rec-db$left.rec),4) , db$prob2.inc_ini)
    
    m2<-ifelse(sum(db$prob2.inc, na.rm = T)<1, 1, round(sum(db$prob2.inc, 
                                                            na.rm = T),0))
    
    db$prob2.inc_f<-db$prob2.inc/m2
    
    db$prob2.inc_f<-ifelse(is.na(db$prob2.inc_f), 0, db$prob2.inc_f)
    
    EST.p22.all<-NULL
    
    for(i in 1:M){
      
      subsample <- sort(sample(1:max(db$id), size = m, replace = TRUE, 
                               prob = db$prob2.inc_f))#; subsample
      
      db.s<-db[subsample,]; db.s; dim(db.s)
      
      #p22:
      
      km22<- survfit(Surv(Stime, event) ~ 1, data= db.s)
      
      est.p22<-as.numeric(as.data.frame(summary(km22, times =  times)[6])[,1])
      
      EST.p22.all<-cbind(EST.p22.all, est.p22)
      
    }
    
    #head(EST.TB.p13.all)  #p13
    #head(EST.all) #p11
    #(rowMeans(EST.TB.p13.all)-rowMeans(EST.all))[1:6]
    #(rowMeans(EST.TB.p13.all)-rowMeans(EST.all))[1:2545]  #p12
    
    lmaj_11<-as.numeric(EST.all)
    lmaj_12<-as.numeric(EST.TB.p13.all)-as.numeric(EST.all)
    lmaj_13<-1-as.numeric(EST.TB.p13.all)
    lmaj_22<-as.numeric(EST.p22.all)
    lmaj_23<-1-as.numeric(EST.p22.all)
    
    res<- data.frame(cbind(times, lmaj_11, lmaj_12, lmaj_13, lmaj_22, lmaj_23))  

    names(res) <- c('t', 'p00', 'p01', 'p02', 'p11', 'p12')
    
    result <- list(est=res)
    
    result$est
    
    res.aj$est
    
    align_locf <- function(est, t_new) {
      idx <- findInterval(t_new, est$t)          # última linha de est com t <= t_new
      out <- est[pmax(idx, 1), , drop = FALSE]   # evita índice 0
      out$t <- t_new
      out[idx == 0, setdiff(names(est), "t")] <- NA  # t_new anterior ao 1.º tempo de est
      rownames(out) <- NULL
      out
    }
    
    result_aligned <- align_locf(result$est, res.aj$est$t)
    
    #result_aligned[1:15,]
    
    #result$est[1:40,]
    
    #result_aligned$t==res.aj$t
    
    #result_aligned #nao preciso de colocar agora $est
    
    #result$est #nao tem $est
    
    #verificar
    
    # 1) Nos tempos comuns, os valores coincidem com os originais
    #common <- intersect(result$est$t, res.aj$est$t)
    #all.equal(result_aligned[result_aligned$t %in% common, -1],
    #          result$est[result$est$t %in% common, -1],
    #          check.attributes = FALSE)
    
    # 2) Mesmos tempos que res.aj
    #identical(result_aligned$t, res.aj$est$t)
    
    # 3) Tempos sem valor anterior em result (devem ser 0 neste caso, pois ambos começam em 365)
    #sum(is.na(result_aligned$p00))
    
    #  
    
    #res <- list(est = res)#object$est
    
    #res$est
    
    #names(res) <- c("t", "00", "01", "02", "11", "12")
    
    #names(res$est) <- c('t', 'p00', 'p01', 'p02', 'p11', 'p12')
    
    #class(result) <- c("LMAJ", "flexIC_tp") 
    
    #return(result)
    
    x.nlevels <- 1
    
    levels <- NULL
    
    res <- list(est = result_aligned,  s = s, t = result_aligned$t, Nlevels=x.nlevels, conf=conf)
  
    class(res)<-c("LMAJ", "survIDM")
    
    return(res)
    
  }
}