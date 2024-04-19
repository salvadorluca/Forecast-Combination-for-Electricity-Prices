
for.perf<-function(serie,previsioni)
{
  y.orig<-serie
  previ<-previsioni
  hor<-length(previsioni)
  n.obs<-length(serie)
  
  oss.out<-y.orig[(n.obs-hor+1):n.obs]   # valori osservati sull'orizzonte previsivo: y_(T-h+1),...,y_(T))
  err.p<-oss.out-previ  # errori di previsione
  EM<-mean(err.p)       # Errore Medio
  VE<-var(err.p)        # Varianza errori di previsione
  EQM<-mean(err.p^2)    # Mean Squared Error 
  # Scomposizione di EQM
  ES<-EM^2              # Errore Sistematico
  sd.previ<-sqrt(var(previ)*((hor-1)/hor))   # Standard deviation previsioni
  sd.oss<-sqrt(var(oss.out)*((hor-1)/hor))  # Standard deviation valori osservati
  EV<-(sd.previ-sd.oss)^2                   # Errore in varianza
  EC<-2*(1-cor(previ,oss.out))*sd.previ*sd.oss    # Errore in covarianza
  
  q.ES<-ES/EQM
  q.EV<-EV/EQM
  q.EC<-EC/EQM  
  
  EAM<-mean(abs(err.p))
  EQM<-mean((err.p)^2)  # RMSE
  RMSE<-sqrt(EQM)  # RMSE
  MSE=EQM
  # Errori percentuali
  p.err<-err.p/oss.out*100
  RMSPE<-sqrt(mean(p.err^2))
  MAPE<-mean(abs(p.err))
  
  # Calcoli per indice di Theil
  p.naive<-y.orig[(n.obs-hor):(n.obs-1)]   # Previsioni naive
  e.naive<-oss.out-p.naive                 # Errori Previsioni naive
  EQM.naive<-mean(e.naive^2)
  
  U.theil<-sqrt(EQM)/sqrt(EQM.naive)
  smape =(1/length(previ) * sum(2*abs(previ-oss.out) / (abs(oss.out)+abs(previ))*100))
  mase=mase(oss.out,previ)
  mase2 <- mean(abs(err.p)) / (mean(abs(diff(oss.out))) / (n.obs - hor))
  
  seasonal_diff <- abs(diff(oss.out, lag = 7))
  mase_denominator <- mean(seasonal_diff, na.rm = TRUE)
  mase3 <- (mean(abs(err.p)) / (mase_denominator))/(hor-7)
  
  
  return(list(EAM=EAM,RMSE=RMSE,q.ES=q.ES,q.EV=q.EV,q.EC=q.EC,RMSPE=RMSPE,MAPE=MAPE,U.theil=U.theil,smape=smape,mase=mase,mase2=mase2, mase3=mase3,MSE=EQM))
  
  }