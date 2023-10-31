library(hash)

index_closed_form_formula=function(x, index_flag, pk_type){
  # x is the sample
  # index_flag can be "qZI" or "qDI"
  # pk_type is the type of the p_k sequence, can be:
  # "E" - empirical, "H" - Hazen, "WG" - Weibull and Gumbel, "HF" - Hyndman and Fan
  x=sort(x)
  n=length(x)
  k=floor((n-1)/2)
  m=k
  if (pk_type=='E'){
    if(index_flag=='qDI'){
      kf=floor((n)/2)
      result=1-2/n*sum(x[1:kf]/rev(x)[1:kf])
      if(n%%2==1){
        result=result-1/n
      }
    }
    if(index_flag=='qZI'){
      kf=floor((n)/2)
      kc=ceiling((n)/2)
      nf=kf+kf
      nc=kc+kf
      result=1-1/n*sum(x[1:kf]/x[(kf+1):nf]+x[1:kf]/x[(kc+1):nc])
      if(n%%2==1){
        result=result-x[kf+1]/x[n]/n
      }
    }
  }
  else{
    if (pk_type=='H'){
      pk_seq=((1:n)-0.5)/(n)
    }
    if (pk_type=='WG'){
      pk_seq=(1:n)/(n+1)
    }
    if (pk_type=='HF'){
      pk_seq=((1:n)-1/3)/(n+1/3)
    }
    delta_p=pk_seq[2]-pk_seq[1]
    if(index_flag=='qDI'){
      ai=x
      bi=(x[2:n]-x[1:(n-1)])/(pk_seq[2:n]-pk_seq[1:(n-1)])
      ci=pk_seq
      di=rev(x)[2:n]
      ei=(rev(x)[1:(n-1)]-rev(x)[2:n])/(rev(pk_seq)[1:(n-1)]-rev(pk_seq)[2:n])
      fi=rev(pk_seq)[2:n]
      i1=x[1]/x[n]*pk_seq[1]
      i2=sum(-bi[1:k]/ei[1:k]*(delta_p+(ai[1:k]/bi[1:k]-ci[1:k]+di[1:k]/ei[1:k]+1-fi[1:k])*log(abs((pk_seq[2:(k+1)]-di[1:k]/ei[1:k]-1+fi[1:k])/(pk_seq[1:k]-di[1:k]/ei[1:k]-1+fi[1:k])))))
      if(n%%2==0){
        i3=-bi[k+1]/ei[k+1]*(1/2-ci[k+1]+(ai[k+1]/bi[k+1]-ci[k+1]+di[k+1]/ei[k+1]+1-fi[k+1])*log(abs((0.5-di[k+1]/ei[k+1]-1+fi[k+1])/(pk_seq[k+1]-di[k+1]/ei[k+1]-1+fi[k+1]))))
      }
      else{
        i3=0
      }
      result=1-2*(i1+i2+i3)
    }
    if(index_flag=='qZI'){
      ai=x
      bi=(x[2:n]-x[1:(n-1)])/(pk_seq[2:n]-pk_seq[1:(n-1)])
      ci=pk_seq
      di=x[(k+1):n]
      ei=(x[(k+2):n]-x[(k+1):(n-1)])/(pk_seq[(k+2):n]-pk_seq[(k+1):(n-1)])
      fi=pk_seq[(k+1):n]
      if(n%%2==0){
        etl2=(x[m+2]-x[m+1])/delta_p
        etl3=(x[m+3]-x[m+2])/delta_p
        tail_left_1=x[1]/etl2*log(abs((x[m+1]/etl2+pk_seq[m+2]-pk_seq[m+1])/(x[m+1]/etl2-pk_seq[m+1]+0.5))) #poprawny tail_1
        tail_left_2=x[1]/etl3*log(abs((x[m+2]/etl3+pk_seq[1]-pk_seq[m+2]+0.5)/(x[m+2]/etl3-pk_seq[m+2]+pk_seq[m+2])))
        tail_left=tail_left_1+tail_left_2
      }
      else{
        etl2=(x[m+2]-x[m+1])/delta_p
        tail_left=x[1]/etl2*log(abs((x[m+1]/etl2+pk_seq[1]+0.5-pk_seq[m+1])/(x[m+1]/etl2-pk_seq[m+1]+0.5)))
      }
      if(n%%2==0){
        etr1=(x[m+1]-x[m])/delta_p
        etr2=(x[m+2]-x[m+1])/delta_p
        tail_right_2=(etr2/2*(0.25-(pk_seq[k+1])^2)+(x[m+1]-pk_seq[m+1]*etr2)*(1/2-pk_seq[k+1]))/x[n]
        tail_right_1=(etr1/2*((pk_seq[k+1])^2-(pk_seq[k+2+k]-1/2)^2)+(x[m]-pk_seq[m]*etr1)*(pk_seq[k+1]-pk_seq[k+2+k]+1/2))/x[n]
        tail_right=tail_right_1+tail_right_2
      }
      else{
        etr1=(x[m+1]-x[m])/delta_p
        etr2=(x[m+2]-x[m+1])/delta_p
        tail_right_2=(etr2/2*(0.25-(pk_seq[k+1])^2)+(x[m+1]-pk_seq[m+1]*etr2)*(1/2-pk_seq[k+1]))/x[n]
        tail_right_1=(etr1/2*((pk_seq[k+1])^2-(pk_seq[n]-1/2)^2)+(x[m]-pk_seq[m]*etr1)*(pk_seq[k+1]-pk_seq[n]+1/2))/x[n]
        tail_right=tail_right_1+tail_right_2
      }
      if(n%%2==0){
        int_1=bi[1:m]/ei[2:(m+1)]*((pk_seq[(m+3):(n)]-1/2-pk_seq[1:(m)])+(fi[2:(m+1)]-di[2:(m+1)]/ei[2:(m+1)]+ai[1:(m)]/bi[1:(m)]-ci[1:(m)]-0.5)*log(abs((pk_seq[(m+3):(n)]-1/2-fi[2:(m+1)]+di[2:(m+1)]/ei[2:(m+1)]+0.5)/(0.5+pk_seq[1:(m)]-fi[2:(m+1)]+di[2:(m+1)]/ei[2:(m+1)]))))
        i1=sum(int_1)
        int_2=bi[1:(m-1)]/ei[3:(m+1)]*((pk_seq[2:(m)]-pk_seq[(m+3):(n-1)]+1/2)+(fi[3:(m+1)]-di[3:(m+1)]/ei[3:(m+1)]+ai[1:(m-1)]/bi[1:(m-1)]-ci[1:(m-1)]-0.5)*log(abs((pk_seq[2:(m)]-fi[3:(m+1)]+di[3:(m+1)]/ei[3:(m+1)]+0.5)/(0.5+pk_seq[(m+3):(n-1)]-1/2-fi[3:(m+1)]+di[3:(m+1)]/ei[3:(m+1)]))))
        i2=sum(int_2)
      }
      else{
        int_1=bi[1:m]/ei[1:m]*((delta_p*(1:m)-pk_seq[1:m])+(fi[1:m]-di[1:m]/ei[1:m]+ai[1:m]/bi[1:m]-ci[1:m]-0.5)*log(abs((delta_p*(1:m)-fi[1:m]+di[1:m]/ei[1:m]+0.5)/(0.5+pk_seq[1:m]-fi[1:m]+di[1:m]/ei[1:m]))))
        i1=sum(int_1)
        int_2=bi[1:(m-1)]/ei[2:(m)]*((pk_seq[2:(m)]-delta_p*(1:(m-1)))+(fi[2:(m)]-di[2:(m)]/ei[2:(m)]+ai[1:(m-1)]/bi[1:(m-1)]-ci[1:(m-1)]-0.5)*log(abs((pk_seq[2:(m)]-fi[2:(m)]+di[2:(m)]/ei[2:(m)]+0.5)/(0.5+delta_p*(1:(m-1))-fi[2:(m)]+di[2:(m)]/ei[2:(m)]))))
        i2=sum(int_2)
      }
      result=1-2*(tail_left+i1+i2+tail_right)
    }
  }
  return(result)
}


estimator_type_labels = hash()
estimator_type_labels['E'] = 1
estimator_type_labels['H'] = 5
estimator_type_labels['WG'] = 6
estimator_type_labels['HF'] = 8
# functions qz and qd compute the values of estimator of qZ(p) and qD(p) in point p for a sample x using method est_method
qz = function(x, p, est_method){
  est_number = estimator_type_labels[[est_method]]
  result = 1-quantile(x,p/2,type=est_number)/quantile(x,(1+p)/2,type=est_number)
  return(result[[1]])
}
qd = function(x, p, est_method){
  est_number = estimator_type_labels[[est_method]]
  result = 1-quantile(x,p/2,type=est_number)/quantile(x,1-p/2,type=est_number)
  return(result[[1]])
}
# functions qzi and qdi compute the values of estimator of qZI and qDI for a sample x using method est_method
qzi = function(x, est_method){
  exception_catcher <- function(e) NA
  est_number = estimator_type_labels[[est_method]]
  qz = function(p){
    result = 1-quantile(x,p/2,type=est_number)/quantile(x,(1+p)/2,type=est_number)
    return(result[[1]])
  }
  i1=tryCatch(expr=(integrate(Vectorize(qz),0,1,subdivisions = 1000)$value),error = exception_catcher)
  return(i1)
}
qdi = function(x, est_method){
  exception_catcher <- function(e) NA
  est_number = estimator_type_labels[[est_method]]
  qd = function(p){
    result = 1-quantile(x,p/2,type=est_number)/quantile(x,1-p/2,type=est_number)
    return(result[[1]])
  }
  i1=tryCatch(expr=(integrate(Vectorize(qd),0,1,subdivisions = 1000)$value),error = exception_catcher)
  return(i1)
}



# example of use
x=rweibull(20,1,3) # sample of data
est='HF' # estimating method
qzi(x,est) # estimate using integrate function
index_closed_form_formula(x,'qZI',est) # estimate using closed-form expression



