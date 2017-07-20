CrossPowerSpectrum2 = function(X, T, m, factor=1){
  ### X: time series by COLUMNs, T: time, m: window smoothing.
  ##factor: PCA component return.
  #h = 1/(2*m + 1) ##uniform weight
  n = dim(X)[2]
  x.omega = matrix(0,n,T)
  
  for(i in 1:n){
    #x.omega[i,] = (1/sqrt(T))*fft(X[,i]);
    x.omega[i,] = fft(X[,i]);
  } 
  
  paddingSpectrum = array(0, dim = c(n,n,T+2*m))
  for(i in 1:(T/2) ){
    paddingSpectrum[,,(1+m):(T+m)] = x.omega[,i+1]%*%t(Conj(x.omega[,i+1]))
  }
  for(i in 1:m){
    paddingSpectrum[,,m+1-i] = Conj(paddingSpectrum[,,i+m+1])
    paddingSpectrum[,,T+m-i] = Conj(paddingSpectrum[,,i+T+m])
  }
  #Smoothing
  Sxx = array(0, dim=c(n,n,T/2))
  for(i in 1:(T/2)){
    Sxx[,,i] = apply(paddingSpectrum[,,(m+i-m):(m+i+m)],1:2,mean)
  }
  
  
  ##Store only store 1st component
  EigenValues = matrix(0,(T),1)
  EigenVectors = matrix(0,(T),n)
  for (k in 1:(T/2)) {
    u = eigen(Sxx[,,k], symmetric=TRUE)
    EigenValues[k,1] = u$values[factor]
    EigenVectors[k,] = u$vectors[,factor]
  }
  return(list(x.omega=x.omega,EigenValues=EigenValues,EigenVectors=EigenVectors))
}

ComputeFactors = function(x.omega, EigenValues, EigenVectors){
  n = dim(x.omega)[1]
  c.omega  = matrix(0,T,n)
  for(i in 1:(T/2)){
    c.omega[i,] = EigenVectors[i,]
  }
  for(i in (T/2+1):(T)){
    c.omega[i,] = Conj(c.omega[T-i+1,])
  }
  y=NULL
  for(i in 1:(T)){
    y = c(y, Conj(c.omega[T-i+1,])%*%(x.omega[,i]))
  }
  
  Y = fft(y,inverse=TRUE)/(T)
  
  return(Y)
}

## Example
T=1000;
u1 = arima.sim(T, model=list(ar=c(0.9)));
plot(u1)

X1 = u1 + rnorm(T,0,(0.1))
X2 = -u1 + rnorm(T,0,(0.1))

plot(X1,type="l",ylim=c(min(X1, X2)-1, max(X1, X2)+1),ylab="")
lines(X2,col="blue")
m = 1/2*(X1+X2);
plot(m)
plot(cbind(mean(X1), mean(X2))); 

temp = CrossPowerSpectrum2(t(rbind(X2,X1)), T, m =8,factor = 1)
X.new = Re(ComputeFactors(temp$x.omega,temp$EigenValues,temp$EigenVectors))
lines(X.new, col="red")
plot(-X.new, type="l");
lines(u1, col="2")
dff = -X.new - u1
plot(dff, type="l")

fftd = fft(dff);
  pd = (abs(fftd))^2
plot(pd, type="l");

## white noise 
noise = rnorm(T);
fftd = fft(noise);
pd = (abs(fftd))^2
plot(pd, type="l");


temp = CrossPowerSpectrum2(t(rbind(X2,X1)), T, m =8,factor = 2)
X.new = Re(ComputeFactors(temp$x.omega,temp$EigenValues,temp$EigenVectors))
lines(X.new, col="orange")

temp = CrossPowerSpectrum2(t(rbind(X2,X1)), T, m =8,factor = 1)
X.new = Re(ComputeFactors(temp$x.omega,temp$EigenValues,temp$EigenVectors))

