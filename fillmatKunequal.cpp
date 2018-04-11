#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fillMat(NumericVector data, NumericVector group, int gsize, int K){
	int n = data.size();
	
	NumericMatrix outmat(gsize, K);
  NumericVector count(K);
  
	for(int j = 0; j < n; j++){
	  for(int k = 0; k < K; k++){
	    if((group[j]-1 != k) && (count[k] < gsize)){
				outmat(count[k],k) = data[j];
	      count[k] = count[k] + 1;
	    }
	  }
	}
		
	return(outmat);
}

// [[Rcpp::export]]
NumericMatrix sumCovsK(NumericVector covariates, NumericVector covsums, NumericVector group, int gsize, int K){
  int n = covariates.size();
  
  NumericMatrix outmat(gsize, K);
  NumericVector count(K);

  for(int k = 0; k < K; k++){
    outmat(0, k) = covsums[k];
    for(int i = 0; i < n; i++){
      if(group[i]-1 != k){
        count[k] = count[k] + 1;
        //outmat(i,count) = outmat(i,count-1) - covariates[j];
        if(count[k] < gsize){
          outmat(count[k],k) = outmat(count[k]-1,k) - covariates[i];
        }
      }
      
    }
  }
  return(outmat);
}

// [[Rcpp::export]]
NumericMatrix haztimeMat(NumericMatrix tm, int K, NumericVector group){
  int n = group.size();
  int countmax = tm.nrow();
  
  NumericMatrix outmat(n, K);
  for(int k = 0; k < K; k++){
    int count = 0;
    for(int j = 0; j < n; j++){
      outmat(j,k) = tm(count,k);
      if((count ==  countmax) && (group[j]-1 == k)){
        outmat(j,k) = tm(count-1, k);
      }
      if(group[j]-1 != k){
        count++;
      }
    }
  }
  return(outmat);
}

// [[Rcpp::export]]
NumericMatrix hazMat(NumericMatrix times, NumericVector scales, NumericVector shapes){
  int n = times.nrow();
  int K = times.ncol();  
  NumericMatrix outmat(n, K);
  for(int i = 0; i < n; i++){
    for(int k = 0; k < K; k++){
      outmat(i,k) = (shapes[k]/scales[k])*pow(times(i,k)/scales[k], shapes[k]-1);

    }
  }
  return(outmat);
}

// [[Rcpp::export]]
List fixHazards(NumericMatrix hazards, NumericMatrix evmat, int K, bool raw){
  int n = hazards.nrow();
  NumericMatrix hazfix(n, K);
  NumericMatrix useEv(n, K);
  for(int k = 0; k < K; k++){
    double currmax = 0.0;
    for(int j = 0; j < n; j++){
      if(hazards(j,k) <= currmax){
        hazfix(j,k) = currmax;
        if(currmax == 0.0){
          useEv(j,k) = 0;
          if(raw){
            hazfix(j,k) = 0;
          }else{
            hazfix(j,k) = 1;
          }
        }else{
          useEv(j,k) = evmat(j,k);
        }
      }else{
        hazfix(j,k) = hazards(j,k);
        currmax = hazards(j,k);
        useEv(j,k) = evmat(j,k);
      }
    }
  }
  
  return List::create(Named("useEv", useEv), Named("hazfix", hazfix));
  
}

// [[Rcpp::export]]
List fixHazards2(NumericMatrix hazards, NumericMatrix evmat, int K, bool raw){
  int n = hazards.nrow();
  NumericMatrix hazfix(n, K);
  NumericMatrix useEv(n, K);
  for(int k = 0; k < K; k++){
    double currmax = 0.0;
    for(int j = 0; j < n; j++){
      if(hazards(j,k) <= currmax){
        hazfix(j,k) = currmax;
        if(currmax == 0.0){
          useEv(j,k) = 0;
          if(raw){
            hazfix(j,k) = 0;
          }else{
            hazfix(j,k) = 0;
          }
        }else{
          useEv(j,k) = evmat(j,k);
        }
      }else{
        hazfix(j,k) = hazards(j,k);
        currmax = hazards(j,k);
        useEv(j,k) = evmat(j,k);
      }
    }
  }
  
  return List::create(Named("useEv", useEv), Named("hazfix", hazfix));
  
}


// [[Rcpp::export]]
double addMatVals(NumericMatrix mat, NumericVector groups, int useable){
  //int n = groups.size();
  double outval = 0.0;
  
  for(int i = 0; i < useable; i++){
    outval = outval + mat(i, groups[i]-1);
  }
  
  
  return(outval);
}

// [[Rcpp::export]]
NumericVector genpexp(NumericVector ExpVals, NumericVector rates, NumericVector times){
  int n = ExpVals.size();
  int R = rates.size();
  NumericVector ret(n);
  
  
  for(int i = 0; i < n; i++){
    double cumVal = 0.0;
    double cumVal2 = 0.0;
    bool stop = FALSE;
    int count = 0;
    while(!stop){
      if(count > 0){
        cumVal = cumVal + rates[count]*(times[count] - times[count-1]); 
      }else{
        cumVal = cumVal + rates[count]*times[count];
      }
      if(ExpVals[i] < cumVal){
        stop = TRUE;
        if(count > 0){
          ret[i] = times[count-1] + (ExpVals[i] - cumVal2)/rates[count];
        }else{
          ret[i] = ExpVals[i]/rates[0];
        }
      }else if((count == (R-1)) && (ExpVals[i] > cumVal)){
        stop = TRUE;
        ret[i] = times[count] + (ExpVals[i] - cumVal)/rates[count];
      }
      cumVal2 = cumVal; 
      count++;
    }
  }
  return(ret);
}

// [[Rcpp::export]]
NumericVector smoothHaz(NumericVector estgrid, NumericVector times, NumericVector baseHaz, double b){
  int npoints = estgrid.size();
  int n = times.size();
  
  NumericVector retvec(npoints);
  
  for(int i = 0; i < npoints; i++){
    //double currsum = 0.0;
    for(int j = 0; j < n; j++){
      double currval = (estgrid[i] - times[j])/b;
      if( fabs(currval) <= 1.0){
        retvec[i] = retvec[i] + (1/b)*(0.75)*(1-pow(currval, 2.0))*baseHaz[j];
      }
    }
    //retvec[i] = currsum/b;
  }
  return(retvec);
}

// [[Rcpp::export]]
NumericVector peHaz(NumericVector estgrid, NumericVector times, NumericVector baseHaz, NumericVector timesout, double b){
  int npoints = estgrid.size();
  int n = times.size();
  int nout = timesout.size();
  
  NumericVector pevec(npoints);
  NumericVector retvec(nout);
  double currtime = 0.0;
  
  for(int i = 0; i < npoints; i++){
    double numevents = 0.0;
    double hazval = 0.0;
    for(int j = 0; j < n; j++){
      if( (times[j] <= currtime + b) & (times[j] > currtime)){
          numevents = numevents + 1.0;
          hazval = hazval + baseHaz[j];
      }

    }
    pevec[i] = hazval/b;
    currtime = currtime + b;
  }
  
  for(int i = 0; i < nout; i++){
    for(int j = 1; j < npoints; j++){
      if((timesout[i] > estgrid[j-1]) & (timesout[i] <= estgrid[j])){
        retvec[i] = pevec[j];
      }
    }
  }
  
  retvec[n-1] = pevec[npoints-1];
  
  return(retvec);
}
