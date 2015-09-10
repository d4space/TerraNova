#include <iostream>
#include <TString.h>
#include <TRandom.h>

static const int NB=18;

const double db[NB]={2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 10, 10, 10, 20, 20, 20, 40, 40, 60, 350}; 

const double nn[NB]={682.757, 1131.67, 1061.95, 789.873, 726.111, 492.814, 460.524, 350.853, 961.589, 
		     532.466, 329.296, 352.747, 145.448, 67.082, 54.0626, 25.0404, 4.40732, 12.7016};

const double dn[NB]={37.7383, 51.9091, 50.8778, 45.6727, 42.9076, 37.0839, 34.9406, 31.1642, 39.2429, 
		     29.216, 23.3506, 22.9646, 15.171, 10.5114, 9.04859, 6.1784, 3.0343, 5.60278};

double ff[NB];
double df[NB];
double dp[NB];

void ErrPropaNormXsec(){

    ///Here is error propagation
    double NN=0;
    for(int i=0; i<NB; ++i) NN+=nn[i];
    for(int i=0; i<NB; ++i) {
	ff[i]=nn[i]/NN;
	dp[i]=dn[i]/NN;
	df[i]=0;
    }
    for(int i=0; i<NB; ++i){
	df[i] = pow( (1-ff[i])*dp[i], 2);
	for(int j=0; j<NB; ++j){
	    if(j!=i) df[i]+=ff[i]*ff[i]*dp[j]*dp[j];
	}
	df[i]=sqrt(df[i]);
    }


    //Now Toy
    double mx2[NB];
    double m2x[NB];
    for(int i=0; i<NB; ++i){
	mx2[i]=0;
	m2x[i]=0;
    }
    double temp[NB];
    for(int e=0; e<100000; ++e){
	for(int i=0; i<NB; ++i) temp[i]=gRandom->Gaus(nn[i], dn[i]);
	double sum=0;
	for(int i=0; i<NB; ++i) sum+=temp[i];
	for(int i=0; i<NB; ++i) {
	    temp[i]/=sum;
	    m2x[i]+=temp[i];
	    mx2[i]+=temp[i]*temp[i];
	}
    }
    double rms[NB];
    for(int i=0; i<NB; ++i) {
	m2x[i]/=100000;
	mx2[i]/=100000;
	rms[i]=sqrt(mx2[i]-m2x[i]*m2x[i]);
    }


    for(int i=0; i<NB; ++i){
	cout << Form("%.8f  prop=%.8f  rms=%.8f", ff[i]/db[i], df[i]/db[i], rms[i]/db[i]) << endl;
    }
    
}
