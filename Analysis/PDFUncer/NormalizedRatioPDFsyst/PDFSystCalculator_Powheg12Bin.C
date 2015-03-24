#include <iostream>
#include <set>
#include <fstream>
#include <iomanip>

#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<THStack.h>
#include<TGraphErrors.h>
#include<TGraphAsymmErrors.h>
#include<TCanvas.h>
#include<TFrame.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>

#include <time.h>
#include "DataFormats.h"

//int DoPDFsyst();

//const int nptBins = 18;
//const double xbins_pt[nptBins+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};
const int nptBins = 12;
const double xbins_pt[nptBins+1] = {0.0, 7.5, 12.5, 17.5, 30, 40, 50, 70, 110, 150, 190, 250, 600};

void PDFSystCalculator_Powheg12Bin()
{

  cout << fixed << setprecision(5);
  //variable related PDF systematic uncertainty
  double weightedSelectedEvents[12][53];
  double weighted2SelectedEvents[12][53];
  
  double weightedSelectedEventsWplus[12][53];
  double weightedSelectedEventsWminus[12][53];
  
  double events_central[nptBins] = {0.0};
  double events2_central[nptBins] = {0.0};
  double wa[nptBins] = {0.0};
  double wb[nptBins] = {0.0};
  double wplus[nptBins] = {0.0};
  double wminus[nptBins] = {0.0};

  char tmpName[30];


  ifstream Wplus_file("InputWeights/WpToMuNu_12BinWeightedSelectedEvents.txt");
  ifstream Wminus_file("InputWeights/WmToMuNu_12BinWeightedSelectedEvents.txt");

  // Initialize variable
  for(int iBin(0); iBin<nptBins; iBin++)
  {
    {
      for(unsigned int j=0; j<53; j++)
      {
	weightedSelectedEvents[iBin][j] = 0;
	weightedSelectedEventsWplus[iBin][j] = 0;
	weightedSelectedEventsWminus[iBin][j] = 0;
//	cout << "initialized weightedSelectedEvents : " << weightedSelectedEvents[iBin][j] << endl;
      }
    }
  }


  TString tmp;
  int bins;
  int Npdf;
  
  int pdf=0;

  double AllWeightsWplus[636];
  double AllWeightsWminus[636];

   //while(Wplus_file >> bins >> tmp >> Npdf >> tmp>>  AllWeights[pdf] )
   while(Wplus_file >>AllWeightsWplus[pdf] )
      {
        // cout << AllWeightsWplus[pdf] << endl;
         pdf++;
         if(pdf==636)break;
       }

   pdf=0;
   while(Wminus_file >>  AllWeightsWminus[pdf] )
      {
         //cout << AllWeightsWminus[pdf] << endl;
         pdf++;
         if(pdf==636)break;
       }

  
  // Save weights each bin
    for(int iBin(0); iBin<nptBins; iBin++)
    {
	//for(int j=0; j<weights_CT10->size(); j++)
	for(unsigned int j=0; j<53; j++)
	{

	  weightedSelectedEventsWplus[iBin][j] = AllWeightsWplus[iBin*53 + j] ;
	  weightedSelectedEventsWminus[iBin][j] = AllWeightsWminus[iBin*53 + j] ;
	   
          //weightedSelectedEventsWplus[iBin][j] = AllWeightsWplus[j] ;
	  //weightedSelectedEventsWminus[iBin][j] = AllWeightsWminus[j] ;
	  
	  
	  //weightedSelectedEvents[iBin][j] = AllWeightsWplus[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWminus[iBin*53 + j] ;
	  weightedSelectedEvents[iBin][j] = AllWeightsWplus[iBin*53 + j] + AllWeightsWminus[iBin*53 + j] ;
	  
	 
	  cout<<iBin<<"\t "<<j<<"\t"<<weightedSelectedEvents[iBin][j]<<endl;
	  
	  //cout << Form("weightedSelectedEvents[%d][%d] : %f \t weighted2SelectedEvents[%d][%d] : %f \t weights_CT10[%d] : %f ",iBin,j,weightedSelectedEvents[iBin][j],iBin,j,weighted2SelectedEvents[iBin][j], j, (*weights_CT10)[j]) << endl;
      }
    }
  

  double weightedSelectedEventsCentTot=0;
  for(int iBin(0); iBin<nptBins; iBin++)
  {
     weightedSelectedEventsCentTot += weightedSelectedEvents[iBin][0];
  }
    // Calculate PDF syst
  for(int iBin(0);iBin<nptBins;iBin++)
  {
    //unsigned int nmembers = weights_CT10->size();
    unsigned int nmembers = 52;
    unsigned int npairs = (nmembers-1)/2;
 
    events_central[iBin] = weightedSelectedEvents[iBin][0]/weightedSelectedEventsCentTot/(xbins_pt[iBin+1]-xbins_pt[iBin]);
    if(npairs>0){
      for (unsigned int j=0; j<npairs; ++j)
      {
	double weightedSelectedEventsXPlusTot=0;
	double weightedSelectedEventsXMinusTot=0;
	for(int iBin(0);iBin<nptBins;iBin++)
	{
	  weightedSelectedEventsXPlusTot += weightedSelectedEvents[iBin][2*j+1];
	  weightedSelectedEventsXMinusTot += weightedSelectedEvents[iBin][2*j+2];
	}

	//cout << "events central : " << events_central[iBin] << endl;
	//wa[iBin] = weightedSelectedEvents[iBin][2*j+1]/events_central[iBin]-1.;
	//wb[iBin] = weightedSelectedEvents[iBin][2*j+2]/events_central[iBin]-1.;
	
	wa[iBin] = (weightedSelectedEvents[iBin][2*j+1]/weightedSelectedEventsXPlusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
	wb[iBin] = (weightedSelectedEvents[iBin][2*j+2]/weightedSelectedEventsXMinusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
	
	if (wa[iBin]>wb[iBin]){
	  if (wa[iBin]<0.) wa[iBin] = 0.;
	  if (wb[iBin]>0.) wb[iBin] = 0.;
	  wplus[iBin] += wa[iBin]*wa[iBin];
	  wminus[iBin] += wb[iBin]*wb[iBin];
	}else{
	  if (wb[iBin]<0.) wb[iBin] = 0.;
	  if (wa[iBin]>0.) wa[iBin] = 0.;
	  wplus[iBin] += wb[iBin]*wb[iBin];
	  wminus[iBin] += wa[iBin]*wa[iBin];
	}
      }
      if (wplus[iBin]>0) wplus[iBin] = sqrt(wplus[iBin]);
      if (wminus[iBin]>0) wminus[iBin] = sqrt(wminus[iBin]);
    }else{
      cout << "\tNO eigenvectors for uncertainty estimation" << endl;
    }
    //cout <<iBin+1<<" Bin: Relative uncertainty with respect to central member: +" << 100.*wplus[iBin] << " / -" <<  100.*wminus[iBin] << " [%]" << endl;
    cout <<iBin+1<<" + " << 100.*wplus[iBin] << " / -" <<  100.*wminus[iBin] << " [%]" << endl;
  }

 // TFile* Hist_out = new TFile("ZptPreFSR_Powheg12Bin.root","recreate");

}
