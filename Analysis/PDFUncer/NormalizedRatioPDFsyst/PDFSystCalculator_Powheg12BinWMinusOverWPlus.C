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

void PDFSystCalculator_Powheg12BinWMinusOverWPlus()
{

  cout << fixed << setprecision(5);
  //variable related PDF systematic uncertainty
  double weightedSelectedEvents[12][53];
  double weighted2SelectedEvents[12][53];
  
  double weightedSelectedEventsWplus[12][53];
  double weightedSelectedEventsWminus[12][53];
  
  double weightedSelectedEventsWplusNorm[12][53];
  double weightedSelectedEventsWminusNorm[12][53];
  
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
      for(int j=0; j<53; j++)
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
	for(int j=0; j<53; j++)
	{

	  weightedSelectedEventsWplus[iBin][j] = AllWeightsWplus[iBin*53 + j] ;
	  weightedSelectedEventsWminus[iBin][j] = AllWeightsWminus[iBin*53 + j] ;
	   
          //weightedSelectedEventsWplus[iBin][j] = AllWeightsWplus[j] ;
	  //weightedSelectedEventsWminus[iBin][j] = AllWeightsWminus[j] ;
	  
	  
	  //weightedSelectedEvents[iBin][j] = AllWeightsWplus[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWminus[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWplus[iBin*53 + j] + AllWeightsWminus[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWminus[iBin*53 + j]/AllWeightsWplus[iBin*53 + j] ;
	  
	 
	  //cout<<iBin<<"\t "<<j<<"\t"<<weightedSelectedEvents[iBin][j]<<endl;
	  
	  //cout << Form("weightedSelectedEvents[%d][%d] : %f \t weighted2SelectedEvents[%d][%d] : %f \t weights_CT10[%d] : %f ",iBin,j,weightedSelectedEvents[iBin][j],iBin,j,weighted2SelectedEvents[iBin][j], j, (*weights_CT10)[j]) << endl;
      }
    }
  

  //double weightedSelectedEventsCentTot=0;
  double weightedSelectedEventsWplusTot_central=0;
  double weightedSelectedEventsWminusTot_central=0;
  for(int iBin(0);iBin<nptBins;iBin++)
  {
	  weightedSelectedEventsWplusTot_central  += weightedSelectedEventsWplus[iBin][0] ; 
	  weightedSelectedEventsWminusTot_central += weightedSelectedEventsWminus[iBin][0];
  }
  
  
  // Calculate PDF syst
  for(int iBin(0);iBin<nptBins;iBin++)
  {
    //unsigned int nmembers = weights_CT10->size();
    unsigned int nmembers = 52;
    unsigned int npairs = (nmembers-1)/2;
 
    //events_central[iBin] = weightedSelectedEvents[iBin][0];
    //events_central[iBin] = weightedSelectedEvents[iBin][0]/weightedSelectedEventsCentTot/(xbins_pt[iBin+1]-xbins_pt[iBin]);
    events_central[iBin] = (weightedSelectedEventsWminus[iBin][0]/weightedSelectedEventsWminusTot_central)/(weightedSelectedEventsWplus[iBin][0]/weightedSelectedEventsWplusTot_central);
    if(npairs>0){
      wplus[iBin]=0.;
      wminus[iBin]=0.;

      for (unsigned int j=0; j<npairs; ++j)
      {
          
        double weightedSelectedEventsWplusTot_2j_1=0;
        double weightedSelectedEventsWplusTot_2j_2=0;
        double weightedSelectedEventsWminusTot_2j_1=0;
        double weightedSelectedEventsWminusTot_2j_2=0;
        for(int iBin(0);iBin<nptBins;iBin++)
        {
          weightedSelectedEventsWplusTot_2j_1 += weightedSelectedEventsWplus[iBin][2*j+1];
          weightedSelectedEventsWplusTot_2j_2 += weightedSelectedEventsWplus[iBin][2*j+2];

          weightedSelectedEventsWminusTot_2j_1 += weightedSelectedEventsWminus[iBin][2*j+1];
          weightedSelectedEventsWminusTot_2j_2 += weightedSelectedEventsWminus[iBin][2*j+2];
        }

        //wa[iBin] = weightedSelectedEvents[iBin][2*j+1]/events_central[iBin]-1.;
        //wb[iBin] = weightedSelectedEvents[iBin][2*j+2]/events_central[iBin]-1.;
        //wa[iBin] = (weightedSelectedEvents[iBin][2*j+1]/weightedSelectedEventsXPlusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
        //wb[iBin] = (weightedSelectedEvents[iBin][2*j+2]/weightedSelectedEventsXMinusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;

        wa[iBin] = (
                        (weightedSelectedEventsWminus[iBin][2*j+1]/weightedSelectedEventsWminusTot_2j_1)/(weightedSelectedEventsWplus[iBin][2*j+1]/weightedSelectedEventsWplusTot_2j_1)
                    )/events_central[iBin]-1.;

        wb[iBin] = (
                        (weightedSelectedEventsWminus[iBin][2*j+2]/weightedSelectedEventsWminusTot_2j_2)/(weightedSelectedEventsWplus[iBin][2*j+2]/weightedSelectedEventsWplusTot_2j_2)
                    )/events_central[iBin]-1.;

	//wa[iBin] = weightedSelectedEvents[iBin][2*j+1]/events_central[iBin]-1.;
	//wb[iBin] = weightedSelectedEvents[iBin][2*j+2]/events_central[iBin]-1.;
	//wa[iBin] = (weightedSelectedEvents[iBin][2*j+1]/weightedSelectedEventsXPlusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
	//wb[iBin] = (weightedSelectedEvents[iBin][2*j+2]/weightedSelectedEventsXMinusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
	
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



    double norDiffXsecP[12];

    norDiffXsecP[0] = 0.0470371954 ;
    norDiffXsecP[1] = 0.0416713258 ;
    norDiffXsecP[2] = 0.0249537889 ;
    norDiffXsecP[3] = 0.0127191793 ;
    norDiffXsecP[4] = 0.0059917004 ;
    norDiffXsecP[5] = 0.0034012780 ;
    norDiffXsecP[6] = 0.0016768343 ;
    norDiffXsecP[7] = 0.0005024832 ;
    norDiffXsecP[8] = 0.0001204026 ;
    norDiffXsecP[9] = 0.0000382847 ;
    norDiffXsecP[10]= 0.0000132836  ;
    norDiffXsecP[11]= 0.0000011298  ;

    double norDiffXsecM[12];
    norDiffXsecM[0] = 0.0452816676 ;
    norDiffXsecM[1] = 0.0414175273 ;
    norDiffXsecM[2] = 0.0253264431 ;
    norDiffXsecM[3] = 0.0130558561 ;
    norDiffXsecM[4] = 0.0063198660 ;
    norDiffXsecM[5] = 0.0036153198 ;
    norDiffXsecM[6] = 0.0017662663 ;
    norDiffXsecM[7] = 0.0005253867 ;
    norDiffXsecM[8] = 0.0001249459 ;
    norDiffXsecM[9] = 0.0000412788 ;
    norDiffXsecM[10]= 0.0000129562  ;
    norDiffXsecM[11]= 0.0000010012  ;


    //// normalized PDF sysy in %
    double Norm_PDFErrP[12];
       Norm_PDFErrP[0]      =0.547639 ;
       Norm_PDFErrP[1]      =0.11084  ;
       Norm_PDFErrP[2]      =0.423104 ;
       Norm_PDFErrP[3]      =0.502498 ;
       Norm_PDFErrP[4]      =0.523514 ;
       Norm_PDFErrP[5]      =0.456424 ;
       Norm_PDFErrP[6]      =0.390481 ;
       Norm_PDFErrP[7]      =0.423705 ;
       Norm_PDFErrP[8]      =1.07401  ;
       Norm_PDFErrP[9]      =1.92522  ;
       Norm_PDFErrP[10]     =2.57608  ;
       Norm_PDFErrP[11]     =3.85306  ;
    
    double Norm_PDFErrM[12];
       Norm_PDFErrM[0]      =0.805695 ;
       Norm_PDFErrM[1]      =0.193853 ;
       Norm_PDFErrM[2]      =0.529425 ;
       Norm_PDFErrM[3]      =0.705878 ;
       Norm_PDFErrM[4]      =0.697846 ;
       Norm_PDFErrM[5]      =0.707969 ;
       Norm_PDFErrM[6]      =0.529504 ;
       Norm_PDFErrM[7]      =0.780084 ;
       Norm_PDFErrM[8]      =1.43238  ;
       Norm_PDFErrM[9]      =2.21681  ;
       Norm_PDFErrM[10]     =3.18584  ;
       Norm_PDFErrM[11]     =3.50682  ;
   
    double RatioSyst[12];
    double Ratio[12];
    for (int i=0; i<12;i++)
    {
      Norm_PDFErrP[i]=0.01*Norm_PDFErrP[i]*norDiffXsecP[i];
      Norm_PDFErrM[i]=0.01*Norm_PDFErrM[i]*norDiffXsecM[i];
      Ratio[i]=norDiffXsecM[i]/norDiffXsecP[i];
      RatioSyst[i]=Ratio[i]*sqrt(Norm_PDFErrM[i]*Norm_PDFErrM[i]/norDiffXsecM[i]/norDiffXsecM[i]  +  Norm_PDFErrP[i]*Norm_PDFErrP[i]/norDiffXsecP[i]/norDiffXsecP[i]);
      cout<<"W-/W+  RatioSyst[i] by error propagation formula\t"<< RatioSyst[i] <<"\t \t"<< RatioSyst[i]/(norDiffXsecM[i]/norDiffXsecP[i])*100 <<"\t %"<<endl;
    }

        






}
