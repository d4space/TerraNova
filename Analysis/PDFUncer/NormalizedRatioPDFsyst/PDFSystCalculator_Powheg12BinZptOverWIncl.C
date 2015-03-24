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

void PDFSystCalculator_Powheg12BinZptOverWIncl()
{

  cout << fixed << setprecision(5);
  //variable related PDF systematic uncertainty
  double weightedSelectedEvents[12][53];
  double weighted2SelectedEvents[12][53];
  
  double weightedSelectedEventsWplus[12][53];
  double weightedSelectedEventsWminus[12][53];
  double weightedSelectedEventsWIncl[12][53];
  double weightedSelectedEventsZpt[12][53];
  
  double weightedSelectedEventsWplusNorm[12][53];
  double weightedSelectedEventsWminusNorm[12][53];
  double weightedSelectedEventsWInclNorm[12][53];
  double weightedSelectedEventsZptNorm[12][53];
  
  double events_central[nptBins] = {0.0};
  double events2_central[nptBins] = {0.0};
  double wa[nptBins] = {0.0};
  double wb[nptBins] = {0.0};
  double wplus[nptBins] = {0.0};
  double wminus[nptBins] = {0.0};

  char tmpName[30];


  ifstream Wplus_file("InputWeights/WpToMuNu_12BinWeightedSelectedEvents.txt");
  ifstream Wminus_file("InputWeights/WmToMuNu_12BinWeightedSelectedEvents.txt");
  ifstream Zpt_file("InputWeights/Zpt_12BinWeightedSelectedEvents.txt");

  // Initialize variable
  for(int iBin(0); iBin<nptBins; iBin++)
  {
    {
      for(int j=0; j<53; j++)
      {
	weightedSelectedEvents[iBin][j] = 0;
	weightedSelectedEventsWplus[iBin][j] = 0;
	weightedSelectedEventsWminus[iBin][j] = 0;
	weightedSelectedEventsWplusNorm[iBin][j] = 0;
	weightedSelectedEventsWminusNorm[iBin][j] = 0;

	weightedSelectedEventsWIncl[iBin][j] = 0;
	weightedSelectedEventsWInclNorm[iBin][j] = 0;
	weightedSelectedEventsZpt[iBin][j] = 0;
	weightedSelectedEventsZptNorm[iBin][j] = 0;
//	cout << "initialized weightedSelectedEvents : " << weightedSelectedEvents[iBin][j] << endl;
      }
    }
  }


  TString tmp;
  int bins;
  int Npdf;
  

  double AllWeightsWplus[636];
  double AllWeightsWminus[636];
  double AllWeightsZpt[636];

  int pdf=0;
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
   pdf=0;
   while(Zpt_file >>  AllWeightsZpt[pdf] )
      {
         //cout << AllWeightsZpt[pdf] << endl;
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
	  
	   
	  weightedSelectedEventsWIncl[iBin][j] = AllWeightsWplus[iBin*53 + j] + AllWeightsWminus[iBin*53 + j] ;
	  
	  weightedSelectedEventsZpt[iBin][j] = AllWeightsZpt[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWplus[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWminus[iBin*53 + j] ;
	  //weightedSelectedEvents[iBin][j] = AllWeightsWminus[iBin*53 + j]/AllWeightsWplus[iBin*53 + j] ;
	  
	 
	  //cout<<iBin<<"\t "<<j<<"\t"<<weightedSelectedEvents[iBin][j]<<endl;
	  cout<<"Zpt\t"<<iBin<<"\t "<<j<<"\t"<<weightedSelectedEventsZpt[iBin][j]<<"\t Wincl\t"<<weightedSelectedEventsWIncl[iBin][j]<<endl;
	  
	  //cout << Form("weightedSelectedEvents[%d][%d] : %f \t weighted2SelectedEvents[%d][%d] : %f \t weights_CT10[%d] : %f ",iBin,j,weightedSelectedEvents[iBin][j],iBin,j,weighted2SelectedEvents[iBin][j], j, (*weights_CT10)[j]) << endl;
      }
    }
  

  //double weightedSelectedEventsCentTot=0;
  double weightedSelectedEventsWInclTot_central=0;
  double weightedSelectedEventsZptTot_central=0;
  for(int iBin(0);iBin<nptBins;iBin++)
  {
	weightedSelectedEventsWInclTot_central  += weightedSelectedEventsWIncl[iBin][0] ; 
	weightedSelectedEventsZptTot_central  += weightedSelectedEventsZpt[iBin][0] ; 
	  cout<<iBin<<"\tweightedSelectedEventsWInclTot_central\t"<<weightedSelectedEventsWInclTot_central<<endl;
  }
  
  
  // Calculate PDF syst
  for(int iBin(0);iBin<nptBins;iBin++)
  {
    //unsigned int nmembers = weights_CT10->size();
    unsigned int nmembers = 52;
    unsigned int npairs = (nmembers-1)/2;
 
    //events_central[iBin] = weightedSelectedEvents[iBin][0];
    events_central[iBin] = (weightedSelectedEventsZpt[iBin][0]/weightedSelectedEventsZptTot_central)/(weightedSelectedEventsWIncl[iBin][0]/weightedSelectedEventsWInclTot_central);
    if(npairs>0){
     wplus[iBin]=0.;
     wminus[iBin]=0.;
      for (unsigned int j=0; j<npairs; ++j)
      {
        double weightedSelectedEventsWInclTot_2j_1=0;
        double weightedSelectedEventsWInclTot_2j_2=0;
        double weightedSelectedEventsZptTot_2j_1=0;
        double weightedSelectedEventsZptTot_2j_2=0;
        for(int iBin(0);iBin<nptBins;iBin++)
        {
          weightedSelectedEventsWInclTot_2j_1 += weightedSelectedEventsWIncl[iBin][2*j+1];
          weightedSelectedEventsWInclTot_2j_2 += weightedSelectedEventsWIncl[iBin][2*j+2];

          weightedSelectedEventsZptTot_2j_1 += weightedSelectedEventsZpt[iBin][2*j+1];
          weightedSelectedEventsZptTot_2j_2 += weightedSelectedEventsZpt[iBin][2*j+2];
        }

	//wa[iBin] = weightedSelectedEvents[iBin][2*j+1]/events_central[iBin]-1.;
	//wb[iBin] = weightedSelectedEvents[iBin][2*j+2]/events_central[iBin]-1.;
	//wa[iBin] = (weightedSelectedEvents[iBin][2*j+1]/weightedSelectedEventsXPlusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
	//wb[iBin] = (weightedSelectedEvents[iBin][2*j+2]/weightedSelectedEventsXMinusTot/(xbins_pt[iBin+1]-xbins_pt[iBin]))/events_central[iBin]-1.;
	
        wa[iBin] = (
			(weightedSelectedEventsZpt[iBin][2*j+1]/weightedSelectedEventsZptTot_2j_1)/(weightedSelectedEventsWIncl[iBin][2*j+1]/weightedSelectedEventsWInclTot_2j_1)
		    )/events_central[iBin]-1.;
	
        wb[iBin] = (
			(weightedSelectedEventsZpt[iBin][2*j+2]/weightedSelectedEventsZptTot_2j_2)/(weightedSelectedEventsWIncl[iBin][2*j+2]/weightedSelectedEventsWInclTot_2j_2)
		    )/events_central[iBin]-1.;


	
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



    double norDiffXsecWincl[12];

    norDiffXsecWincl[0] =0.0463305056  ;
    norDiffXsecWincl[1] =0.0415691589  ;
    norDiffXsecWincl[2] =0.0251038013  ;
    norDiffXsecWincl[3] =0.0128547090  ;
    norDiffXsecWincl[4] =0.0061238039  ;
    norDiffXsecWincl[5] =0.0034874407  ;
    norDiffXsecWincl[6] =0.0017128352  ;
    norDiffXsecWincl[7] =0.0005117030  ;
    norDiffXsecWincl[8] =0.0001222315  ;
    norDiffXsecWincl[9] =0.0000394900  ;
    norDiffXsecWincl[10]=0.0000131518   ;
    norDiffXsecWincl[11]=0.0000010780   ;

    double norDiffXsecZpt[12];
    norDiffXsecZpt[0] =0.0441183071  ;
    norDiffXsecZpt[1] =0.0405764549  ;
    norDiffXsecZpt[2] =0.0253380971  ;
    norDiffXsecZpt[3] =0.0131232996  ;
    norDiffXsecZpt[4] =0.0064591970  ;
    norDiffXsecZpt[5] =0.0037984680  ;
    norDiffXsecZpt[6] =0.0019372698  ;
    norDiffXsecZpt[7] =0.0006126328  ;
    norDiffXsecZpt[8] =0.0001503083  ;
    norDiffXsecZpt[9] =0.0000509369  ;
    norDiffXsecZpt[10]=0.0000184618   ;
    norDiffXsecZpt[11]=0.0000014680   ;
  
    //// normalized PDF sysy in %
    double Norm_PDFErrWincl[12];
       Norm_PDFErrWincl[0]      =0.71170 ;
       Norm_PDFErrWincl[1]      =0.16135 ;
       Norm_PDFErrWincl[2]      =0.49645 ;
       Norm_PDFErrWincl[3]      =0.63787 ;
       Norm_PDFErrWincl[4]      =0.63341 ;
       Norm_PDFErrWincl[5]      =0.58773 ;
       Norm_PDFErrWincl[6]      =0.42823 ;
       Norm_PDFErrWincl[7]      =0.56209 ;
       Norm_PDFErrWincl[8]      =1.26718 ;
       Norm_PDFErrWincl[9]      =2.07182 ;
       Norm_PDFErrWincl[10]     =2.97034 ;
       Norm_PDFErrWincl[11]     =3.63010 ;
    
    double Norm_PDFErrZpt[12];
       Norm_PDFErrZpt[0]      =0.280629 ;
       Norm_PDFErrZpt[1]      =0.205053 ;
       Norm_PDFErrZpt[2]      =0.272395 ;
       Norm_PDFErrZpt[3]      =0.297853 ;
       Norm_PDFErrZpt[4]      =0.315051 ;
       Norm_PDFErrZpt[5]      =0.356816 ;
       Norm_PDFErrZpt[6]      =0.619421 ;
       Norm_PDFErrZpt[7]      =1.289640 ;
       Norm_PDFErrZpt[8]      =2.291440 ;
       Norm_PDFErrZpt[9]      =2.923860 ;
       Norm_PDFErrZpt[10]     =4.046270 ;
       Norm_PDFErrZpt[11]     =6.014940 ;
   
    double RatioSyst[12];
    double Ratio[12];
    for (int i=0; i<12;i++)
    {
      Norm_PDFErrWincl[i]=0.01*Norm_PDFErrWincl[i]*norDiffXsecWincl[i];
      Norm_PDFErrZpt[i]=0.01*Norm_PDFErrZpt[i]*norDiffXsecZpt[i];
      Ratio[i]=norDiffXsecZpt[i]/norDiffXsecWincl[i];
      RatioSyst[i]=Ratio[i]*sqrt(Norm_PDFErrZpt[i]*Norm_PDFErrZpt[i]/norDiffXsecZpt[i]/norDiffXsecZpt[i]  +  Norm_PDFErrWincl[i]*Norm_PDFErrWincl[i]/norDiffXsecWincl[i]/norDiffXsecWincl[i]);
      cout<<"Z/WIncl  RatioSyst[i] by error propagation formula\t"<< RatioSyst[i] <<"\t \t"<< RatioSyst[i]/(norDiffXsecZpt[i]/norDiffXsecWincl[i])*100 <<"\t %"<<endl;
    }

        






}
