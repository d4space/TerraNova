#include <iostream>
#include <set>

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

const int nptBins = 18;
const double xbins_pt[nptBins+1] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 30, 40, 50, 70, 90, 110, 150, 190, 250, 600};

void PDFSystCalculator()
{
  typedef std::pair<double,double> pairOfDouble;
  set<pairOfDouble> uniqueEventsTrue;

  typedef std::pair<int,int> pairOfInt;
  set<pairOfInt> uniqueEventsReco;

  double mLow = 60.; // Z mass
  double mHigh = 120.; // Z mass
  double ptLow = 0.;
  double ptHigh = 4000.;

  //variable related PDF systematic uncertainty
  double weightedSelectedEvents[18][53];
  double weighted2SelectedEvents[18][53];
  double events_central[nptBins] = {0.0};
  double events2_central[nptBins] = {0.0};
  double wa[nptBins] = {0.0};
  double wb[nptBins] = {0.0};
  double wplus[nptBins] = {0.0};
  double wminus[nptBins] = {0.0};


  // read input root file
  TChain* tree = new TChain("tree");
  tree->AddFile("DYJetsToLL_nocut_1.root");  // Madgraph Ntuple
  tree->AddFile("DYJetsToLL_nocut_2.root");  // Madgraph Ntuple
  tree->AddFile("DYJetsToLL_nocut_3.root");  // Madgraph Ntuple

  TH1D* hPtMCTot = new TH1D("hPtMCTot", "", nptBins, xbins_pt); // no cut
  TH1D* hPtMCAcc = new TH1D("hPtMCAcc", "", nptBins, xbins_pt); // acc.cuts
  TH1D* hPtMCEff = new TH1D("hPtMCEff", "", nptBins, xbins_pt); // acc + all selection cuts


  //Hammid MC truth information pre FSR
  _TrackInfo true1PreFSR, true2PreFSR;
  tree->SetBranchAddress("true1PreFSR", &true1PreFSR);
  tree->SetBranchAddress("true2PreFSR", &true2PreFSR);

  float trueMassPreFSR, truePtPreFSR;
  tree->SetBranchAddress("trueMassPreFSR",	&trueMassPreFSR);
  tree->SetBranchAddress("truePtPreFSR",	&truePtPreFSR);

  vector<float> *weights_CT10;
  TBranch *b_weights_CT10;
  weights_CT10 = 0;
  tree->SetBranchAddress("weights_CT10", &weights_CT10, &b_weights_CT10);


  // Initialize variable
  for(int iBin(0); iBin<nptBins; iBin++)
  {
    {
      for(int j=0; j<53; j++)
      {
	weightedSelectedEvents[iBin][j] = 0;
	weighted2SelectedEvents[iBin][j] = 0;
//	cout << "initialized weightedSelectedEvents : " << weightedSelectedEvents[iBin][j] << endl;
      }
    }
  }

  // Print # of events
  cout << "Loop over the " << tree->GetEntries() << " entries ...\n";

  // Fill histogram
  for(int iEvt(0); iEvt<tree->GetEntries(); iEvt++)
  //for(int iEvt(0); iEvt<20; iEvt++)
  {
    if ( (iEvt % 100000)==0 ) cout << "event " << iEvt << endl;
    tree -> GetEntry(iEvt);

    //nocut start
    pairOfDouble trueMassPtPreFSR(trueMassPreFSR,truePtPreFSR);
    if( !uniqueEventsTrue.insert( trueMassPtPreFSR ).second ) continue;

    if (trueMassPreFSR <  mLow) continue;
    if (trueMassPreFSR > mHigh) continue;

    hPtMCTot -> Fill( truePtPreFSR );

    // acc cut start
    if (true1PreFSR.charge == -999) continue;
    if (true2PreFSR.charge == -999) continue;
    if (true1PreFSR.charge*true2PreFSR.charge>0) continue;

    if (trueMassPreFSR <  mLow) continue;
    if (trueMassPreFSR > mHigh) continue;
    if (true1PreFSR.pt < 20         || true2PreFSR.pt < 20        ) continue;
    if (fabs(true1PreFSR.eta) > 2.1 || fabs(true2PreFSR.eta) > 2.1) continue;
    hPtMCAcc -> Fill( truePtPreFSR );

    //cout << "trueMassPreFSR : " << trueMassPreFSR << endl;
    //cout << "true1PreFSR : " << true1PreFSR.pt << endl;
    //cout << "truePtPreFSR : " << truePtPreFSR << endl;
    //cout << "weights size : " << weights_CT10->size() << endl;


    // Check weight numbers are normal. 
    //cout << Form("iEvt : %d ============================= ",iEvt) << endl;
    //for(int i(0); i<weights_CT10->size(); i++)
    //{
    //  cout << "weights_CT10 : " << (*weights_CT10)[i] << endl;
    //}


    // Save weights each bin
    for(int iBin(0); iBin<nptBins; iBin++)
    {
      if(truePtPreFSR >= xbins_pt[iBin] && truePtPreFSR < xbins_pt[iBin+1])
      {
	for(int j=0; j<weights_CT10->size(); j++)
	{
	  weightedSelectedEvents[iBin][j] += (*weights_CT10)[j];
	  weighted2SelectedEvents[iBin][j] += (*weights_CT10)[j] * (*weights_CT10)[j];
	  //cout << Form("weightedSelectedEvents[%d][%d] : %f \t weighted2SelectedEvents[%d][%d] : %f \t weights_CT10[%d] : %f ",iBin,j,weightedSelectedEvents[iBin][j],iBin,j,weighted2SelectedEvents[iBin][j], j, (*weights_CT10)[j]) << endl;
	}
      }
    }
  }

  // Calculate PDF syst
  for(int iBin(0);iBin<nptBins;iBin++)
  {
    unsigned int nmembers = weights_CT10->size();
    unsigned int npairs = (nmembers-1)/2;
    events_central[iBin] = weightedSelectedEvents[iBin][0];
    events2_central[iBin] = weighted2SelectedEvents[iBin][0];
    if(npairs>0){
      for (unsigned int j=0; j<npairs; ++j)
      {
	//cout << "events central : " << events_central[iBin] << endl;
	wa[iBin] = weightedSelectedEvents[iBin][2*j+1]/events_central[iBin]-1.;
	wb[iBin] = weightedSelectedEvents[iBin][2*j+2]/events_central[iBin]-1.;
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
    cout <<iBin+1<<" Bin: Relative uncertainty with respect to central member: +" << 100.*wplus[iBin] << " / -" <<  100.*wminus[iBin] << " [%]" << endl;
  }

  for(int i(0); i<nptBins; i++)
  {
    cout << i+1 << " bin Events number : " << hPtMCAcc->GetBinContent(i+1) << endl;
  }

  TFile* Hist_out = new TFile("preFSR_Madgraph.root","recreate");

  hPtMCTot->Write();
  hPtMCAcc->Write();
}
