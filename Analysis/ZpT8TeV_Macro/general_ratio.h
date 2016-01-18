#include<iostream>
#include<TLatex.h>
#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<THStack.h>
#include<TCanvas.h>
#include<TLegend.h>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>

//#include "/afs/cern.ch/user/d/digiovan/scripts/Init/modifiedStyle.C"

// ---------------------------------------------------------------------------
// general variables
double const PDG_MASS_Z  = 91.1876;//GeV/c2
double const PDG_WIDTH_Z = 2.4952; //GeV/c2
//const double Pi = TMath::Pi();

// useful strings
char cmspreliminary[128];
char cmspreliminaryextd[128];
char cmssimulation[128];

char lumiString[128];
char energyString[128];

TLatex *lumi = new TLatex();
TLatex *energy = new TLatex();

void loadGeneral() {
  sprintf(cmspreliminary, "CMS Preliminary");
  sprintf(cmspreliminaryextd, "CMS Preliminary, %3.1f pb^{-1} at #sqrt{s}=8 TeV", 18.424);
  sprintf(cmssimulation, "CMS Simulation");

  // intLumi is in include/datasets.h
  //sprintf(lumiString, "#int L=%4.2f pb^{-1}", intLumi_invpb);
  //sprintf(lumiString, "#int L=%3.1f pb^{-1}", 18.424);
  sprintf(lumiString, "L=%3.1f pb^{-1}", 18.424);
  //sprintf(lumiString, "18.4 pb^{-1}");

  lumi->SetTextFont (  42);
  lumi->SetTextSize(0.040);
  lumi->SetTextAlign  (22);

  energy->SetTextFont (  42);
  energy->SetTextSize(0.040);
  energy->SetTextAlign  (22);

  sprintf(energyString, "#sqrt{s} = 8 TeV");
}


void PrintIt(TPad *pad, TString title)//char *title)
{
  TLatex *latex = new TLatex();
  latex->SetTextFont(  42);
  //latex->SetTextSize(0.05);
  latex->SetTextSize(0.07);

  // Get the most recent changes
  pad->Update();


  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();

  //double xpos = xmin + 0.50*(xmax - xmin);
  //double ypos = ymax + 0.05*(ymax - ymin);
  //double xpos = xmin + 0.005*(xmax - xmin);
  double xpos = xmin;
  double ypos = ymax + 0.02*(ymax - ymin);

  if (pad->GetLogy())  ypos = pow(10,ypos);
  if (pad->GetLogx())  xpos = pow(10,xpos);

  //latex->SetTextAlign(22);
  latex->DrawLatex(xpos,ypos,title);
}

void PrintIt2(TPad *pad, TString title)//char *title)
{
  TLatex *latex = new TLatex();
  latex->SetTextFont(  42);
  //latex->SetTextSize(0.05);
  latex->SetTextSize(0.07);

  // Get the most recent changes
  pad->Update();


  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();

  double xpos = xmax - 0.38*(xmax - xmin);
  double ypos = ymax + 0.02*(ymax - ymin);

  if (pad->GetLogy())  ypos = pow(10,ypos);
  if (pad->GetLogx())  xpos = pow(10,xpos);

  //latex->SetTextAlign(22);
  latex->DrawLatex(xpos,ypos,title);
}

void PrintIt3(TPad *pad, TString title)//char *title)
{
  TLatex *latex = new TLatex();
  latex->SetTextFont(  42);
  //latex->SetTextSize(0.05);
  latex->SetTextSize(0.07);

  // Get the most recent changes
  pad->Update();


  double xmin = pad->GetUxmin();
  double xmax = pad->GetUxmax();
  double ymin = pad->GetUymin();
  double ymax = pad->GetUymax();

  double xpos = xmax - 0.58*(xmax - xmin);
  double ypos = ymax - 0.12*(ymax - ymin);

  if (pad->GetLogy())  ypos = pow(10,ypos);
  if (pad->GetLogx())  xpos = pow(10,xpos);

  //latex->SetTextAlign(22);
  latex->DrawLatex(xpos,ypos,title);
}
//------------------------------------------------------------------------------
// Draw projections and residuals
//------------------------------------------------------------------------------
void axis1F(TH1F  *histo,
            TAxis *xaxis,
            TAxis *yaxis,
            char  *xtitle,
            char  *ytitle)
{
  histo->SetMarkerSize(0.5);
  histo->SetMarkerStyle(kFullCircle);

  histo->SetTitle("");

  xaxis = histo->GetXaxis();
  yaxis = histo->GetYaxis();

  xaxis->SetLabelFont(42);
  yaxis->SetLabelFont(42);
  xaxis->SetLabelOffset(0.005);
  yaxis->SetLabelOffset(0.005);
  //xaxis->SetLabelSize(0.04);
  yaxis->SetLabelSize(0.04);

  xaxis->SetNdivisions(505);
  yaxis->SetNdivisions(505);

  xaxis->SetTitle(xtitle);
  yaxis->SetTitle(ytitle);
  xaxis->SetTitleColor(kBlack);
  yaxis->SetTitleColor(kBlack);
  xaxis->SetTitleFont(42);
  yaxis->SetTitleFont(42);
  xaxis->SetTitleOffset(1.0);
  yaxis->SetTitleOffset(1.3);
  //xaxis->SetTitleSize(0.045);
  //yaxis->SetTitleSize(0.045);
  yaxis->CenterTitle(kTRUE);
}

void DrawWithRes(TCanvas* canvas, char* title, char* title2, 
                 TH1F *hData, TH1F *hMC, THStack *sMC, 
                 TLegend *tl=NULL, bool isLogScaleY=false, bool isLogScaleX=false) // From bin1
                 //TLegend *tl=NULL, bool isLogScaleY=true, bool isLogScaleX=false) // After bin9
{

  // sanity check
  if (hData->GetNbinsX() != hMC->GetNbinsX()){
    std::cout<<" *** Error: binning not consistent between data and MC -> Exit!\n";
    return;
  }


  TH1F *hPull = (TH1F*)hData ->Clone("hPull");
  hPull->Reset();

  for (int i=1; i<=hPull->GetNbinsX(); i++) {
    double data  = hData->GetBinContent(i);
    double mc    = hMC->GetBinContent(i);
  
 Double_t error_highM, error_lowM, error_highD, error_lowD; 
   double sigmaM, sigmaD;
     error_highM = hMC->GetBinErrorUp(i);
     error_lowM  = hMC->GetBinErrorLow(i);
     if (error_highM > error_lowM) sigmaM = error_highM;
     else sigmaM = error_lowM;

     error_highD = hData->GetBinErrorUp(i);
     error_lowD  = hData->GetBinErrorLow(i);
     if (error_highD > error_lowD) sigmaD = error_highD;
     else sigmaD = error_lowD;

    double error = hData->GetBinError(i);
    double resid = (error) ? ((data-mc)/error) : 0.0;
    if (fabs(resid) > 100) resid = 10;
  
    hPull->SetBinContent(i,resid); // Chi2 style
    //double ratio = mc/data; // MC/Data style

     double sigmaR = sqrt(((1/data)*(1/data))*(sigmaD*sigmaD)+ ((data/(mc*mc))*(data/(mc*mc)))*(sigmaM*sigmaM));
   //  std::cout << iBin << ") " << yNum << "/" << yDen << ", Ratio= " << ratio << std::endl;
     hPull->SetBinContent(i,resid); // Chi2 style 
     hPull->SetBinError(i,sigmaR);  // Chi2 styl
     //hPull->SetBinContent(i,ratio); // MC/Data style
}

  //----------------------------------------------------------------------------
  // Create the pads
  //----------------------------------------------------------------------------
  TPad* pad1; 
  TPad* pad2; 
  
  pad1 = new TPad("pad1","This is pad1",0.02,0.30,0.99,0.98,0); 
  pad2 = new TPad("pad2","This is pad2",0.02,0.01,0.99,0.29,0); 
  
  if (isLogScaleY) pad1->SetLogy();
  if (isLogScaleX) pad1->SetLogx();
  pad1->SetBottomMargin(0.01);
  pad1->SetLeftMargin(0.12);
  pad2->SetTopMargin   (0.10); // Chi2 style
  //pad2->SetTopMargin   (0.05); // MC/Data style
  pad2->SetBottomMargin(0.33);
  pad2->SetLeftMargin(0.12);

  pad1->SetTicky(1);
  pad1->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetTickx(1);

  pad1->Draw(); // Projections pad
  pad2->Draw(); // Residuals   pad

  pad1->cd();
  //hMC  ->Draw();
  hData->Draw("pe");
  sMC  ->Draw("same");
  hData->Draw("pe same");
  hData->Draw("AXIS same");

  lumi->SetNDC();
  energy->SetNDC();
  lumi  ->SetTextSize(0.060);
  energy->SetTextSize(0.060);
   
/*   TString name(hMC->GetName()); */
/*   if (name.Contains("Mass")) { */
/*     lumi  ->DrawLatex(0.35,0.80, lumiString); */
/*     energy->DrawLatex(0.35,0.70, energyString); */
/*   }    */

/*   if (name.Contains("Pt")) { */
/*     lumi  ->DrawLatex(0.57,0.80, lumiString); */
/*     energy->DrawLatex(0.57,0.68, energyString); */
/*   }    */

  PrintIt(pad1, title);
  PrintIt2(pad1, title2);
  PrintIt3(pad1, "p_{T}^{Z} < 30 GeV"); // From bin1
  //PrintIt3(pad1, "p_{T}^{Z} #geq 30 GeV"); // After bin9

  if (tl) tl ->Draw("same");

  // // Bin range text
  //TPaveText *bintext = new TPaveText(0.7,0.9,0.96,0.95);
  //bintext->AddBox(0.7,0.9,0.96,0.95);
  //bintext->SetFillStyle(2);
  //bintext->SetBorderSize(2);
  //bintext->AddText("p_{T}^{Z} < 30 GeV");
  //bintext->Draw("same");

  
  //----------------------------------------------------------------------------
  // Residuals pad
  //----------------------------------------------------------------------------
  pad2->cd();
  if (isLogScaleX) pad2->SetLogx();

  TH1F* hData_Clone = (TH1F*)hData->Clone("hData_Clone");
  for (int i=1; i<=hPull->GetNbinsX(); i++)
  {
    hData_Clone->SetBinContent(i,1);
    hData_Clone->SetBinError(i,hData->GetBinErrorUp(i) / hData->GetBinContent(i));
    cout << Form("Error : %.8f ", hData_Clone->GetBinErrorUp(i)) << endl; 
  }

  TAxis *xPull = NULL;
  TAxis *yPull = NULL;
  char xAxisName[200];
  sprintf(xAxisName,"%s",hData->GetXaxis()->GetTitle());
  axis1F(hPull,xPull,yPull,xAxisName,"(Data-MC)/#sigma_{Data}"); // Chi2 style
  //axis1F(hPull,xPull,yPull,xAxisName,"MC/Data"); // MC/Data style 
  //if (hPull->GetMaximum() > 100) {
  //std::cout << "Setting the hPull boundaries!\n";
  hPull->SetMinimum(-5); // Chi2 style 
  hPull->SetMaximum(5);  // Chi2 style
  //hPull->SetMinimum(0.6); // MC/Data style , From bin1  
  //hPull->SetMaximum(1.2); // MC/Data style , From bin1
  //hPull->SetMinimum(-1); //  MC/Data style , After bin9
  //hPull->SetMaximum(3); //   MC/Data style , After bin9
  //}


  hPull->GetXaxis()->SetLabelOffset(0.005);
  hPull->GetXaxis()->SetLabelSize  (0.11);
  //hPull->GetXaxis()->CenterTitle(1);
  hPull->GetXaxis()->SetTitleOffset(1.0);
  hPull->GetXaxis()->SetTitleSize  (0.12);
  hPull->GetXaxis()->SetNdivisions(510);

  hPull->GetYaxis()->SetLabelSize(0.12);
  hPull->GetYaxis()->SetNdivisions(104);
  hPull->GetYaxis()->CenterTitle(1);
  //hPull->GetYaxis()->SetTitleOffset(0.3);
  hPull->GetYaxis()->SetTitleOffset(0.4);
  hPull->GetYaxis()->SetTitleSize  (0.12);

  hPull->SetFillColor(856);
  hPull->SetFillStyle(1001);
  hPull->Draw("histo b"); //
  //hPull->SetMarkerColor(2);  // MC/Data style  
  //hPull->SetMarkerStyle(20); // MC/Data style
  //hPull->SetMarkerSize(1.0); // MC/Data style
 
  //TColor *colBlue = gROOT->GetColor(kBlue+1);
  //colBlue->SetAlpha(0.3);

  TGraphAsymmErrors* tgDataErrBand = new TGraphAsymmErrors(hData_Clone);
  tgDataErrBand->SetFillStyle(3001); 
  tgDataErrBand->SetFillStyle(0); 
  tgDataErrBand->SetFillColor(kBlack); 
  //tgDataErrBand->SetFillColor(kBlue+1); // MC/Data style 
  //tgDataErrBand->SetLineColor(kBlack);  // MC/Data style

  //hPull->Draw("P");
  //tgDataErrBand->Draw("2"); 

  //pad2->Update();
  //pad2->GetFrame()->DrawClone();

}


void outLatex(TH1F* histo, 
              TString name) {

    std::cout << name << " & $" << histo->Integral() << " \\pm " 
            << TMath::Sqrt(histo->Integral())      << "$\\\\\n";

    std::cout << "\\hline\n";
}
