{
#include "../../../Utils/const.h"
  // Make Incl
  //*
  TFile *fp = new TFile("../../ResultWpToEleNu/Result_WpToEleNu.root");
  TFile *fm = new TFile("../../ResultWmToEleNu/Result_WmToEleNu.root");

  TH1D *h_data_p;
  TH1D *h_data_m;
  TH1D *h_MC_p;
  TH1D *h_MC_m;
  TH1D *h_MC_p_Unweighted;
  TH1D *h_MC_m_Unweighted;
  TH1D *h_dataRec_p;
  TH1D *h_dataRec_m;

  h_data_p = (TH1D*)fp->Get("BornEffCorr")->Clone("h_data_p");
  h_data_m = (TH1D*)fm->Get("BornEffCorr")->Clone("h_data_m");
  h_MC_p = (TH1D*)fp->Get("SVD_Born.Gen")->Clone("h_MC_p");
  h_MC_m = (TH1D*)fm->Get("SVD_Born.Gen")->Clone("h_MC_m");
  h_MC_p_Unweighted = (TH1D*)fp->Get("SVD_Born.Gen")->Clone("h_MC_p_Unweighted");
  h_MC_m_Unweighted = (TH1D*)fm->Get("SVD_Born.Gen")->Clone("h_MC_m_Unweighted");

  // To make SVD_Born.Gen to h1_Born_BornFid(Not LumiWeighted)
  h_MC_p_Unweighted->Scale(1./LumiWeight_Ele_WpToEleNu_S8);
  h_MC_m_Unweighted->Scale(1./LumiWeight_Ele_WmToEleNu_S8);

  h_dataRec_p = (TH1D*)fp->Get("data_Rec")->Clone("h_dataRec_p");
  h_dataRec_m = (TH1D*)fm->Get("data_Rec")->Clone("h_dataRec_m");
  cout << "Inclusive Cross-section" << endl;
  cout << "bin\tW+\t\tW-\t\t Wincl \t\t PowhegIncl(Unweighted)"<<endl;

  for( int ipt(1);ipt<14;ipt++)
  {
    cout<<ipt<<"\t"<<h_data_p->GetBinContent(ipt)<<"\t\t"<<h_data_m->GetBinContent(ipt)<< "\t\t" << h_data_p->GetBinContent(ipt)+h_data_m->GetBinContent(ipt) << "\t\t" << h_MC_p_Unweighted->GetBinContent(ipt) + h_MC_m_Unweighted->GetBinContent(ipt) << endl;
  }

  TString resultDir = "ResultWInclToEleNu";
  gSystem->mkdir(resultDir,kTRUE);

  TFile f_out(resultDir+"/Result_WInclToEleNu.root","recreate");
  TH1D* BornEffCorr = new TH1D("BornEffCorr","BornEffCorr",13,0,13);
  TH1D* PowhegErr = new TH1D("PowhegErr","PowhegErr",13,0,13);
  TH1D* SVD_BornGen = new TH1D("SVD_Born.Gen","SVD_Born.Gen",13,0,13);
  TH1D* SVD_BornGen_Unweighted = new TH1D("SVD_Born.Gen_Unweighted","SVD_Born.Gen_Unweighted",13,0,13);
  TH1D* data_Rec = new TH1D("data_Rec","data_Rec",13,0,13);

  //h_MC_p->Scale(1./18.429);
  //h_MC_m->Scale(1./18.429);

  for(int i(1);i<14;i++)
  {
    BornEffCorr->SetBinContent(i,h_data_p->GetBinContent(i) + h_data_m->GetBinContent(i));
    SVD_BornGen->SetBinContent(i,h_MC_p->GetBinContent(i) + h_MC_m->GetBinContent(i));
    SVD_BornGen_Unweighted->SetBinContent(i,h_MC_p_Unweighted->GetBinContent(i) + h_MC_m_Unweighted->GetBinContent(i));
    data_Rec->SetBinContent(i,h_dataRec_p->GetBinContent(i) + h_dataRec_m->GetBinContent(i));
  }
  BornEffCorr->Write();
  SVD_BornGen->Write();
  SVD_BornGen_Unweighted->Write();
  data_Rec->Write();

  // */
}
