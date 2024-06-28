#include<iostream>
#include<TH1D.h>
#include<TF1.h>
#include<TCanvas.h>
#include<TRandom.h>
#include<TStyle.h>
#include<TLegend.h>
#include<TROOT.h>
#include<TPaveText.h>

const double DCT_wire_length=450.0; // same as geometry toml file
const double DCT_wire_resistance=2200.0; // in Ohms
const double MIN_CHARGEDIV=0.0;
const double MAX_CHARGEDIV=1.0;
const double fXPos = 0.0;
double convert_qdv_to_r_n(double in){
  return (1.0-in)*DCT_wire_resistance;
}

double convert_qdv_to_pos(double in){
  double fChargeDivPos = fXPos + DCT_wire_length*((in-MIN_CHARGEDIV)/(MAX_CHARGEDIV-MIN_CHARGEDIV) - 0.5);
  return fChargeDivPos;
}
void make_histo_pretty(TH1D* h,int setblue){
  h->SetLineWidth(2);
  if(setblue==0) h->SetLineColor(kBlue);
  else if(setblue==1) h->SetLineColor(kRed);
  else if(setblue==2) h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle("idk");
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitle("entries");
}
void make_histo_pretty_qdvpos(TH1D* h){
  h->SetLineWidth(2);
  h->SetLineColor(kBlue);
  h->GetXaxis()->SetTitle("South end                                    North end");
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitle("entries");
  h->Draw();
}

void test_random()
{
  double termination_resN=1.5*DCT_wire_resistance;
  double termination_resS=1.1*DCT_wire_resistance;
  double Gain_ratio=0.5;
  double Gain_ratio2=1.3;
  int num_entries=100000;
  gROOT->Reset();
  TStyle * plain = new TStyle("plain","plain");
  plain->SetCanvasBorderMode(0);
  plain->SetPadBorderMode(0);
  plain->SetPadColor(0);
  plain->SetCanvasColor(0);
  plain->SetTitleColor(1);
  plain->SetStatColor(0);
  plain->SetTitleFillColor(0);
  gROOT->SetStyle("plain");
  gStyle->SetPalette(1);

  // min and max chargediv can be used to first sample for fraction
  TH1D * h_div = new TH1D("Charge Division", "", 200,-0.05,1.05);  
  //create histogram for true position
  TH1D * h_true_pos = new TH1D("True Position", "", 200,-1.1*DCT_wire_length/2.0,1.1*DCT_wire_length/2.0);
  TH1D * h_calcqdv = new TH1D("Charge Division calculated", "", 200,-0.05,1.05);  
  TH1D * h_calcqdv2 = new TH1D("Charge Division calculated2", "", 200,-0.05,1.05);  

  //disable display of histogram statistics
  h_div->SetStats(false);
  h_calcqdv->SetStats(false);
  h_true_pos->SetStats(false);
  //fill with true hit positions
  for(double i = 0; i < num_entries; i++){
    double qdv_true; 
    if(i<0.01*num_entries) qdv_true = gRandom->Uniform(0,1);
    else qdv_true=gRandom->Gaus(0.5,0.1);
    //double qdv_true =gRandom->Gaus(0.5,0.1);
    h_div->Fill(qdv_true);
    double qdv_true_pos= convert_qdv_to_pos(qdv_true);
    h_true_pos->Fill(qdv_true_pos);
    //h_true_pos->Fill(gRandom->Uniform(-1.0*DCT_wire_length/2.0,DCT_wire_length/2.0));

    // if qdv is close to 1 then the hit position was close to north side and we should get a larger voltage reading for north side and less resistance from wire contributing!
    double res_north=(1.0-qdv_true)*DCT_wire_resistance;
    double res_south=qdv_true*DCT_wire_resistance;
    double total_resN=res_north+termination_resN;
    double total_resS=res_south+termination_resS;
    //double qdv_new=(total_resS)/(total_resN+total_resS);
    //double qdv_new=res_south/(res_south+(Gain_ratio*res_north));
    double qdv_new=total_resS/(total_resS+(Gain_ratio2*total_resN));
    h_calcqdv->Fill(qdv_new);
    //qdv_new=res_south/(res_south+(Gain_ratio2*res_north));
    //h_calcqdv2->Fill(qdv_new);
    // now we can calculate the 
    //h_true_pos->Fill(gRandom->Uniform(-1.0*DCT_wire_length/2.0,DCT_wire_length/2.0));
  }
  // now apply a voltage to this


/*
  //define fit functions
  TF1 * FitFunc1 = new TF1("FitFunc1","[0]*TMath::Gaus(x,[1],[2])",0,40);
  TF1 * FitFunc2 = new TF1("FitFunc2","[0]*TMath::Landau(x,[1],[2])",0,40);
  TF1 * FitFuncCombined = new TF1("FitFunc2","[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Landau(x,[4],[5])",0,40);
  //fit both peaks individually with reasonable initial parameters and fitting range
  FitFunc1->SetParameters(1,3,4);
  h->Fit(FitFunc1,"0","",0,10);
  FitFunc2->SetParameters(1,17,7);
  h->Fit(FitFunc2,"0","",10,40);
  //use fit parameters as initial parameters for combined fit
  FitFuncCombined->SetParameters(FitFunc1->GetParameter(0), FitFunc1->GetParameter(1), FitFunc1->GetParameter(2), FitFunc2->GetParameter(0), FitFunc2->GetParameter(1), FitFunc2->GetParameter(2));
  h->Fit(FitFuncCombined,"0","");
  //display what we did
*/
  TCanvas * c = new TCanvas("c_ref","c_title", 200,10,600,600);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.04);
  c->SetTopMargin(0.04);
  //Legend
  TLegend* Leg = new TLegend(0.3,0.8,0.99,0.99);
  Leg->SetFillColor(0);
  Leg->SetTextFont(62);
  make_histo_pretty(h_div,0);
  make_histo_pretty(h_calcqdv,1);
  //make_histo_pretty(h_calcqdv2,2);

  h_calcqdv->Draw();
  //h_calcqdv2->Draw("SAME");
  h_div->Draw("SAME");
  c->SaveAs("GainRatio_resistance_model.png");

  TCanvas * c2 = new TCanvas("c_ref2","c_title2", 200,50,600,600);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.04);
  c2->SetTopMargin(0.04);
  make_histo_pretty_qdvpos(h_true_pos);

  /*char text[400];
  sprintf(text,"N=%5.0f Mean=%5.1f RMS=%5.1f", h->GetEntries(), h->GetMean(), h->GetRMS());
  Leg->AddEntry(h,text,"l");
  FitFunc1->SetLineStyle(2);
  FitFunc1->SetLineColor(kRed);
  FitFunc1->Draw("same");
  sprintf(text,"Gaus: Mean=%5.1f#pm%5.1f, #sigma=%5.1f#pm%5.1f Landau: MOP=%5.1f#pm%5.1f, #sigma=%5.1f#pm%5.1f", FitFuncCombined->GetParameter(1), 
  FitFuncCombined->GetParError(1), FitFuncCombined->GetParameter(2), FitFuncCombined->GetParError(2), FitFuncCombined->GetParameter(4), FitFuncCombined->GetParError(4), FitFuncCombined->GetParameter(5), FitFuncCombined->GetParError(5));
  Leg->AddEntry(FitFuncCombined,text,"l");
  FitFunc2->SetLineStyle(2);
  FitFunc2->SetLineColor(kRed);
  FitFunc2->Draw("same");
  FitFuncCombined->SetLineColor(kRed);
  FitFuncCombined->Draw("same");
  Leg->Draw();
  //Save canvas
  c->SaveAs("ex1.eps");
  c->SaveAs("ex1.png");
  c->SaveAs("ex1.root");*/
}
