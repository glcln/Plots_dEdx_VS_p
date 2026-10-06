#include <iostream>
#include <TFile.h>
#include <TH2F.h>


void Return_K_C_params(double &K, double &C, std::string filename)
{
  std::ifstream infile(filename);
  infile >> K >> C;
  infile.close();

  return;
}

void Return_AtlasParams(double &p1, double &p2, double &p3, double &p4, double &p5, std::string filename)
{
  std::ifstream infile(filename);
  infile >> p1 >> p2 >> p3 >> p4 >> p5;
  infile.close();

  return;
}

TF1 *Fit_MassParametrization(double mass, double p_start, double p_end, double C, double K)
{
  TF1 *fit = new TF1("fit", Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  fit->SetParameter(0, K);
  fit->FixParameter(1, C);

  return fit;
}

Double_t AtlasFunction(Double_t *x, Double_t *par)
{
  Double_t m = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t p3 = par[3];
  Double_t p4 = par[4];
  Double_t p5 = par[5];

  Float_t p = x[0];

  Float_t term1 = pow( ( sqrt(pow(p/m,4) + 4*(p/m)*(p/m)) - (p/m)*(p/m) )/2 , p2/2);

  return p1 * term1 * log(1 + pow(p3 * p/m, p4)) - p5;
}

TF1 *Fit_MassParametrization(double mass, double p_start, double p_end, double p1, double p2, double p3, double p4, double p5)
{
  TF1 *fit = new TF1("fit", AtlasFunction, p_start, p_end, 6);
  fit->FixParameter(0, mass);
  fit->FixParameter(1, p1);
  fit->FixParameter(2, p2);
  fit->FixParameter(3, p3);
  fit->FixParameter(4, p4);
  fit->FixParameter(5, p5);
  
  return fit;
}

void Display_TH2_Fit(TH2F *h2, TF1* fit_original_pion, TF1* fit_original_proton, TF1 *fit_pion, TF1 *fit_kaon, TF1 *fit_proton, TF1 *fit_deuteron, std::string filename_PDF, std::string filename_C, std::string filename_ROOT, std::string filename_PNG)
{
  gROOT->SetBatch(kTRUE);

  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("Track momentum [GeV/c]");
  h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  h2->GetZaxis()->SetTitle("Number of tracks");
  h2->SetTitle(" ");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetRangeUser(0, 14);

  TCanvas *c = new TCanvas("c","c",2500,1773);
  c->SetRightMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  h2->Draw("colz");
  c->SetLogz();
  fit_pion->Draw("same");
  fit_kaon->Draw("same");
  fit_proton->Draw("same");
  fit_deuteron->Draw("same");
  fit_original_pion->SetLineColor(kBlack);
  fit_original_proton->SetLineColor(kBlack);
  fit_original_pion->Draw("same");
  fit_original_proton->Draw("same");

  TLatex * tex = new TLatex(0.85,0.915,"2025 (13.6 TeV)");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.10,0.915,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.08);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.22,0.96,"Preliminary");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(52);
  tex->SetTextSize(0.0608);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.12,0.3,"#pi");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.12,0.86,"K");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.16,0.86,"p");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.23,0.86,"D");
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  TLegend *leg = new TLegend(0.57,0.7,0.82,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(fit_original_pion,"Fit to reference data","l");
  leg->AddEntry(fit_kaon,"Extrapolation","l");
  leg->Draw();


  c->SaveAs(filename_PDF.c_str());
  c->SaveAs(filename_C.c_str());
  c->SaveAs(filename_ROOT.c_str());
  c->SaveAs(filename_PNG.c_str());

  delete c;

  return;
}

void Display_TH2_NoFit(TH2F *h2, std::string filename_PDF, std::string filename_C, std::string filename_ROOT, std::string filename_PNG, bool charge)
{
  gROOT->SetBatch(kTRUE);

  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("Track momentum [GeV/c]");
  if (charge) h2->GetXaxis()->SetTitle("Charge sign x track momentum [GeV/c]");
  h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  h2->GetZaxis()->SetTitle("Number of tracks");
  h2->SetTitle(" ");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetRangeUser(0, 14);

  TCanvas *c = new TCanvas("c","c",2500,1773);
  c->SetRightMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  h2->Draw("colz");
  c->SetLogz();

  TLatex * tex = new TLatex(0.85,0.915,"2025 (13.6 TeV)");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.10,0.915,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.08);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();
  tex = new TLatex(0.22,0.96,"Preliminary");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(52);
  tex->SetTextSize(0.0608);
  tex->SetLineWidth(2);
  c->cd();
  tex->Draw();

  c->SaveAs(filename_PDF.c_str());
  c->SaveAs(filename_C.c_str());
  c->SaveAs(filename_ROOT.c_str());
  c->SaveAs(filename_PNG.c_str());

  delete c;

  return;
}


void dEdx_VS_p_Display(bool isNewCorr)
{
  std::string corrType;
  if (isNewCorr) corrType = "NewCorr";
  else corrType = "OldCorr";

  // Setup
  TFile *ifile = new TFile("ROOT_histograms/dEdx_output.root", "READ");
    

  TH2F *dEdX0stripVsP_lowp, *dEdX0stripVsP_charge;
  if (isNewCorr) {
    dEdX0stripVsP_lowp = (TH2F*)ifile->Get("dEdX0stripVsP_lowp_NewCorr");
    dEdX0stripVsP_charge = (TH2F*)ifile->Get("dEdX0stripVsP_charge_NewCorr");
  }
  else {
    dEdX0stripVsP_lowp = (TH2F*)ifile->Get("dEdX0stripVsP_lowp_OldCorr");
    dEdX0stripVsP_charge = (TH2F*)ifile->Get("dEdX0stripVsP_charge_OldCorr");
  }

  // K and C fit
  cout << endl;
  cout << "     K and C fit:" << endl;
  double K = -1, C = -1;
  Return_K_C_params(K, C, Form("Results/K_C_fit_%s.txt",corrType.c_str()));

  TF1* Fit_pion_original = Fit_MassParametrization(0.13957, 2.5, 5, C, K);
  TF1* Fit_proton_original = Fit_MassParametrization(0.93827, 0.575, 1, C, K);
  
  TF1* Fit_pion = Fit_MassParametrization(0.13957, 0.2, 5, C, K);
  TF1* Fit_kaon = Fit_MassParametrization(0.49368, 0.2, 5, C, K);
  TF1* Fit_proton = Fit_MassParametrization(0.93827, 0.2, 5, C, K);
  TF1* Fit_deuteron = Fit_MassParametrization(1.87561, 0.5, 5, C, K);

  Display_TH2_NoFit(dEdX0stripVsP_charge, Form("Results/dEdX0stripVsP_charge_%s.pdf",corrType.c_str()), Form("Results/dEdX0stripVsP_charge_%s.C",corrType.c_str()), Form("Results/dEdX0stripVsP_charge_%s.root",corrType.c_str()), Form("Results/dEdX0stripVsP_charge_%s.png",corrType.c_str()), true);
  Display_TH2_NoFit(dEdX0stripVsP_lowp, Form("Results/dEdX0stripVsP_lowp_%s.pdf",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_%s.C",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_%s.root",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_%s.png",corrType.c_str()), false);
  Display_TH2_Fit(dEdX0stripVsP_lowp, Fit_pion_original, Fit_proton_original, Fit_pion, Fit_kaon, Fit_proton, Fit_deuteron, Form("Results/dEdX0stripVsP_lowp_FIT_%s.pdf",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_FIT_%s.C",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_FIT_%s.root",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_FIT_%s.png",corrType.c_str()));


  // Atlas fit
  cout << endl;
  cout << "     Atlas fit: in " << Form("Results/Atlas_fit_%s.txt",corrType.c_str()) << endl;
  double p1 = -1, p2 = -1, p3 = -1, p4 = -1, p5 = -1;
  Return_AtlasParams(p1, p2, p3, p4, p5, Form("Results/Atlas_fit_%s.txt",corrType.c_str()));
  
  TF1* Fit_proton_original_A = Fit_MassParametrization(0.93827, 0.55, 1.4, p1, p2, p3, p4, p5);
  TF1* Fit_pion_original_A = Fit_MassParametrization(0.13957, 2, 5, p1, p2, p3, p4, p5);
  
  TF1* Fit_pion_A = Fit_MassParametrization(0.13957, 0.2, 5, p1, p2, p3, p4, p5);
  TF1* Fit_kaon_A = Fit_MassParametrization(0.4937, 0.2, 5, p1, p2, p3, p4, p5); //0.55
  TF1* Fit_proton_A = Fit_MassParametrization(0.93827, 0.2, 5, p1, p2, p3, p4, p5);
  TF1* Fit_deuteron_A = Fit_MassParametrization(1.87561, 0.5, 5, p1, p2, p3, p4, p5);

  // save fit results
  TFile *ofile = new TFile("Results/debug.root", "RECREATE");
  ofile->cd();
  Fit_proton_A->Write("Fit_proton_A");
  Fit_deuteron_A->Write("Fit_deuteron_A");
  Fit_pion_A->Write("Fit_pion_A");
  ofile->Close();

  cout << "Atlas fit parameters: " << p1 << " " << p2 << " " << p3 << " " << p4 << " " << p5 << endl;
  Display_TH2_Fit(dEdX0stripVsP_lowp, Fit_proton_original_A, Fit_pion_original_A, Fit_pion_A, Fit_kaon_A, Fit_proton_A, Fit_deuteron_A, Form("Results/dEdX0stripVsP_lowp_FIT_A_%s.pdf",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_FIT_A_%s.C",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_FIT_A_%s.root",corrType.c_str()), Form("Results/dEdX0stripVsP_lowp_FIT_A_%s.png",corrType.c_str()));
  
  return;
}


void doDisplay() {

  dEdx_VS_p_Display(true);
  dEdx_VS_p_Display(false);

  return;
}
