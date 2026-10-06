#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include <iostream>


void SaveCanvasFit (TFile *file, TString name, TF1 *ffit, TH1D *h1, TH2 *h2, bool logscale, bool ifTH2)
{
  gROOT->SetBatch(kTRUE);

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  gPad->SetGrid();
  gStyle->SetOptStat(0);
  if (logscale) c->SetLogy();
  gPad->SetGrid(0, 0);
  h1->SetMaximum(h1->GetMaximum()*1.2);
  h1->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h1->GetYaxis()->SetTitle("Entries");
  ffit->SetLineColor(kRed);
  
  if (ifTH2)
  {
    h2->GetXaxis()->SetTitle("p [GeV/c]");
    h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
    h2->GetZaxis()->SetTitle("Entries");
    c->SetRightMargin(0.15);
    h2->Draw("COLZ");
    c->SetLogz();
  }
  h1->Draw("SAME E0");
  ffit->Draw("SAME");
  
  file->cd();
  c->Write();

  cout << "Canvas " << name << " saved in " << file->GetName() << endl;

  delete c;

  return;
}

void SaveCanvasFit (TFile *file, TString name, TF1 *ffit, TH1D *h1, TH2 *h2, bool logscale, bool ifTH2, std::string legend_text)
{
  gROOT->SetBatch(kTRUE);

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  gPad->SetGrid();
  gStyle->SetOptStat(0);
  if (logscale) c->SetLogy();
  gPad->SetGrid(0, 0);
  h1->SetMaximum(h1->GetMaximum()*1.2);
  h1->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h1->GetYaxis()->SetTitle("Entries");
  ffit->SetLineColor(kRed);
  
  if (ifTH2)
  {
    h2->GetXaxis()->SetTitle("p [GeV/c]");
    h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
    h2->GetZaxis()->SetTitle("Entries");
    c->SetRightMargin(0.15);
    h2->Draw("COLZ");
    c->SetLogz();
  }
  h1->Draw("SAME E0");
  ffit->Draw("SAME");
  
  TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->AddEntry(ffit, legend_text.c_str(), "l");
  legend->Draw();
  
  file->cd();
  c->Write();

  cout << "Canvas " << name << " saved in " << file->GetName() << endl;

  delete c;

  return;
}

void SaveCanvasFit (TFile *file, TString name, TF1 *ffit, TH1D *h1)
{
  gROOT->SetBatch(kTRUE);

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  gStyle->SetOptStat(0);
  ffit->SetLineColor(kRed);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);
  h1->GetYaxis()->SetRangeUser(0, 14);
  h1->Draw("SAME P0");
  ffit->Draw("SAME");
  
  file->cd();
  c->Write();
  cout << "Canvas " << name << " saved in " << file->GetName() << endl;

  return;
}

void SaveCanvasProtonPion (TFile *file, TString name, TH1D *h1, TH2 *h2, float start_fit, float end_fit)
{
  gROOT->SetBatch(kTRUE);

  // cloner h1 et mettre que les bins avec p entre start et end
  TH1D *h1_clone = (TH1D*)h1->Clone("h1_clone");
  h1_clone->Reset();
  for (int i=1; i<=h1->GetNbinsX(); ++i)
  {
    if (h1->GetBinCenter(i) >= start_fit && h1->GetBinCenter(i) <= end_fit)
    {
      h1_clone->SetBinContent(i, h1->GetBinContent(i));
      h1_clone->SetBinError(i, h1->GetBinError(i));
    }
  }

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);
  h1->SetMarkerColor(kBlack);
  h1_clone->SetMarkerStyle(20);
  h1_clone->SetMarkerSize(0.5);
  h1_clone->SetMarkerColor(kRed);
  h2->Draw("COLZ");
  h1->Draw("SAME P0");
  h1_clone->Draw("SAME P0");
  c->SetLogz();

  TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->AddEntry(h1_clone, "Points used for the f(#beta#gamma) fit", "p");
  legend->Draw();
  
  file->cd();
  c->Write();
  cout << "Canvas " << name << " saved in " << file->GetName() << endl;

  return;
}

TF1 *FitC (TH1D *h1, double p_start, double p_end, double &C, double &Cerr)
{
  TF1 *fit_C = new TF1("fit_C", "pol0", p_start, p_end);
  h1->Fit("fit_C", "RM0Q");

  C = fit_C->GetParameter(0);
  Cerr = fit_C->GetParError(0);

  return fit_C;
}

TF1 *FitK (TH1D *h1, double mass, double p_start, double p_end, double C, double &K, double &Kerr)
{
  TF1 *fit_K = new TF1("fit_K", Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  fit_K->FixParameter(1, C);
  h1->Fit("fit_K", "RM0Q");

  K = fit_K->GetParameter(0);
  Kerr = fit_K->GetParError(0);

  return fit_K;
}

Double_t AtlasFunction(Double_t *x, Double_t *par)
{
  Double_t p1 = par[0];
  Double_t p2 = par[1];
  Double_t p3 = par[2];
  Double_t p4 = par[3];
  Double_t p5 = par[4];

  Float_t bg = x[0];    // beta*gamma

  Float_t term1 = pow( ( sqrt(pow(bg,4) + 4*bg*bg) - bg*bg )/2 , p2/2);

  return p1 * term1 * log(1 + pow(p3 * bg, p4)) - p5;
}

TF1 *Fit_Atlas (TH1D *h1, double p_start, double p_end, double &p1, double &p2, double &p3, double &p4, double &p5, bool fixParams)
{
  TF1 *fit_A = new TF1("fit_A", AtlasFunction, p_start, p_end, 5);

  if (fixParams)
  {
    fit_A->FixParameter(0, p1);
    fit_A->FixParameter(1, p2);
    fit_A->FixParameter(2, p3);
    fit_A->FixParameter(3, p4);
    fit_A->FixParameter(4, p5);
  }
  else fit_A->SetParameters(p1, p2, p3, p4, p5);
  
  TFitResultPtr r = h1->Fit("fit_A", "R0QS");

  p1 = fit_A->GetParameter(0);
  p2 = fit_A->GetParameter(1);
  p3 = fit_A->GetParameter(2);
  p4 = fit_A->GetParameter(3);
  p5 = fit_A->GetParameter(4);

  cout << "Covariance matrix:" << endl;
  TMatrixDSym cov = r->GetCovarianceMatrix();
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) cout << cov(i,j) << " ";
    cout << endl;
  }

  cout << "Atlas fit:" << endl;
  cout << "           p1 = " << p1 << " +/- " << fit_A->GetParError(0) << endl;
  cout << "           p2 = " << p2 << " +/- " << fit_A->GetParError(1) << endl;
  cout << "           p3 = " << p3 << " +/- " << fit_A->GetParError(2) << endl;
  cout << "           p4 = " << p4 << " +/- " << fit_A->GetParError(3) << endl;
  cout << "           p5 = " << p5 << " +/- " << fit_A->GetParError(4) << endl;
  cout << endl;
  cout << "chi2/NDF = " << fit_A->GetChisquare() << " / " << fit_A->GetNDF() << endl;

  cout << "parameters value while imposing chi2/NDF to be equal to 1: ---> errors*sqrt(chi2/NDF) = errors*" << sqrt(fit_A->GetChisquare()/fit_A->GetNDF()) << endl;
  cout << "           p1 = " << p1 << " +/- " << fit_A->GetParError(0)*sqrt(fit_A->GetChisquare()/fit_A->GetNDF()) << endl;
  cout << "           p2 = " << p2 << " +/- " << fit_A->GetParError(1)*sqrt(fit_A->GetChisquare()/fit_A->GetNDF()) << endl;
  cout << "           p3 = " << p3 << " +/- " << fit_A->GetParError(2)*sqrt(fit_A->GetChisquare()/fit_A->GetNDF()) << endl;
  cout << "           p4 = " << p4 << " +/- " << fit_A->GetParError(3)*sqrt(fit_A->GetChisquare()/fit_A->GetNDF()) << endl;
  cout << "           p5 = " << p5 << " +/- " << fit_A->GetParError(4)*sqrt(fit_A->GetChisquare()/fit_A->GetNDF()) << endl;
  cout << endl;

  return fit_A;
}

double langaufun(double *x, double *par)
{
 
   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.
 
      // Numeric constants
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      double mpshift  = -0.22278298;       // Landau maximum location
 
      // Control constants
      double np = 100.0;      // number of convolution steps
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
 
      // Variables
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;
 
 
      // MP shift correction
      mpc = par[1] - mpshift * par[0];
 
      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];
 
      step = (xupp-xlow) / np;
 
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
 
         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }
 
      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit1peak(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF,TFitResultPtr &fitResultPtr)
{
   int i;
   char FunName[100];
   sprintf(FunName,"Fitfcn_%s",his->GetName());
   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);

   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MPV","Area","GSigma");
   for (i=0; i<4; i++) ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);

   fitResultPtr = his->Fit(FunName,"RB0SLQ");
   ffit->GetParameters(fitparams);
   for (i=0; i<4; i++) fiterrors[i] = ffit->GetParError(i);
  
   ChiSqr[0] = ffit->GetChisquare();
   NDF[0] = ffit->GetNDF();
   return (ffit);
}

TF1 *langaufit2peaks(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF,TFitResultPtr &fitResultPtr)
{
  int i;
  char FunName[100];
  sprintf(FunName,"Fitfcn_%s",his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName, [fitparams](double *x, double *p) {
  return langaufun(x, p) + langaufun(x, &p[4]);
  }, fitrange[0], fitrange[1], 8);

  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width1","MP1","Area1","GSigma1","Width2","MP2","Area2","GSigma2");
  for (i=0; i<8; i++) ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);

  fitResultPtr = his->Fit(FunName,"RB0SQ");
  ffit->GetParameters(fitparams);
  for (i=0; i<8; i++) fiterrors[i] = ffit->GetParError(i);
  
  ChiSqr[0] = ffit->GetChisquare();
  NDF[0] = ffit->GetNDF();

  return (ffit);
}

std::vector <TH1D*> ProjectTH2_eachBinX(const TH2 * h2, bool rebin)
{
  std::vector <TH1D*> h1s;
  for (int i=1; i<=h2->GetNbinsX(); ++i)
  {
    TH1D *h1 = h2->ProjectionY(Form("py_bin%d", i), i, i);
    
    if (rebin) h1 = (TH1D*)h1->Rebin(2, h1->GetName());

    h1s.push_back(h1);
  }

  return h1s;
}

void Compare_Fit(TH1D *h1, double p_start, double p_end, std::string filename)
{
  gROOT->SetBatch(kTRUE);

  TF1 *fit_gaus = new TF1("Gaussian fit", "gaus", p_start, p_end);
  fit_gaus->SetParameters(h1->Integral(), h1->GetMean(), h1->GetStdDev());
  fit_gaus->SetParNames("Normalization","Mean","StdDev");
  h1->Fit("Gaussian fit", "RM0Q");

  TF1 *fit_langaus = new TF1("Langaus fit", langaufun, p_start, p_end, 4);
  fit_langaus->SetParameters(h1->GetStdDev()/2, h1->GetMean(), h1->Integral()*0.08, h1->GetStdDev()*0.1);
  fit_langaus->SetParLimits(0, 0, 0.6);
  fit_langaus->SetParLimits(1, h1->GetMean()-0.5, h1->GetMean()+0.5);
  fit_langaus->SetParLimits(2, h1->Integral()*0.01, h1->Integral()*0.25);
  fit_langaus->SetParLimits(3, 0, h1->GetStdDev());
  fit_langaus->SetParNames("Width","MPV","Area","GSigma");
  h1->Fit("Langaus fit", "RM0Q");

  TCanvas *c = new TCanvas("Comparison of different fits", "Comparison of different fits", 2500, 1773);
  c->cd();
  gPad->SetGrid();
  gStyle->SetOptStat(0);
  h1->SetMaximum(h1->GetMaximum()*1.2);
  h1->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h1->GetYaxis()->SetTitle("Entries");
  h1->GetYaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetRangeUser(0, 7);
  h1->Draw("SAME E0");
  fit_gaus->SetLineColor(kRed);
  fit_gaus->Draw("SAME");
  fit_langaus->SetLineColor(kBlue);
  fit_langaus->Draw("SAME");
  c->SetLogy();
  gPad->SetGrid(0, 0);

  TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->AddEntry(fit_gaus, Form("Gaussian fit : #mu=%3.2f", fit_gaus->GetParameter(1)), "l");
  legend->AddEntry(fit_langaus, Form("Langaus fit : MPV=%3.2f", fit_langaus->GetParameter(1)), "l");
  legend->Draw();

  c->SaveAs(filename.c_str());
}

void FitFunction(bool isNewCorr)
{
  std::string corrType;
  if (isNewCorr) corrType = "NewCorr";
  else corrType = "OldCorr";


  // Setup
  TFile *ifile = new TFile("ROOT_histograms/dEdx_output.root", "READ");
  TFile *ofile = new TFile(Form("Results/dEdx_fit_%s.root",corrType.c_str()), "RECREATE");
  ofstream Save_K_C(Form("Results/K_C_fit_%s.txt",corrType.c_str()), std::ofstream::out);
  ofstream Save_AtlasParams(Form("Results/Atlas_fit_%s.txt",corrType.c_str()), std::ofstream::out);

  TH2F *dEdX0stripVsP_lowp;
  if (isNewCorr) dEdX0stripVsP_lowp = (TH2F*)ifile->Get("dEdX0stripVsP_lowp_NewCorr");
  else dEdX0stripVsP_lowp = (TH2F*)ifile->Get("dEdX0stripVsP_lowp_OldCorr");

  TH2F *CLONE_Rebin_dEdX0stripVsP_lowp = (TH2F*)dEdX0stripVsP_lowp->Clone("CLONE_Rebin_dEdX0stripVsP_lowp");
  std::vector <TH1D*> dEdx_lowP_C = ProjectTH2_eachBinX(dEdX0stripVsP_lowp, false);
  std::vector <TH1D*> dEdx_lowP_C_CLONE = ProjectTH2_eachBinX(CLONE_Rebin_dEdX0stripVsP_lowp, true);

  TH1D* histo_pion = new TH1D("histo_pion", "histo_pion", dEdX0stripVsP_lowp->GetNbinsX(), 0, dEdX0stripVsP_lowp->GetXaxis()->GetXmax());
  TH1D* histo_kaon = new TH1D("histo_kaon", "histo_kaon", dEdX0stripVsP_lowp->GetNbinsX(), 0, dEdX0stripVsP_lowp->GetXaxis()->GetXmax());
  TH1D* histo_proton = new TH1D("histo_proton", "histo_proton", dEdX0stripVsP_lowp->GetNbinsX(), 0, dEdX0stripVsP_lowp->GetXaxis()->GetXmax());

  int bound_2peak_1peak = 20;
  int bound_2peak = 36;
  int bound_1peak = 58;
  int bound_end = dEdx_lowP_C.size();

  std::vector <double> table_pion;
  std::vector <double> table_kaon;
  std::vector <double> table_proton;
  std::vector <double> table_p;
  for (int i=0; i<bound_2peak_1peak; i++) {table_pion.push_back(0); table_kaon.push_back(0); table_proton.push_back(0);}
  for (int i=0; i<dEdX0stripVsP_lowp->GetNbinsX(); i++) table_p.push_back(dEdX0stripVsP_lowp->GetXaxis()->GetBinCenter(i+1));


  // Iterative fit ---- START
  double startvalues_1peak[4] = {0.2, 11, 2*110000*0.40, 0.5};
  double parlimitslo_1peak[4] = {0.05, 10, 2*110000*0.005, 0};
  double parlimitshi_1peak[4] = {0.6, 12, 2*110000*1.5, 2};
  double fitrange_1peak[2] = {8, 16};
  double fitparams_1peak[4];
  double fiterrors_1peak[4];
  double ChiSqr_1peak;
  int NDF_1peak;
  TFitResultPtr fitResultPtr_1peak;

  double startvalues_1peak_pk[8] = {0.1, 3.2, 1000000*0.25, 0.2, 0.1, 4.5, 100000*0.25, 0.5};
  double parlimitslo_1peak_pk[8] = {0, 3, 1000000*0.01, 0, 0, 4.3, 100000*0.01, 0};
  double parlimitshi_1peak_pk[8] = {0.6, 3.3, 1000000, 0.5, 0.6, 5.6, 100000, 0.8};
  double fitrange_1peak_pk[2] = {2.5, 7};
  double fitparams_1peak_pk[8];
  double fiterrors_1peak_pk[8];
  double ChiSqr_1peak_pk;
  int NDF_1peak_pk;
  TFitResultPtr fitResultPtr_1peak_pk;

  double fitrange_2peak[2] = {1.5, 7};
  double startvalues_2peak[8] = {0.00379035, 3.07077, 1000000, 0.31263, 0.6, 5.03734, 100000, 1.17822e-07};
  double parlimitslo_2peak[8] = {0, 2.8, 1000000*0.1, 0, 0, 5.7, 1000000*0.05, 0};
  double parlimitshi_2peak[8] = {0.6, 3.2, 10000000, 2, 0.6, 6.2, 1000000, 2};
  double fitparams_2peak[8];
  double fiterrors_2peak[8];
  double ChiSqr_2peak;
  int NDF_2peak;
  TFitResultPtr fitResultPtr_2peak;

  // FIT
  for (int i=bound_2peak_1peak; i<bound_2peak; i++)
  {
    if (i <= bound_2peak_1peak+2) table_proton.push_back(0);
    if (i > bound_2peak_1peak+2)
    {
      TF1 *fit = langaufit1peak(dEdx_lowP_C_CLONE[i], fitrange_1peak, startvalues_1peak, parlimitslo_1peak, parlimitshi_1peak, fitparams_1peak, fiterrors_1peak, &ChiSqr_1peak, &NDF_1peak, fitResultPtr_1peak);
    
      table_proton.push_back(fitparams_1peak[1]);
      histo_proton->SetBinContent(i+1, fitparams_1peak[1]);
      histo_proton->SetBinError(i+1, fiterrors_1peak[1]);

      //SaveCanvasFit(ofile, Form("fit_%d", i), fit, dEdx_lowP_C_CLONE[i], dEdX0stripVsP_lowp, true, false);
      
      fitrange_1peak[0] -= 0.27;
      if (fitrange_1peak[0] < 0) fitrange_1peak[0] = 0;
      fitrange_1peak[1] -= 0.5;

      startvalues_1peak[0] = fitparams_1peak[0];   // Width
      parlimitslo_1peak[0] = fitparams_1peak[0]*0.95;
      parlimitshi_1peak[0] = fitparams_1peak[0]*1.05;

      startvalues_1peak[1] = fitparams_1peak[1]; // MPV
      parlimitslo_1peak[1] = fitparams_1peak[1]*0.9;
      parlimitshi_1peak[1] = fitparams_1peak[1]*1.1;

      startvalues_1peak[2] = fitparams_1peak[2];   // Integral
      parlimitslo_1peak[2] = fitparams_1peak[2]*0.9;
      parlimitshi_1peak[2] = fitparams_1peak[2]*1.4;

      startvalues_1peak[3] = fitparams_1peak[3];   // Sigma
      parlimitslo_1peak[3] = fitparams_1peak[3]*0.9;
      parlimitshi_1peak[3] = fitparams_1peak[3]*1.1;
    }
    
    
    TF1 *fit_pk = langaufit2peaks(dEdx_lowP_C[i], fitrange_1peak_pk, startvalues_1peak_pk, parlimitslo_1peak_pk, parlimitshi_1peak_pk, fitparams_1peak_pk, fiterrors_1peak_pk, &ChiSqr_1peak_pk, &NDF_1peak_pk, fitResultPtr_1peak_pk);
    table_pion.push_back(fitparams_1peak_pk[1]);
    table_kaon.push_back(fitparams_1peak_pk[5]);
    histo_pion->SetBinContent(i, fitparams_1peak_pk[1]);
    histo_pion->SetBinError(i, fiterrors_1peak_pk[1]);
    histo_kaon->SetBinContent(i, fitparams_1peak_pk[5]);
    histo_kaon->SetBinError(i, fiterrors_1peak_pk[5]);

    //SaveCanvasFit(ofile, Form("fit_%d", i), fit_pk, dEdx_lowP_C_CLONE[i], dEdX0stripVsP_lowp, true, false);

    for (int i = 0; i < 8; ++i)
    {
      startvalues_1peak_pk[i] = fitparams_1peak_pk[i];
      parlimitslo_1peak_pk[i] = fitparams_1peak_pk[i]*0.9;
      parlimitshi_1peak_pk[i] = fitparams_1peak_pk[i]*1.5;
    }

    if (i > 22) fitrange_1peak_pk[1] -= 0.3;
    
  }


  for (int i=bound_2peak; i<bound_1peak; ++i)
  {
    TF1 *fit = langaufit2peaks(dEdx_lowP_C[i], fitrange_2peak, startvalues_2peak, parlimitslo_2peak, parlimitshi_2peak, fitparams_2peak, fiterrors_2peak, &ChiSqr_2peak, &NDF_2peak, fitResultPtr_2peak);
    table_pion.push_back(fitparams_2peak[1]);
    table_kaon.push_back(fitparams_2peak[1]);
    table_proton.push_back(fitparams_2peak[5]);
    histo_pion->SetBinContent(i+1, fitparams_2peak[1]);
    histo_pion->SetBinError(i+1, fiterrors_2peak[1]);
    histo_kaon->SetBinContent(i+1, fitparams_2peak[1]);
    histo_kaon->SetBinError(i+1, fiterrors_2peak[1]);
    histo_proton->SetBinContent(i+1, fitparams_2peak[5]);
    histo_proton->SetBinError(i+1, fiterrors_2peak[5]);

    if (i < 55)
    {
      //if (i == 43) continue;
      for (int i = 0; i < 8; ++i)
      {
        startvalues_2peak[i] = fitparams_2peak[i];
        parlimitslo_2peak[i] = fitparams_2peak[i]*0.9;
        parlimitshi_2peak[i] = fitparams_2peak[i]*1.1;
      }
    }
    
    if (i < 55) fitrange_2peak[1] -= 0.1;

    //SaveCanvasFit(ofile, Form("fit_%d", i), fit, dEdx_lowP_C[i], dEdX0stripVsP_lowp, true, false);
  }


  double mean = dEdx_lowP_C[bound_1peak]->GetMean();
  double stddev = dEdx_lowP_C[bound_1peak]->GetStdDev();
  double integral = dEdx_lowP_C[bound_1peak]->Integral();
  double fitrange_1peakend[2] = {mean-2*stddev, mean+2*stddev};
  double startvalues_1peakend[4] = {stddev/2, mean, integral*0.08, stddev*0.1};
  double parlimitslo_1peakend[4] = {0, mean-0.5, integral*0.01, 0};
  double parlimitshi_1peakend[4] = {0.6, mean+0.5, integral*0.25, stddev};
  double fitparams_1peakend[4];
  double fiterrors_1peakend[4];
  double ChiSqr_1peakend;
  int NDF_1peakend;
  TFitResultPtr fitResultPtr_1peakend;
  for (int i=bound_1peak; i<bound_end; ++i)
  {
    TF1 *fit = langaufit1peak(dEdx_lowP_C[i], fitrange_1peakend, startvalues_1peakend, parlimitslo_1peakend, parlimitshi_1peakend, fitparams_1peakend, fiterrors_1peakend, &ChiSqr_1peakend, &NDF_1peakend, fitResultPtr_1peakend);
    table_pion.push_back(fitparams_1peakend[1]);
    table_kaon.push_back(fitparams_1peakend[1]);
    table_proton.push_back(fitparams_1peakend[1]);
    histo_pion->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_pion->SetBinError(i+1, fiterrors_1peakend[1]);
    histo_kaon->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_kaon->SetBinError(i+1, fiterrors_1peakend[1]);
    histo_proton->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_proton->SetBinError(i+1, fiterrors_1peakend[1]);

    for (int i = 0; i < 4; ++i)
    {
      startvalues_1peakend[i] = fitparams_1peakend[i];
      parlimitslo_1peakend[i] = fitparams_1peakend[i]*0.9;
      parlimitshi_1peakend[i] = fitparams_1peakend[i]*1.1;
    }

    //SaveCanvasFit(ofile, Form("fit_%d", i), fit, dEdx_lowP_C[i], dEdX0stripVsP_lowp, true, false);
  }


  TGraph *graph_pion = new TGraph(table_p.size(), &table_p[0], &table_pion[0]);
  TGraph *graph_kaon = new TGraph(table_p.size(), &table_p[0], &table_kaon[0]);
  TGraph *graph_proton = new TGraph(table_p.size(), &table_p[0], &table_proton[0]);

  gROOT->SetBatch(kTRUE);
  TCanvas *c_test = new TCanvas("c_test","c_test",2500,1500);
  c_test->cd();
  graph_pion->SetMarkerColor(kRed);
  graph_kaon->SetMarkerColor(kBlue);
  graph_proton->SetMarkerColor(kGreen);
  dEdX0stripVsP_lowp->Draw("COLZ");
  graph_pion->Draw("SAME *");
  graph_kaon->Draw("SAME *");
  graph_proton->Draw("SAME *");

  double C = -1, K = -1;
  double Cerr = -1, Kerr = -1;
  TF1 *fit_C = FitC(histo_pion, 2.5, 5, C, Cerr);
  TF1 *fit_K = FitK(histo_proton, 0.938, 0.55, 0.8, C, K, Kerr);   // [0.53 ; 0.73] : New correction
  cout << "Fit with K and C:" <<endl;
  cout << "                  C = " << C << " +/- " << Cerr << endl;
  cout << "                  K = " << K << " +/- " << Kerr << endl;
  SaveCanvasFit(ofile, "C_fit", fit_C, histo_pion, dEdX0stripVsP_lowp, false, true, "Fit to reference data");
  SaveCanvasFit(ofile, "K_fit", fit_K, histo_proton, dEdX0stripVsP_lowp, false, true, "Fit to reference data");
  Save_K_C << K << " " << C << endl;


  // ATLAS fit
  TH2 *dEdx_Rebin4 = dEdX0stripVsP_lowp->RebinX(8,"dEdx_Rebin4");
  std::vector <TH1D*> dEdx_Rebin4_proj = ProjectTH2_eachBinX(dEdx_Rebin4, false);
  std::vector <double> table_Atlas;

  TH1D *pion_Atlas = new TH1D("pion_Atlas", "pion_Atlas", dEdx_Rebin4_proj.size(), 0, 5);
  pion_Atlas->GetXaxis()->SetTitle("#beta#gamma");
  pion_Atlas->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  for (int i=8; i<dEdx_Rebin4_proj.size(); i++)
  {
    double mean = dEdx_lowP_C[i]->GetMean();
    double stddev = dEdx_lowP_C[i]->GetStdDev();
    double integral = dEdx_lowP_C[i]->Integral();
    double fitrange_1peakend[2] = {mean-2*stddev, mean+2*stddev};
    double startvalues_1peakend[4] = {stddev/2, mean, integral*0.08, stddev*0.1};
    double parlimitslo_1peakend[4] = {0, mean-0.5, integral*0.01, 0};
    double parlimitshi_1peakend[4] = {0.6, mean+0.5, integral*0.25, stddev};
    double fitparams_1peakend[4];
    double fiterrors_1peakend[4];
    double ChiSqr_1peakend;
    int NDF_1peakend;
    TFitResultPtr fitResultPtr_1peakend;

    TF1 *fit = langaufit1peak(dEdx_Rebin4_proj[i], fitrange_1peakend, startvalues_1peakend, parlimitslo_1peakend, parlimitshi_1peakend, fitparams_1peakend, fiterrors_1peakend, &ChiSqr_1peakend, &NDF_1peakend, fitResultPtr_1peakend);
    pion_Atlas->SetBinContent(i+1, fitparams_1peakend[1]);
    pion_Atlas->SetBinError(i+1, fiterrors_1peakend[1]);

    //SaveCanvasFit(ofile, Form("fit_%d", i), fit, dEdx_Rebin4_proj[i], dEdx_Rebin4, true, false);
  }

  double table_p_Atlas[60+168];
  for (int i = 0; i < 60; i++) table_p_Atlas[i] = 0.025*i;
  for (int i = 0; i < 168; i++) table_p_Atlas[60 + i] = 0.2*i + 1.5;

  TH1D *h_AtlasFit = new TH1D("h_AtlasFit", "h_AtlasFit", sizeof(table_p_Atlas)/sizeof(table_p_Atlas[0]) - 1 , table_p_Atlas);
  h_AtlasFit->GetXaxis()->SetTitle("#beta#gamma");
  h_AtlasFit->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  for (int i = 1; i <= histo_proton->GetNbinsX(); i++)
  {
    double p = histo_proton->GetXaxis()->GetBinCenter(i);

    if (p >= 0.5 && p <= 1.4)
    {
      int bin = h_AtlasFit->FindBin(p/0.93827);
      h_AtlasFit->SetBinContent(bin, histo_proton->GetBinContent(i));
      h_AtlasFit->SetBinError(bin, histo_proton->GetBinError(i));
    }
  }
  for (int i = 1; i <= pion_Atlas->GetNbinsX(); i++)
  {
    double p = pion_Atlas->GetXaxis()->GetBinCenter(i);

    if (p >= 2 && p <= 5)
    {
      int bin = h_AtlasFit->FindBin(p/0.13957);
      h_AtlasFit->SetBinContent(bin, pion_Atlas->GetBinContent(i));
      h_AtlasFit->SetBinError(bin, pion_Atlas->GetBinError(i));
    }
  }

  //double p1 = 0.0063056, p2 = -26.6638, p3 = 0.989189, p4 = 6.88310, p5 = -2.97688;
  double p1 = 0.0400404, p2 = -8.99996, p3 = 3.35448, p4 = 7.12258, p5 = -2.29656;
  bool fixParams = false;
  TF1 *fit_Atlas = Fit_Atlas(h_AtlasFit, 0, 35, p1, p2, p3, p4, p5, fixParams);
  cout << endl;
  SaveCanvasFit(ofile, "Atlas_fit", fit_Atlas, h_AtlasFit);

  SaveCanvasProtonPion(ofile, "proton_histo_onTH2", histo_proton, dEdX0stripVsP_lowp, 0.5, 1.45);
  SaveCanvasProtonPion(ofile, "pion_histo_onTH2", histo_pion, dEdX0stripVsP_lowp, 2, 5);
  Save_AtlasParams << p1 << " " << p2 << " " << p3 << " " << p4 << " " << p5 << endl;


  // Draw on a canvas the Atlas histo + the fit
  TCanvas* c_testb = new TCanvas("c_test_Atlas","c_test_Atlas",2500,1500);
  c_testb->cd();
  h_AtlasFit->SetMaximum(h_AtlasFit->GetMaximum()*1.2);
  h_AtlasFit->GetXaxis()->SetTitle("#beta#gamma");
  h_AtlasFit->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
  h_AtlasFit->GetXaxis()->SetRangeUser(0, 5);
  h_AtlasFit->SetMarkerStyle(20);
  h_AtlasFit->SetMarkerSize(1);
  h_AtlasFit->Draw("E1");
  fit_Atlas->SetLineColor(kRed);
  fit_Atlas->Draw("SAME");
  c_testb->SaveAs(Form("Results/Atlas_fit_%s.pdf",corrType.c_str()));

  TCanvas* c_test2 = new TCanvas("c_test2_Atlas","c_test2_Atlas",2500,1500);
  c_test2->cd();
  h_AtlasFit->GetXaxis()->SetRangeUser(5, 35);
  h_AtlasFit->GetYaxis()->SetRangeUser(2.9, 3.1);
  h_AtlasFit->Draw("E1");
  fit_Atlas->Draw("SAME");
  c_test2->SaveAs(Form("Results/Atlas_fit2_%s.pdf",corrType.c_str()));
  


  ofile->cd();
  c_test->Write();
  histo_pion->Write();
  histo_kaon->Write();
  histo_proton->Write();

  pion_Atlas->Write();
  h_AtlasFit->Write();
  fit_Atlas->Write();
  ofile->Close();

  dEdx_lowP_C.shrink_to_fit();
  dEdx_lowP_C_CLONE.shrink_to_fit();
}



void doFitOndEdx() {

  FitFunction(true);
  FitFunction(false);


  return;
}
