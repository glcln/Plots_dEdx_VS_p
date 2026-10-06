// ============================================================================
//  doFitOndEdx.C
// ----------------------------------------------------------------------------
//  Fits the "dE/dx estimator versus momentum" histogram filled by code_dedx.C
//  and extracts the parameters of the two dE/dx parametrisations that
//  doDisplay.C uses to draw the mass lines:
//
//    - "K and C" : dE/dx = K * m^2/p^2 + C
//    - "Atlas"   : dE/dx = p1 * X^(p2/2) * ln(1 + (p3*bg)^p4) - p5
//                  with bg = p/m and X = (sqrt(bg^4 + 4*bg^2) - bg^2) / 2
//
//  Method, applied once per saturation correction (NewCorr, OldCorr):
//    1. The 2D histogram is cut in momentum slices (one per x-bin). In each
//       slice, the dE/dx peak of each visible species is fitted with a Landau
//       convoluted with a Gaussian ("langaus"), and its most probable value
//       (MPV) is taken as the dE/dx of the species at that momentum. Each fit
//       starts from the result of the previous slice. Three momentum regions
//       are treated differently, see the SETTINGS section.
//    2. C is the constant fitted to the pion MPVs at high momentum; K is then
//       fitted to the proton MPVs at low momentum, with C fixed.
//    3. The proton MPVs at low momentum and the pion MPVs at high momentum are
//       gathered as a function of bg = p/m and fitted with the Atlas function.
//
//  Usage, from the root of the repository:
//      root -l -b -q doFitOndEdx.C
//
//  Input : ROOT_histograms/dEdx_output.root (written by code_dedx.C).
//  Output, in Results/ (the directory must exist), <corr> = NewCorr or OldCorr:
//      K_C_fit_<corr>.txt     "K C", read by doDisplay.C
//      Atlas_fit_<corr>.txt   "p1 p2 p3 p4 p5", read by doDisplay.C
//      dEdx_fit_<corr>.root   MPV histograms, Atlas fit and control canvases
//      Atlas_fit_<corr>.pdf   Atlas fit, view of the proton points
//      Atlas_fit2_<corr>.pdf  Atlas fit, zoom on the pion points
// ============================================================================

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


// ----------------------------------------------------------------------------
//                                  SETTINGS
// ----------------------------------------------------------------------------

// Input file (written by code_dedx.C) and output directory.
const TString kInputFile = "ROOT_histograms/dEdx_output.root";
const TString kOutputDir = "Results";

// Set to true to also save the fit of every momentum slice in the output ROOT
// file (one canvas per fit).
const bool kSaveSliceFits = false;

// Particle masses [GeV].
const double kMassPion = 0.13957;
const double kMassProton = 0.93827;
const double kMassProtonKFit = 0.938;   // value used in the fit of K

// Momentum slices. Slice i is the x-bin i+1 of the 2D histogram, that is
// 0.025*i <= p < 0.025*(i+1) GeV with the binning of code_dedx.C (200 bins
// from 0 to 5 GeV). These boundaries, and the starting values, limits and
// ranges of the slice fits in FitFunction(), are tuned for that binning and
// for the 2025 data.
const int kSliceFirst = 20;         // p = 0.500 GeV: first slice fitted (the tracks have pT > 0.5 GeV)
const int kSliceFirstProton = 23;   // p = 0.575 GeV: first slice where the proton peak is fitted
const int kSlicePiKMerged = 36;     // p = 0.900 GeV: from here on, pions and kaons are fitted as one peak
const int kSliceAllMerged = 58;     // p = 1.450 GeV: from here on, the three species are fitted as one peak

// Momentum ranges [GeV] of the MPVs used to fit the parametrisations.
const double kCFitPmin = 2.5, kCFitPmax = 5;                   // pions, fit of C
const double kKFitPmin = 0.55, kKFitPmax = 0.8;                // protons, fit of K   // [0.53 ; 0.73] : New correction
const double kAtlasProtonPmin = 0.5, kAtlasProtonPmax = 1.4;   // protons, Atlas fit
const double kAtlasPionPmin = 2, kAtlasPionPmax = 5;           // pions, Atlas fit

// Atlas fit: the pion MPVs come from momentum slices 8 times wider (0.2 GeV),
// fitted from the wide slice kAtlasFirstSlice on (p = 1.6 GeV). The fit is done
// for bg between kAtlasBgMin and kAtlasBgMax, starting from kAtlasStart
// (p1 ... p5); kAtlasFixParams = true keeps the parameters at these values.
const int kAtlasRebinP = 8;
const int kAtlasFirstSlice = 8;
const double kAtlasBgMin = 0, kAtlasBgMax = 35;
const double kAtlasStart[5] = {0.0400404, -8.99996, 3.35448, 7.12258, -2.29656};
const bool kAtlasFixParams = false;


// ----------------------------------------------------------------------------
//                               FIT FUNCTIONS
// ----------------------------------------------------------------------------

// Landau convoluted with a Gaussian (from the ROOT "langaus" tutorial).
//   par[0] : width (scale) parameter of the Landau density
//   par[1] : most probable value (MPV, location) of the Landau density
//   par[2] : total area (integral from -inf to inf, normalisation constant)
//   par[3] : width (sigma) of the convoluted Gaussian
// In the Landau distribution (CERNLIB approximation), the maximum is located at
// x = -0.22278298 for a location parameter of 0. This shift is corrected here,
// so that the actual maximum is identical to the MPV parameter.
double langaufun(double *x, double *par) {
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

// Fits one langaus to the histogram `his` (binned likelihood fit).
//   fitrange                 : {lower, upper} bound of the fit range
//   startvalues              : starting values {Width, MPV, Area, GSigma}
//   parlimitslo, parlimitshi : lower and upper limits of the 4 parameters
//   fitparams, fiterrors     : outputs, fitted values and their uncertainties
// Returns the fitted function, which is not drawn.
TF1 *langaufit1peak(TH1D *his,
                    double *fitrange,
                    double *startvalues,
                    double *parlimitslo,
                    double *parlimitshi,
                    double *fitparams,
                    double *fiterrors) {
  TString FunName = Form("Fitfcn_%s", his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName, langaufun, fitrange[0], fitrange[1], 4);

  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");
  for (int i=0; i<4; i++) ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);

  his->Fit(ffit, "RB0LQ");
  ffit->GetParameters(fitparams);
  for (int i=0; i<4; i++) fiterrors[i] = ffit->GetParError(i);

  return ffit;
}

// Starting point of a one-langaus fit taken from the histogram itself: fit range
// of +-2 RMS around the mean, and starting values and limits of {Width, MPV,
// Area, GSigma} derived from the mean, the RMS and the integral.
void SetStartFromHistogram(TH1D *his,
                           double *fitrange,
                           double *startvalues,
                           double *parlimitslo,
                           double *parlimitshi) {

  double mean = his->GetMean();
  double stddev = his->GetStdDev();
  double integral = his->Integral();

  fitrange[0] = mean-2*stddev;
  fitrange[1] = mean+2*stddev;

  startvalues[0] = stddev/2;        parlimitslo[0] = 0;               parlimitshi[0] = 0.6;             // Width
  startvalues[1] = mean;            parlimitslo[1] = mean-0.5;        parlimitshi[1] = mean+0.5;        // MPV
  startvalues[2] = integral*0.08;   parlimitslo[2] = integral*0.01;   parlimitshi[2] = integral*0.25;   // Area
  startvalues[3] = stddev*0.1;      parlimitslo[3] = 0;               parlimitshi[3] = stddev;          // GSigma
}

// Fits the sum of two langaus to the histogram `his` (chi2 fit). Same arguments
// as langaufit1peak, with 8 parameters: {Width1, MP1, Area1, GSigma1, Width2,
// MP2, Area2, GSigma2}.
TF1 *langaufit2peaks(TH1D *his,
                     double *fitrange,
                     double *startvalues,
                     double *parlimitslo,
                     double *parlimitshi,
                     double *fitparams,
                     double *fiterrors) {

  TString FunName = Form("Fitfcn_%s", his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName, [](double *x, double *p) {
    return langaufun(x, p) + langaufun(x, &p[4]);
  }, fitrange[0], fitrange[1], 8);

  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width1","MP1","Area1","GSigma1","Width2","MP2","Area2","GSigma2");
  for (int i=0; i<8; i++) ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);

  his->Fit(ffit, "RB0Q");
  ffit->GetParameters(fitparams);
  for (int i=0; i<8; i++) fiterrors[i] = ffit->GetParError(i);

  return ffit;
}

// Fits the constant C to the MPVs of h1 for p_start < p < p_end.
TF1 *FitC (TH1D *h1, double p_start, double p_end, double &C, double &Cerr) {

  TF1 *fit_C = new TF1("fit_C", "pol0", p_start, p_end);
  h1->Fit(fit_C, "RM0Q");

  C = fit_C->GetParameter(0);
  Cerr = fit_C->GetParError(0);

  return fit_C;
}

// Fits K in K * mass^2/p^2 + C to the MPVs of h1 for p_start < p < p_end, with
// C fixed.
TF1 *FitK (TH1D *h1, double mass, double p_start, double p_end, double C, double &K, double &Kerr) {

  TF1 *fit_K = new TF1("fit_K", Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  fit_K->FixParameter(1, C);
  h1->Fit(fit_K, "RM0Q");

  K = fit_K->GetParameter(0);
  Kerr = fit_K->GetParError(0);

  return fit_K;
}

// Atlas parametrisation of dE/dx as a function of x[0] = beta*gamma, with
// par = {p1, p2, p3, p4, p5}. doDisplay.C holds the same function written as a
// function of the momentum: keep the two in sync.
double AtlasFunction(double *x, double *par) {

  double p1 = par[0];
  double p2 = par[1];
  double p3 = par[2];
  double p4 = par[3];
  double p5 = par[4];

  double bg = x[0];    // beta*gamma

  double term1 = pow( ( sqrt(pow(bg,4) + 4*bg*bg) - bg*bg )/2 , p2/2);

  return p1 * term1 * log(1 + pow(p3 * bg, p4)) - p5;
}

// Fits the Atlas function to h1 (MPVs versus beta*gamma) between p_start and
// p_end, and prints the result. p1 ... p5 hold the starting values on input and
// the fitted values on output; fixParams = true keeps them fixed.
TF1 *Fit_Atlas (TH1D *h1,
                double p_start,
                double p_end,
                double &p1,
                double &p2,
                double &p3,
                double &p4,
                double &p5,
                bool fixParams) {

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

  TFitResultPtr r = h1->Fit(fit_A, "R0QS");

  p1 = fit_A->GetParameter(0);
  p2 = fit_A->GetParameter(1);
  p3 = fit_A->GetParameter(2);
  p4 = fit_A->GetParameter(3);
  p5 = fit_A->GetParameter(4);

  std::cout << "Covariance matrix:" << std::endl;
  TMatrixDSym cov = r->GetCovarianceMatrix();
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) std::cout << cov(i,j) << " ";
    std::cout << std::endl;
  }

  std::cout << "Atlas fit:" << std::endl;
  for (int i = 0; i < 5; i++) std::cout << "           p" << i+1 << " = " << fit_A->GetParameter(i) << " +/- " << fit_A->GetParError(i) << std::endl;
  std::cout << std::endl;
  std::cout << "chi2/NDF = " << fit_A->GetChisquare() << " / " << fit_A->GetNDF() << std::endl;

  // Same uncertainties, scaled so that chi2/NDF = 1
  double scale = sqrt(fit_A->GetChisquare()/fit_A->GetNDF());
  std::cout << "parameters value while imposing chi2/NDF to be equal to 1: ---> errors*sqrt(chi2/NDF) = errors*" << scale << std::endl;
  for (int i = 0; i < 5; i++) std::cout << "           p" << i+1 << " = " << fit_A->GetParameter(i) << " +/- " << fit_A->GetParError(i)*scale << std::endl;
  std::cout << std::endl;

  return fit_A;
}

// Returns the projection on the y axis of each x-bin of h2: element i of the
// vector is the x-bin i+1 and is named "<prefix>_<i>". With rebin > 1, the bins
// of each projection are merged by groups of `rebin`.
// The prefix must be different in each call: ROOT reuses (resets and refills)
// an existing histogram when a projection is asked with a name already in use.
std::vector <TH1D*> ProjectTH2_eachBinX(const TH2 *h2, TString prefix, int rebin = 1) {

  std::vector <TH1D*> h1s;
  for (int i=1; i<=h2->GetNbinsX(); ++i)
  {
    TH1D *h1 = h2->ProjectionY(Form("%s_%d", prefix.Data(), i-1), i, i);

    if (rebin > 1) h1->Rebin(rebin);

    h1s.push_back(h1);
  }

  return h1s;
}


// ----------------------------------------------------------------------------
//                               CONTROL PLOTS
// ----------------------------------------------------------------------------

// Saves in `file` a canvas named `name` with the function ffit drawn on top of
// the histogram h1.
//   h2 == nullptr : h1 is the dE/dx distribution of a momentum slice;
//   h2 != nullptr : h1 holds MPVs versus momentum, drawn over the 2D
//                   histogram h2.
// A legend entry for ffit is added if legend_text is not empty.
void SaveCanvasFit (TFile *file,
                    TString name,
                    TF1 *ffit,
                    TH1D *h1,
                    TH2 *h2 = nullptr,
                    bool logscale = false,
                    std::string legend_text = "") {

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  if (logscale) c->SetLogy();
  ffit->SetLineColor(kRed);

  if (h2)
  {
    c->SetRightMargin(0.15);
    h2->Draw("COLZ");
    c->SetLogz();
  }
  else
  {
    h1->SetMaximum(h1->GetMaximum()*1.2);
    h1->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    h1->GetYaxis()->SetTitle("Entries");
  }
  h1->Draw("SAME E0");
  ffit->Draw("SAME");

  if (!legend_text.empty())
  {
    TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
    legend->SetLineColor(0);
    legend->SetFillColor(0);
    legend->AddEntry(ffit, legend_text.c_str(), "l");
    legend->Draw();
  }

  file->cd();
  c->Write();

  std::cout << "Canvas " << name << " saved in " << file->GetName() << std::endl;

  delete c;
}

// Saves in `file` a canvas named `name` with the Atlas fit drawn on top of the
// MPVs versus beta*gamma (h1).
void SaveCanvasAtlasFit (TFile *file, TString name, TF1 *ffit, TH1D *h1) {

  TCanvas *c = new TCanvas(name, name, 2500, 1500);
  c->cd();
  ffit->SetLineColor(kRed);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(0.5);
  h1->GetYaxis()->SetRangeUser(0, 14);
  h1->Draw("P0");
  ffit->Draw("SAME");

  file->cd();
  c->Write();
  std::cout << "Canvas " << name << " saved in " << file->GetName() << std::endl;

  delete c;
}

// Saves in `file` a canvas named `name` with the MPVs of h1 drawn over the 2D
// histogram h2; the points with start_fit <= p <= end_fit, used in the Atlas
// fit, are drawn in red.
void SaveCanvasProtonPion (TFile *file, TString name, TH1D *h1, TH2 *h2, double start_fit, double end_fit) {

  // Copy of h1 restricted to the points used in the fit
  TH1D *h1_clone = (TH1D*)h1->Clone(Form("%s_fitrange", h1->GetName()));
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
  std::cout << "Canvas " << name << " saved in " << file->GetName() << std::endl;

  delete c;
}

// Saves as `filename` the Atlas fit drawn on top of the MPVs versus beta*gamma
// (h1), for xmin < beta*gamma < xmax and ymin < dE/dx < ymax. h1 is not
// modified.
void SaveAtlasFitView (TH1D *h1, TF1 *ffit, TString filename, double xmin, double xmax, double ymin, double ymax) {

  TH1D *h1_view = (TH1D*)h1->Clone(Form("%s_view", h1->GetName()));

  TCanvas *c = new TCanvas("c_Atlas_view", "c_Atlas_view", 2500, 1500);
  c->cd();
  h1_view->GetXaxis()->SetRangeUser(xmin, xmax);
  h1_view->GetYaxis()->SetRangeUser(ymin, ymax);
  h1_view->SetMarkerStyle(20);
  h1_view->SetMarkerSize(1);
  h1_view->Draw("E1");
  ffit->SetLineColor(kRed);
  ffit->Draw("SAME");
  c->SaveAs(filename);

  delete c;
  delete h1_view;
}

// Compares a Gaussian and a langaus fit of the dE/dx distribution h1 between
// p_start and p_end, and saves the plot as `filename`. Not called by default.
void Compare_Fit(TH1D *h1, double p_start, double p_end, std::string filename) {

  gROOT->SetBatch(kTRUE);

  TF1 *fit_gaus = new TF1("Gaussian fit", "gaus", p_start, p_end);
  fit_gaus->SetParameters(h1->Integral(), h1->GetMean(), h1->GetStdDev());
  fit_gaus->SetParNames("Normalization","Mean","StdDev");
  h1->Fit(fit_gaus, "RM0Q");

  TF1 *fit_langaus = new TF1("Langaus fit", langaufun, p_start, p_end, 4);
  fit_langaus->SetParameters(h1->GetStdDev()/2, h1->GetMean(), h1->Integral()*0.08, h1->GetStdDev()*0.1);
  fit_langaus->SetParLimits(0, 0, 0.6);
  fit_langaus->SetParLimits(1, h1->GetMean()-0.5, h1->GetMean()+0.5);
  fit_langaus->SetParLimits(2, h1->Integral()*0.01, h1->Integral()*0.25);
  fit_langaus->SetParLimits(3, 0, h1->GetStdDev());
  fit_langaus->SetParNames("Width","MPV","Area","GSigma");
  h1->Fit(fit_langaus, "RM0Q");

  TCanvas *c = new TCanvas("Comparison of different fits", "Comparison of different fits", 2500, 1773);
  c->cd();
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

  TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
  legend->SetLineColor(0);
  legend->SetFillColor(0);
  legend->AddEntry(fit_gaus, Form("Gaussian fit : #mu=%3.2f", fit_gaus->GetParameter(1)), "l");
  legend->AddEntry(fit_langaus, Form("Langaus fit : MPV=%3.2f", fit_langaus->GetParameter(1)), "l");
  legend->Draw();

  c->SaveAs(filename.c_str());
}


// ----------------------------------------------------------------------------
//                                    MAIN
// ----------------------------------------------------------------------------

// Runs the whole procedure for one saturation correction: the new one if
// isNewCorr is true, the old one otherwise.
void FitFunction(bool isNewCorr) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  std::string corrType;
  if (isNewCorr) corrType = "NewCorr";
  else corrType = "OldCorr";


  // ---------------------------------------------------------------- Setup ---

  TFile *ifile = TFile::Open(kInputFile, "READ");
  if (!ifile || ifile->IsZombie())
  {
    std::cerr << "cannot open input file " << kInputFile << std::endl;
    return;
  }

  TString histName = Form("dEdX0stripVsP_lowp_%s", corrType.c_str());
  TH2 *dEdX0stripVsP_lowp = dynamic_cast<TH2*>(ifile->Get(histName));
  if (!dEdX0stripVsP_lowp)
  {
    std::cerr << "cannot find histogram " << histName << " in " << kInputFile << std::endl;
    return;
  }
  dEdX0stripVsP_lowp->GetXaxis()->SetTitle("p [GeV/c]");
  dEdX0stripVsP_lowp->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  dEdX0stripVsP_lowp->GetZaxis()->SetTitle("Entries");

  // The histograms created below belong to this file until it is closed.
  TFile *ofile = new TFile(Form("%s/dEdx_fit_%s.root", kOutputDir.Data(), corrType.c_str()), "RECREATE");
  if (ofile->IsZombie())
  {
    std::cerr << "cannot open output file " << ofile->GetName() << " (does the directory " << kOutputDir << " exist?)" << std::endl;
    return;
  }

  // dE/dx distribution of each momentum slice, with the binning of the 2D
  // histogram and with dE/dx bins twice as wide (for the proton peak at low p)
  std::vector <TH1D*> dEdx_lowP_C = ProjectTH2_eachBinX(dEdX0stripVsP_lowp, "proj");
  std::vector <TH1D*> dEdx_lowP_C_Rebin2 = ProjectTH2_eachBinX(dEdX0stripVsP_lowp, "proj_rebin2", 2);

  // MPV of each species versus momentum, one bin per slice (slice i in bin i+1)
  int nSlices = dEdX0stripVsP_lowp->GetNbinsX();
  double pMax = dEdX0stripVsP_lowp->GetXaxis()->GetXmax();
  TH1D* histo_pion = new TH1D("histo_pion", "histo_pion;p [GeV/c];dE/dx MPV [MeV/cm]", nSlices, 0, pMax);
  TH1D* histo_kaon = new TH1D("histo_kaon", "histo_kaon;p [GeV/c];dE/dx MPV [MeV/cm]", nSlices, 0, pMax);
  TH1D* histo_proton = new TH1D("histo_proton", "histo_proton;p [GeV/c];dE/dx MPV [MeV/cm]", nSlices, 0, pMax);

  // Same MPVs as arrays for the TGraphs, 0 where a species is not fitted
  std::vector <double> table_pion;
  std::vector <double> table_kaon;
  std::vector <double> table_proton;
  std::vector <double> table_p;
  for (int i=0; i<kSliceFirst; i++) {table_pion.push_back(0); table_kaon.push_back(0); table_proton.push_back(0);}
  for (int i=0; i<nSlices; i++) table_p.push_back(dEdX0stripVsP_lowp->GetXaxis()->GetBinCenter(i+1));


  // ------------------------------------------------------------ Slice fits ---
  // Parameter order: {Width, MPV, Area, GSigma}, twice for the two-peak fits.

  // Proton peak below kSlicePiKMerged (one langaus)
  double startvalues_1peak[4] = {0.2, 11, 2*110000*0.40, 0.5};
  double parlimitslo_1peak[4] = {0.05, 10, 2*110000*0.005, 0};
  double parlimitshi_1peak[4] = {0.6, 12, 2*110000*1.5, 2};
  double fitrange_1peak[2] = {8, 16};
  double fitparams_1peak[4];
  double fiterrors_1peak[4];

  // Pion and kaon peaks below kSlicePiKMerged (two langaus)
  double startvalues_1peak_pk[8] = {0.1, 3.2, 1000000*0.25, 0.2, 0.1, 4.5, 100000*0.25, 0.5};
  double parlimitslo_1peak_pk[8] = {0, 3, 1000000*0.01, 0, 0, 4.3, 100000*0.01, 0};
  double parlimitshi_1peak_pk[8] = {0.6, 3.3, 1000000, 0.5, 0.6, 5.6, 100000, 0.8};
  double fitrange_1peak_pk[2] = {2.5, 7};
  double fitparams_1peak_pk[8];
  double fiterrors_1peak_pk[8];

  // Pion+kaon peak and proton peak from kSlicePiKMerged to kSliceAllMerged (two langaus)
  double fitrange_2peak[2] = {1.5, 7};
  double startvalues_2peak[8] = {0.00379035, 3.07077, 1000000, 0.31263, 0.6, 5.03734, 100000, 1.17822e-07};
  double parlimitslo_2peak[8] = {0, 2.8, 1000000*0.1, 0, 0, 5.7, 1000000*0.05, 0};
  double parlimitshi_2peak[8] = {0.6, 3.2, 10000000, 2, 0.6, 6.2, 1000000, 2};
  double fitparams_2peak[8];
  double fiterrors_2peak[8];

  // --- Slices kSliceFirst to kSlicePiKMerged-1: the three peaks are separated ---
  for (int i=kSliceFirst; i<kSlicePiKMerged; i++)
  {
    // Proton peak: one langaus on the projection with the wider dE/dx bins
    if (i < kSliceFirstProton) table_proton.push_back(0);
    else
    {
      TF1 *fit = langaufit1peak(dEdx_lowP_C_Rebin2[i], fitrange_1peak, startvalues_1peak, parlimitslo_1peak, parlimitshi_1peak, fitparams_1peak, fiterrors_1peak);

      table_proton.push_back(fitparams_1peak[1]);
      histo_proton->SetBinContent(i+1, fitparams_1peak[1]);
      histo_proton->SetBinError(i+1, fiterrors_1peak[1]);

      if (kSaveSliceFits) SaveCanvasFit(ofile, Form("fit_proton_%d", i), fit, dEdx_lowP_C_Rebin2[i], nullptr, true);

      // Next slice: the fit range follows the peak towards lower dE/dx, and
      // the parameters start from this result and stay close to it
      fitrange_1peak[0] -= 0.27;
      if (fitrange_1peak[0] < 0) fitrange_1peak[0] = 0;
      fitrange_1peak[1] -= 0.5;

      startvalues_1peak[0] = fitparams_1peak[0];   // Width
      parlimitslo_1peak[0] = fitparams_1peak[0]*0.95;
      parlimitshi_1peak[0] = fitparams_1peak[0]*1.05;

      startvalues_1peak[1] = fitparams_1peak[1];   // MPV
      parlimitslo_1peak[1] = fitparams_1peak[1]*0.9;
      parlimitshi_1peak[1] = fitparams_1peak[1]*1.1;

      startvalues_1peak[2] = fitparams_1peak[2];   // Area
      parlimitslo_1peak[2] = fitparams_1peak[2]*0.9;
      parlimitshi_1peak[2] = fitparams_1peak[2]*1.4;

      startvalues_1peak[3] = fitparams_1peak[3];   // GSigma
      parlimitslo_1peak[3] = fitparams_1peak[3]*0.9;
      parlimitshi_1peak[3] = fitparams_1peak[3]*1.1;
    }

    // Pion and kaon peaks: two langaus
    TF1 *fit_pk = langaufit2peaks(dEdx_lowP_C[i], fitrange_1peak_pk, startvalues_1peak_pk, parlimitslo_1peak_pk, parlimitshi_1peak_pk, fitparams_1peak_pk, fiterrors_1peak_pk);
    table_pion.push_back(fitparams_1peak_pk[1]);
    table_kaon.push_back(fitparams_1peak_pk[5]);
    histo_pion->SetBinContent(i+1, fitparams_1peak_pk[1]);
    histo_pion->SetBinError(i+1, fiterrors_1peak_pk[1]);
    histo_kaon->SetBinContent(i+1, fitparams_1peak_pk[5]);
    histo_kaon->SetBinError(i+1, fiterrors_1peak_pk[5]);

    if (kSaveSliceFits) SaveCanvasFit(ofile, Form("fit_pion_kaon_%d", i), fit_pk, dEdx_lowP_C[i], nullptr, true);

    // Next slice: start from this result
    for (int j = 0; j < 8; ++j)
    {
      startvalues_1peak_pk[j] = fitparams_1peak_pk[j];
      parlimitslo_1peak_pk[j] = fitparams_1peak_pk[j]*0.9;
      parlimitshi_1peak_pk[j] = fitparams_1peak_pk[j]*1.5;
    }

    // From slice 23 on, the upper bound of the fit range is lowered at each slice
    if (i > 22) fitrange_1peak_pk[1] -= 0.3;
  }

  // --- Slices kSlicePiKMerged to kSliceAllMerged-1: pions and kaons merged ---
  // First peak: pions and kaons, second peak: protons.
  for (int i=kSlicePiKMerged; i<kSliceAllMerged; ++i)
  {
    TF1 *fit = langaufit2peaks(dEdx_lowP_C[i], fitrange_2peak, startvalues_2peak, parlimitslo_2peak, parlimitshi_2peak, fitparams_2peak, fiterrors_2peak);
    table_pion.push_back(fitparams_2peak[1]);
    table_kaon.push_back(fitparams_2peak[1]);
    table_proton.push_back(fitparams_2peak[5]);
    histo_pion->SetBinContent(i+1, fitparams_2peak[1]);
    histo_pion->SetBinError(i+1, fiterrors_2peak[1]);
    histo_kaon->SetBinContent(i+1, fitparams_2peak[1]);
    histo_kaon->SetBinError(i+1, fiterrors_2peak[1]);
    histo_proton->SetBinContent(i+1, fitparams_2peak[5]);
    histo_proton->SetBinError(i+1, fiterrors_2peak[5]);

    if (kSaveSliceFits) SaveCanvasFit(ofile, Form("fit_%d", i), fit, dEdx_lowP_C[i], nullptr, true);

    // Next slice: start from this result and shorten the fit range. The last
    // slices of the region (from 55 on) keep the settings of slice 54.
    if (i < 55)
    {
      for (int j = 0; j < 8; ++j)
      {
        startvalues_2peak[j] = fitparams_2peak[j];
        parlimitslo_2peak[j] = fitparams_2peak[j]*0.9;
        parlimitshi_2peak[j] = fitparams_2peak[j]*1.1;
      }

      fitrange_2peak[1] -= 0.1;
    }
  }

  // --- Slices from kSliceAllMerged on: a single peak for the three species ---
  // The first fit starts from the mean, RMS and integral of its slice.
  double fitrange_1peakend[2];
  double startvalues_1peakend[4];
  double parlimitslo_1peakend[4];
  double parlimitshi_1peakend[4];
  double fitparams_1peakend[4];
  double fiterrors_1peakend[4];
  SetStartFromHistogram(dEdx_lowP_C[kSliceAllMerged], fitrange_1peakend, startvalues_1peakend, parlimitslo_1peakend, parlimitshi_1peakend);
  for (int i=kSliceAllMerged; i<nSlices; ++i)
  {
    TF1 *fit = langaufit1peak(dEdx_lowP_C[i], fitrange_1peakend, startvalues_1peakend, parlimitslo_1peakend, parlimitshi_1peakend, fitparams_1peakend, fiterrors_1peakend);
    table_pion.push_back(fitparams_1peakend[1]);
    table_kaon.push_back(fitparams_1peakend[1]);
    table_proton.push_back(fitparams_1peakend[1]);
    histo_pion->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_pion->SetBinError(i+1, fiterrors_1peakend[1]);
    histo_kaon->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_kaon->SetBinError(i+1, fiterrors_1peakend[1]);
    histo_proton->SetBinContent(i+1, fitparams_1peakend[1]);
    histo_proton->SetBinError(i+1, fiterrors_1peakend[1]);

    if (kSaveSliceFits) SaveCanvasFit(ofile, Form("fit_%d", i), fit, dEdx_lowP_C[i], nullptr, true);

    // Next slice: start from this result
    for (int j = 0; j < 4; ++j)
    {
      startvalues_1peakend[j] = fitparams_1peakend[j];
      parlimitslo_1peakend[j] = fitparams_1peakend[j]*0.9;
      parlimitshi_1peakend[j] = fitparams_1peakend[j]*1.1;
    }
  }


  // -------------------------------------------------------- K and C fit ---

  double C = -1, K = -1;
  double Cerr = -1, Kerr = -1;
  TF1 *fit_C = FitC(histo_pion, kCFitPmin, kCFitPmax, C, Cerr);
  TF1 *fit_K = FitK(histo_proton, kMassProtonKFit, kKFitPmin, kKFitPmax, C, K, Kerr);
  std::cout << "Fit with K and C:" << std::endl;
  std::cout << "                  C = " << C << " +/- " << Cerr << std::endl;
  std::cout << "                  K = " << K << " +/- " << Kerr << std::endl;


  // ----------------------------------------------------------- Atlas fit ---

  // Pion MPVs in momentum slices kAtlasRebinP times wider, each fitted with
  // one langaus starting from the mean, RMS and integral of the slice
  TH2 *dEdx_RebinP = dEdX0stripVsP_lowp->RebinX(kAtlasRebinP, "dEdx_RebinP");
  std::vector <TH1D*> dEdx_RebinP_proj = ProjectTH2_eachBinX(dEdx_RebinP, "proj_rebinP");
  int nWideSlices = dEdx_RebinP_proj.size();

  TH1D *pion_Atlas = new TH1D("pion_Atlas", "pion_Atlas;p [GeV/c];dE/dx MPV [MeV/cm]", nWideSlices, 0, pMax);
  for (int i=kAtlasFirstSlice; i<nWideSlices; i++)
  {
    double fitrange[2];
    double startvalues[4];
    double parlimitslo[4];
    double parlimitshi[4];
    double fitparams[4];
    double fiterrors[4];
    SetStartFromHistogram(dEdx_RebinP_proj[i], fitrange, startvalues, parlimitslo, parlimitshi);

    TF1 *fit = langaufit1peak(dEdx_RebinP_proj[i], fitrange, startvalues, parlimitslo, parlimitshi, fitparams, fiterrors);
    pion_Atlas->SetBinContent(i+1, fitparams[1]);
    pion_Atlas->SetBinError(i+1, fiterrors[1]);

    if (kSaveSliceFits) SaveCanvasFit(ofile, Form("fit_atlas_%d", i), fit, dEdx_RebinP_proj[i], nullptr, true);
  }

  // MPVs versus beta*gamma = p/m. Bin edges: 0.025 wide up to 1.5 (proton
  // points), 0.2 wide from 1.5 to 34.9 (pion points).
  double table_p_Atlas[60+168];
  for (int i = 0; i < 60; i++) table_p_Atlas[i] = 0.025*i;
  for (int i = 0; i < 168; i++) table_p_Atlas[60 + i] = 0.2*i + 1.5;

  TH1D *h_AtlasFit = new TH1D("h_AtlasFit", "h_AtlasFit;#beta#gamma;dE/dx [MeV/cm]", sizeof(table_p_Atlas)/sizeof(table_p_Atlas[0]) - 1 , table_p_Atlas);
  for (int i = 1; i <= histo_proton->GetNbinsX(); i++)
  {
    double p = histo_proton->GetXaxis()->GetBinCenter(i);

    if (p >= kAtlasProtonPmin && p <= kAtlasProtonPmax)
    {
      int bin = h_AtlasFit->FindBin(p/kMassProton);
      h_AtlasFit->SetBinContent(bin, histo_proton->GetBinContent(i));
      h_AtlasFit->SetBinError(bin, histo_proton->GetBinError(i));
    }
  }
  for (int i = 1; i <= pion_Atlas->GetNbinsX(); i++)
  {
    double p = pion_Atlas->GetXaxis()->GetBinCenter(i);

    if (p >= kAtlasPionPmin && p <= kAtlasPionPmax)
    {
      int bin = h_AtlasFit->FindBin(p/kMassPion);
      h_AtlasFit->SetBinContent(bin, pion_Atlas->GetBinContent(i));
      h_AtlasFit->SetBinError(bin, pion_Atlas->GetBinError(i));
    }
  }

  double p1 = kAtlasStart[0], p2 = kAtlasStart[1], p3 = kAtlasStart[2], p4 = kAtlasStart[3], p5 = kAtlasStart[4];
  TF1 *fit_Atlas = Fit_Atlas(h_AtlasFit, kAtlasBgMin, kAtlasBgMax, p1, p2, p3, p4, p5, kAtlasFixParams);
  std::cout << std::endl;


  // ------------------------------------------------------------- Outputs ---

  // Parameters read by doDisplay.C
  std::ofstream Save_K_C(Form("%s/K_C_fit_%s.txt", kOutputDir.Data(), corrType.c_str()), std::ofstream::out);
  Save_K_C << K << " " << C << std::endl;
  Save_K_C.close();

  std::ofstream Save_AtlasParams(Form("%s/Atlas_fit_%s.txt", kOutputDir.Data(), corrType.c_str()), std::ofstream::out);
  Save_AtlasParams << p1 << " " << p2 << " " << p3 << " " << p4 << " " << p5 << std::endl;
  Save_AtlasParams.close();

  // MPV histograms and Atlas fit, written before the plots below change the
  // way they are drawn
  ofile->cd();
  histo_pion->Write();
  histo_kaon->Write();
  histo_proton->Write();
  pion_Atlas->Write();
  h_AtlasFit->Write();
  fit_Atlas->Write();

  // MPVs of the three species over the 2D histogram
  TGraph *graph_pion = new TGraph(table_p.size(), &table_p[0], &table_pion[0]);
  TGraph *graph_kaon = new TGraph(table_p.size(), &table_p[0], &table_kaon[0]);
  TGraph *graph_proton = new TGraph(table_p.size(), &table_p[0], &table_proton[0]);

  TCanvas *c_test = new TCanvas("c_test","c_test",2500,1500);
  c_test->cd();
  graph_pion->SetMarkerColor(kRed);
  graph_kaon->SetMarkerColor(kBlue);
  graph_proton->SetMarkerColor(kGreen);
  dEdX0stripVsP_lowp->Draw("COLZ");
  graph_pion->Draw("SAME *");
  graph_kaon->Draw("SAME *");
  graph_proton->Draw("SAME *");
  ofile->cd();
  c_test->Write();
  delete c_test;

  // Fits of the parametrisations
  SaveCanvasFit(ofile, "C_fit", fit_C, histo_pion, dEdX0stripVsP_lowp, false, "Fit to reference data");
  SaveCanvasFit(ofile, "K_fit", fit_K, histo_proton, dEdX0stripVsP_lowp, false, "Fit to reference data");
  SaveCanvasAtlasFit(ofile, "Atlas_fit", fit_Atlas, h_AtlasFit);

  // Points entering the Atlas fit
  SaveCanvasProtonPion(ofile, "proton_histo_onTH2", histo_proton, dEdX0stripVsP_lowp, kAtlasProtonPmin, kAtlasProtonPmax);
  SaveCanvasProtonPion(ofile, "pion_histo_onTH2", histo_pion, dEdX0stripVsP_lowp, kAtlasPionPmin, kAtlasPionPmax);

  // Atlas fit: view of the proton points, and zoom on the pion points
  SaveAtlasFitView(h_AtlasFit, fit_Atlas, Form("%s/Atlas_fit_%s.pdf", kOutputDir.Data(), corrType.c_str()), 0, 5, 0, 16.8);
  SaveAtlasFitView(h_AtlasFit, fit_Atlas, Form("%s/Atlas_fit2_%s.pdf", kOutputDir.Data(), corrType.c_str()), 5, 35, 2.9, 3.1);

  // Closing the files deletes the histograms that belong to them
  ofile->Close();
  ifile->Close();
}


void doFitOndEdx()
{
  FitFunction(true);
  FitFunction(false);
}