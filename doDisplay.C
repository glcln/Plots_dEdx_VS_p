// ============================================================================
//  doDisplay.C
// ----------------------------------------------------------------------------
//  Draws and saves the "dE/dx estimator versus momentum" plots. For each
//  saturation correction (<corr> = NewCorr or OldCorr), four plots are made:
//
//    dEdX0stripVsP_charge_<corr>      versus charge sign times momentum
//    dEdX0stripVsP_lowp_<corr>        versus momentum
//    dEdX0stripVsP_lowp_FIT_<corr>    same, with the "K and C" mass lines
//    dEdX0stripVsP_lowp_FIT_A_<corr>  same, with the "Atlas" mass lines
//
//  The mass lines are the dE/dx expected for pions, kaons, protons and
//  deuterons, from the parameters fitted by doFitOndEdx.C:
//
//    - "K and C" : dE/dx = K * m^2/p^2 + C
//    - "Atlas"   : dE/dx = p1 * X^(p2/2) * ln(1 + (p3*bg)^p4) - p5
//                  with bg = p/m and X = (sqrt(bg^4 + 4*bg^2) - bg^2) / 2
//
//  On the pion and proton lines, a momentum range is drawn in black and
//  labelled "Fit to reference data"; the rest is labelled "Extrapolation".
//
//  Usage, from the root of the repository, after code_dedx.C and
//  doFitOndEdx.C:
//      root -l -b -q doDisplay.C
//
//  Inputs : ROOT_histograms/dEdx_output.root (written by code_dedx.C),
//           Results/K_C_fit_<corr>.txt and Results/Atlas_fit_<corr>.txt
//           (written by doFitOndEdx.C).
//  Outputs: each plot in Results/, as .pdf, .C, .root and .png.
// ============================================================================

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>


// ----------------------------------------------------------------------------
//                                  SETTINGS
// ----------------------------------------------------------------------------

// Input file (written by code_dedx.C) and directory holding the fitted
// parameters (written by doFitOndEdx.C) and the plots.
const TString kInputFile = "ROOT_histograms/dEdx_output.root";
const TString kOutputDir = "Results";

// Text written at the top right of the plots.
const char *kDataLabel = "2025 (13.6 TeV)";

// Upper edge of the dE/dx axis [MeV/cm].
const double kDeDxMax = 14;

// Particle masses [GeV].
const double kMassPion = 0.13957;
const double kMassKaon = 0.49368;
const double kMassProton = 0.93827;
const double kMassDeuteron = 1.87561;

// Momentum range [GeV] over which the mass lines are drawn.
const double kLinePmin = 0.2;
const double kLinePminDeuteron = 0.5;
const double kLinePmax = 5;

// Momentum ranges [GeV] drawn in black as "Fit to reference data". They are
// set here by hand: the ranges actually used by the fits are in the SETTINGS
// of doFitOndEdx.C.
const double kRefKCPionPmin = 2.5, kRefKCPionPmax = 5;            // "K and C", pion line
const double kRefKCProtonPmin = 0.575, kRefKCProtonPmax = 1;      // "K and C", proton line
const double kRefAtlasPionPmin = 2, kRefAtlasPionPmax = 5;        // "Atlas", pion line
const double kRefAtlasProtonPmin = 0.55, kRefAtlasProtonPmax = 1.4;   // "Atlas", proton line

// Set to false to stop writing Results/debug.root, which holds the pion,
// proton and deuteron "Atlas" lines of the last correction processed.
const bool kWriteDebugFile = true;


// ----------------------------------------------------------------------------
//                                 MASS LINES
// ----------------------------------------------------------------------------

// The lines drawn on a plot.
struct MassLines
{
  TF1 *pion, *kaon, *proton, *deuteron;   // one line per species
  TF1 *pion_ref, *proton_ref;             // parts drawn as "Fit to reference data"
};

// Reads the n numbers stored in the text file `filename` into `values`.
// Returns false, with an error message, if they cannot be read.
bool ReadParameters(TString filename, double *values, int n) {

  std::ifstream infile(filename.Data());
  for (int i=0; i<n; i++)
  {
    if (!(infile >> values[i]))
    {
      std::cerr << "cannot read " << n << " parameters from " << filename << " (written by doFitOndEdx.C)" << std::endl;
      return false;
    }
  }

  return true;
}

// "K and C" mass line: K * mass^2/p^2 + C, for p_start < p < p_end.
TF1 *MassLine_KC(TString name, double mass, double p_start, double p_end, double K, double C) {

  TF1 *line = new TF1(name, Form("[0]* %f/(x*x) + [1]", mass*mass), p_start, p_end);
  line->SetParameters(K, C);

  return line;
}

// Atlas parametrisation of dE/dx as a function of the momentum x[0], with
// par = {mass, p1, p2, p3, p4, p5}. doFitOndEdx.C holds the same function
// written as a function of beta*gamma: keep the two in sync.
double AtlasFunction(double *x, double *par) {

  double m = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double p5 = par[5];

  double p = x[0];

  double term1 = pow( ( sqrt(pow(p/m,4) + 4*(p/m)*(p/m)) - (p/m)*(p/m) )/2 , p2/2);

  return p1 * term1 * log(1 + pow(p3 * p/m, p4)) - p5;
}

// "Atlas" mass line for p_start < p < p_end, with par = {p1, p2, p3, p4, p5}.
TF1 *MassLine_Atlas(TString name, double mass, double p_start, double p_end, const double *par) {

  TF1 *line = new TF1(name, AtlasFunction, p_start, p_end, 6);
  line->SetParameter(0, mass);
  for (int i=0; i<5; i++) line->SetParameter(i+1, par[i]);

  return line;
}


// ----------------------------------------------------------------------------
//                                   PLOTS
// ----------------------------------------------------------------------------

// Axis titles, sizes and range of the 2D histogram.
void SetHistogramStyle(TH2 *h2, const char *xtitle) {

  h2->SetStats(0);
  h2->SetTitle(" ");
  h2->GetXaxis()->SetTitle(xtitle);
  h2->GetYaxis()->SetTitle("dE/dx estimator [MeV/cm]");
  h2->GetZaxis()->SetTitle("Number of tracks");
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetTitleOffset(0.6);
  h2->GetZaxis()->SetTitleOffset(0.8);
  h2->GetYaxis()->SetRangeUser(0, kDeDxMax);
}

// Creates the canvas and draws the 2D histogram, with a logarithmic z axis.
TCanvas *DrawHistogram(TH2 *h2) {

  TCanvas *c = new TCanvas("c","c",2500,1773);
  c->SetRightMargin(0.15);
  h2->Draw("colz");
  c->SetLogz();

  return c;
}

// Draws the "CMS Preliminary" and data labels above the frame.
void DrawCMSLabels() {

  TLatex *tex = new TLatex(0.85,0.915,kDataLabel);
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->Draw();

  tex = new TLatex(0.10,0.915,"CMS");
  tex->SetNDC();
  tex->SetTextFont(61);
  tex->SetTextSize(0.08);
  tex->Draw();

  tex = new TLatex(0.22,0.96,"Preliminary");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(52);
  tex->SetTextSize(0.0608);
  tex->Draw();
}

// Draws the name of a species in red at the position (x, y), given as fractions
// of the canvas.
void DrawParticleLabel(double x, double y, const char *text) {

  TLatex *tex = new TLatex(x,y,text);
  tex->SetNDC();
  tex->SetTextColor(kRed);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->Draw();
}

// Saves the canvas as <basename>.pdf, .C, .root and .png.
void SaveCanvas(TCanvas *c, TString basename) {

  c->SaveAs(basename + ".pdf");
  c->SaveAs(basename + ".C");
  c->SaveAs(basename + ".root");
  c->SaveAs(basename + ".png");
}

// Plot of the 2D histogram alone. `charge` tells whether its x axis is the
// charge sign times the momentum or the momentum.
void Display_TH2_NoFit(TH2 *h2, TString basename, bool charge) {

  SetHistogramStyle(h2, charge ? "Charge sign x track momentum [GeV/c]" : "Track momentum [GeV/c]");

  TCanvas *c = DrawHistogram(h2);
  DrawCMSLabels();

  SaveCanvas(c, basename);

  delete c;
}

// Plot of the 2D histogram with the mass lines.
void Display_TH2_Fit(TH2 *h2, const MassLines &lines, TString basename) {

  SetHistogramStyle(h2, "Track momentum [GeV/c]");

  TCanvas *c = DrawHistogram(h2);
  lines.pion->Draw("same");
  lines.kaon->Draw("same");
  lines.proton->Draw("same");
  lines.deuteron->Draw("same");
  lines.pion_ref->SetLineColor(kBlack);
  lines.proton_ref->SetLineColor(kBlack);
  lines.pion_ref->Draw("same");
  lines.proton_ref->Draw("same");

  DrawCMSLabels();

  // Names of the species, at fixed positions on the canvas
  DrawParticleLabel(0.12, 0.3, "#pi");
  DrawParticleLabel(0.12, 0.86, "K");
  DrawParticleLabel(0.16, 0.86, "p");
  DrawParticleLabel(0.23, 0.86, "D");

  TLegend *leg = new TLegend(0.57,0.7,0.82,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(lines.pion_ref,"Fit to reference data","l");
  leg->AddEntry(lines.kaon,"Extrapolation","l");
  leg->Draw();

  SaveCanvas(c, basename);

  delete c;
}


// ----------------------------------------------------------------------------
//                                    MAIN
// ----------------------------------------------------------------------------

// Makes the four plots of one saturation correction: the new one if isNewCorr
// is true, the old one otherwise.
void dEdx_VS_p_Display(bool isNewCorr) {

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  TString corrType = isNewCorr ? "NewCorr" : "OldCorr";


  // ---------------------------------------------------------------- Setup ---

  TFile *ifile = TFile::Open(kInputFile, "READ");
  if (!ifile || ifile->IsZombie())
  {
    std::cerr << "cannot open input file " << kInputFile << std::endl;
    return;
  }

  TH2 *dEdX0stripVsP_lowp = dynamic_cast<TH2*>(ifile->Get("dEdX0stripVsP_lowp_" + corrType));
  TH2 *dEdX0stripVsP_charge = dynamic_cast<TH2*>(ifile->Get("dEdX0stripVsP_charge_" + corrType));
  if (!dEdX0stripVsP_lowp || !dEdX0stripVsP_charge)
  {
    std::cerr << "cannot find the " << corrType << " histograms in " << kInputFile << std::endl;
    return;
  }

  TString plotName = kOutputDir + "/dEdX0stripVsP_";   // start of the name of every plot


  // ----------------------------------------------------- Histograms alone ---

  Display_TH2_NoFit(dEdX0stripVsP_charge, plotName + "charge_" + corrType, true);
  Display_TH2_NoFit(dEdX0stripVsP_lowp, plotName + "lowp_" + corrType, false);


  // ---------------------------------------------------- "K and C" lines ---

  TString file_KC = kOutputDir + "/K_C_fit_" + corrType + ".txt";
  double KC[2] = {-1, -1};   // K, C
  if (ReadParameters(file_KC, KC, 2))
  {
    double K = KC[0], C = KC[1];
    std::cout << std::endl;
    std::cout << "     K and C fit: in " << file_KC << std::endl;
    std::cout << "K and C fit parameters: " << K << " " << C << std::endl;

    MassLines lines;
    lines.pion = MassLine_KC("KC_pion", kMassPion, kLinePmin, kLinePmax, K, C);
    lines.kaon = MassLine_KC("KC_kaon", kMassKaon, kLinePmin, kLinePmax, K, C);
    lines.proton = MassLine_KC("KC_proton", kMassProton, kLinePmin, kLinePmax, K, C);
    lines.deuteron = MassLine_KC("KC_deuteron", kMassDeuteron, kLinePminDeuteron, kLinePmax, K, C);
    lines.pion_ref = MassLine_KC("KC_pion_ref", kMassPion, kRefKCPionPmin, kRefKCPionPmax, K, C);
    lines.proton_ref = MassLine_KC("KC_proton_ref", kMassProton, kRefKCProtonPmin, kRefKCProtonPmax, K, C);

    Display_TH2_Fit(dEdX0stripVsP_lowp, lines, plotName + "lowp_FIT_" + corrType);
  }


  // ------------------------------------------------------- "Atlas" lines ---

  TString file_Atlas = kOutputDir + "/Atlas_fit_" + corrType + ".txt";
  double par[5] = {-1, -1, -1, -1, -1};   // p1 ... p5
  if (ReadParameters(file_Atlas, par, 5))
  {
    std::cout << std::endl;
    std::cout << "     Atlas fit: in " << file_Atlas << std::endl;
    std::cout << "Atlas fit parameters: " << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << std::endl;

    MassLines lines;
    lines.pion = MassLine_Atlas("Atlas_pion", kMassPion, kLinePmin, kLinePmax, par);
    lines.kaon = MassLine_Atlas("Atlas_kaon", kMassKaon, kLinePmin, kLinePmax, par);
    lines.proton = MassLine_Atlas("Atlas_proton", kMassProton, kLinePmin, kLinePmax, par);
    lines.deuteron = MassLine_Atlas("Atlas_deuteron", kMassDeuteron, kLinePminDeuteron, kLinePmax, par);
    lines.pion_ref = MassLine_Atlas("Atlas_pion_ref", kMassPion, kRefAtlasPionPmin, kRefAtlasPionPmax, par);
    lines.proton_ref = MassLine_Atlas("Atlas_proton_ref", kMassProton, kRefAtlasProtonPmin, kRefAtlasProtonPmax, par);

    if (kWriteDebugFile)
    {
      TFile *ofile = new TFile(kOutputDir + "/debug.root", "RECREATE");
      ofile->cd();
      lines.proton->Write("Fit_proton_A");
      lines.deuteron->Write("Fit_deuteron_A");
      lines.pion->Write("Fit_pion_A");
      ofile->Close();
    }

    Display_TH2_Fit(dEdX0stripVsP_lowp, lines, plotName + "lowp_FIT_A_" + corrType);
  }

  ifile->Close();
}


void doDisplay() {
  
  dEdx_VS_p_Display(true);
  dEdx_VS_p_Display(false);
}