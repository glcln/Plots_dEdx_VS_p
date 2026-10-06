// ============================================================================
//  code_dedx.C
// ----------------------------------------------------------------------------
//  Fills the "dE/dx estimator versus track momentum" histograms from the
//  ntuples (TTree "stage/ttree"), once for each saturation correction of the
//  SiStrip cluster charge:
//    - OldCorr : the former correction, OldCorrection() below;
//    - NewCorr : the correction of the SaturationCorrection submodule,
//                ReturnCorrVec() in SaturationCorrection/CorrFunctions.h.
//
//  The dE/dx estimator of a track is the harmonic-2 mean of the charge per
//  unit path length of its strip clusters, see getdEdX() below.
//
//  Usage (see Script_dEdx_vs_p.sh), from the root of the repository:
//      .L code_dedx.C
//      code_dedx t(chain);     // chain of "stage/ttree" trees
//      t.Loop();
//
//  Output: the file kOutputFile, read by doFitOndEdx.C and doDisplay.C.
//
//  code_dedx.h is the TTree::MakeClass() header describing the tree.
// ============================================================================

#define code_dedx_cxx
#include "code_dedx.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// Included after the standard headers it needs.
#include "SaturationCorrection/CorrFunctions.h"


// ----------------------------------------------------------------------------
//                                  SETTINGS
// ----------------------------------------------------------------------------

// Output file, read by doFitOndEdx.C and doDisplay.C.
const TString kOutputFile = "ROOT_histograms/dEdx_output.root";

// Number of events to process: -1 for all of them, a positive value for a
// quick test.
const Long64_t kMaxEntries = -1;

// ADC count -> MeV conversion for the strips: 3.61 eV per electron-hole pair in
// silicon times ~265 electrons per ADC count (MeVperADCStrip in CMSSW).
const double kMeVperADCStrip = 3.61e-06*265;

// Minimum charge of the neighbouring strips in the old correction [ADC counts].
const float kOldCorrThreshold = 25;

// Minimum number of strip clusters on a track, counted before cluster cleaning.
const unsigned int kMinStripClusters = 10;


// ----------------------------------------------------------------------------
//                               dE/dx ESTIMATOR
// ----------------------------------------------------------------------------

// dE/dx estimator of a track [MeV/cm]: harmonic-2 mean of the charge per unit
// path length of its clusters,
//     Ih = ( 1/N * sum_i (dE/dx)_i^-2 )^(-1/2).
// The vectors hold one entry per strip cluster of the track:
//   charge        : cluster charge after saturation correction [ADC counts];
//   pathlength    : path length of the track in the sensor [cm];
//   bool_cleaning : false to leave the cluster out (cluster cleaning);
//   mustBeInside  : false to leave the cluster out (dedx_insideTkMod flag).
// Returns -1 if no cluster is left.
double code_dedx::getdEdX(const std::vector<float>& charge,
							 const std::vector<float>& pathlength,
							 const std::vector<bool>& bool_cleaning,
							 const std::vector<bool>& mustBeInside) {
	std::vector<double> vect;   // dE/dx of each cluster kept [MeV/cm]

	for (unsigned int h=0; h<charge.size(); h++)
	{
		if (!bool_cleaning[h]) continue;
		if (!mustBeInside[h]) continue;

		vect.push_back(kMeVperADCStrip * charge[h] / pathlength[h]);
	}

	const int size = vect.size();
	if (size == 0) return -1;

	const double expo = -2;
	double result = 0;
	for (int i=0; i<size; i++) result += pow(vect[i], expo);

	return pow(result/size, 1./expo);
}


// ----------------------------------------------------------------------------
//                          OLD SATURATION CORRECTION
// ----------------------------------------------------------------------------

// Former saturation correction, kept for the comparison with the new one.
// It concerns the clusters of 2 to 8 strips whose maximum is saturated (254 or
// 255 ADC): if the two strips next to the maximum are both above thresholdSat,
// not saturated, and differ by less than 40 ADC, the cluster is replaced by a
// single charge equal to 10 times the mean of these two strips. In every other
// case the cluster is returned unchanged.
//   Q            : strip charges of the cluster [ADC counts];
//   thresholdSat : minimum charge of the two neighbouring strips [ADC counts].
std::vector<uint16_t> code_dedx::OldCorrection(const std::vector<uint16_t>&  Q, float thresholdSat) {
  std::vector<uint16_t> QII;

	if(Q.size()<2 || Q.size()>8)
	{
		for (unsigned int i=0;i<Q.size();i++) QII.push_back((uint16_t) Q[i]);   // too large --> no correction
		return QII;
	}

	vector<uint16_t>::const_iterator mQ = max_element(Q.begin(), Q.end());

	if(*mQ>253){
		if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253) return Q;  //multiple peaks --> no correction
		if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1))<40)
		{
			QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2);
			return QII;	// saturation --> correction
		}
	}
	else return Q; // no saturation --> no x-talk inversion


 return Q; 
}


// ----------------------------------------------------------------------------
//                                 EVENT LOOP
// ----------------------------------------------------------------------------

void code_dedx::Loop() {
	if (fChain == 0) return;

	// Read only the branches used below: the tree holds many more.
	fChain->SetBranchStatus("*", 0);
	const char* usedBranches[] = {
		"npv",
		"ntracks", "track_qual2", "track_nvalidhits", "track_npixhits", "track_validfraction",
		"track_pt", "track_pterr", "track_p", "track_charge", "track_chi2", "track_dxy", "track_dz",
		"track_index_hit", "track_nhits", "track_prescale",
		"ndedxhits", "dedx_detid", "dedx_subdetid", "dedx_pathlength", "dedx_insideTkMod",
		"sclus_index_strip", "sclus_nstrip", "sclus_clusclean2",
		"nstrips", "strip_ampl"
	};
	for (const char* name : usedBranches) fChain->SetBranchStatus(name, 1);

	// --- Histograms: dE/dx estimator [MeV/cm] versus track momentum [GeV] ---
	//   _lowp   : p < 5 GeV, used by the fits and the displayed plots;
	//   _charge : charge sign times momentum, p < 5 GeV;
	//   (none)  : p up to 50 GeV;
	//   __large : p up to 4 TeV and wider dE/dx range (new correction only).

	// Old correction
	TH2D* dEdX0stripVsP_lowp_OldCorr = new TH2D("dEdX0stripVsP_lowp_OldCorr", "dEdX:Momentum [GeV]", 200,0,5, 400, 0.,20.);
	TH2D* dEdX0stripVsP_OldCorr = new TH2D("dEdX0stripVsP_OldCorr", "dEdX:Momentum [GeV]", 250,0,50, 400, 0.,20.);
	TH2D* dEdX0stripVsP_charge_OldCorr = new TH2D("dEdX0stripVsP_charge_OldCorr", "dEdX:Charge*Momentum [GeV]", 400,-5,5, 400, 0.,20.);
	dEdX0stripVsP_lowp_OldCorr->Sumw2();
	dEdX0stripVsP_OldCorr->Sumw2();
	dEdX0stripVsP_charge_OldCorr->Sumw2();

	// New correction
	TH2D* dEdX0stripVsP_lowp_NewCorr = new TH2D("dEdX0stripVsP_lowp_NewCorr", "dEdX:Momentum [GeV]", 200,0,5, 400, 0.,20.);
	TH2D* dEdX0stripVsP_NewCorr = new TH2D("dEdX0stripVsP_NewCorr", "dEdX:Momentum [GeV]", 250,0,50, 400, 0.,20.);
	TH2D* dEdX0stripVsP_charge_NewCorr = new TH2D("dEdX0stripVsP_charge_NewCorr", "dEdX:Charge*Momentum [GeV]", 400,-5,5, 400, 0.,20.);
	TH2D* dEdX0stripVsP_NewCorr__large = new TH2D("dEdX0stripVsP_NewCorr__large", "dEdX:Momentum [GeV]", 1000, 0, 4000, 500, 0.,50.);
	dEdX0stripVsP_lowp_NewCorr->Sumw2();
	dEdX0stripVsP_NewCorr->Sumw2();
	dEdX0stripVsP_charge_NewCorr->Sumw2();
	dEdX0stripVsP_NewCorr__large->Sumw2();

	// --- Output file, opened before the loop to fail early if it cannot be written ---
	TFile* OutputHisto = new TFile(kOutputFile, "RECREATE");
	if (OutputHisto->IsZombie())
	{
		std::cerr << "cannot open output file " << kOutputFile << std::endl;
		return;
	}

	Long64_t nentries = fChain->GetEntries();

	std::cout << std::endl;
	std::cout << "    Output file: " << kOutputFile << std::endl;
	std::cout << "--- Total event: " << nentries << std::endl;

	if (kMaxEntries >= 0 && kMaxEntries < nentries)
	{
		nentries = kMaxEntries;
		std::cout << "--- Test mode, events processed: " << nentries << std::endl;
	}

	for (Long64_t jentry=0; jentry<nentries; jentry++)
	{
		fChain->GetEntry(jentry);

		if (jentry%100000 == 0 && jentry!=0) std::cout << jentry << " (" << (100.*jentry)/(1.*nentries) << " %)" << std::endl;

		if (npv<1) continue;   // at least one primary vertex

		for (int itr=0; itr<ntracks; itr++)
		{
			// --- Track selection ---
			if (!track_qual2[itr]) continue;                       // track quality flag
			if (track_nvalidhits[itr] < 8) continue;               // number of valid hits
			if (track_npixhits[itr] < 2) continue;                 // number of pixel hits
			if (track_validfraction[itr] < 0.8) continue;          // fraction of valid hits
			if (track_pt[itr] < 0.5) continue;                     // transverse momentum [GeV]
			if (std::fabs(track_dxy[itr]) > 0.5) continue;         // transverse impact parameter [cm]
			if (std::fabs(track_dz[itr]) > 0.5) continue;          // longitudinal impact parameter [cm]
			if (track_pterr[itr]/track_pt[itr] > 0.25) continue;   // relative pT uncertainty
			if (track_chi2[itr] > 5) continue;                     // track chi2 variable of the ntuple

			// --- Strip clusters of the track: one entry per cluster ---
			std::vector<float> charge_corr_Old, charge_corr_New;       // charge after correction [ADC counts]
			std::vector<float> pathlength;                             // path length in the sensor [cm]
			std::vector<bool> bool_cleaning_Old, bool_cleaning_New;    // is the cluster kept in the estimator?
			std::vector<bool> mustBeInside;                            // dedx_insideTkMod flag

			for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++)
			{
				// Strips only: subdetid 3 (TIB), 4 (TID), 5 (TOB), 6 (TEC); 1 and 2 are the pixels.
				if (dedx_subdetid[iclu] < 3) continue;

				// Strip charges of the cluster, before any correction
				std::vector<uint16_t> Quncor;
				int cpt_sat = 0;   // number of saturated strips (254 or 255 ADC)
				for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++)
				{
					Quncor.push_back(strip_ampl[istrip]);
					if (strip_ampl[istrip]>253) cpt_sat++;
				}

				// Old saturation correction
				std::vector<uint16_t> QcorOld = OldCorrection(Quncor, kOldCorrThreshold);
				float SumQcorrOld = std::accumulate(QcorOld.begin(), QcorOld.end(), 0);

				// New saturation correction. IsSameCluster is set to true when the
				// cluster is returned unchanged.
				bool IsSameCluster = false;
				int layer = FindLayer(dedx_subdetid[iclu], dedx_detid[iclu]);
				std::vector<uint16_t> QcorNew = ReturnCorrVec(Quncor, layer, IsSameCluster);
				float SumQcorrNew = std::accumulate(QcorNew.begin(), QcorNew.end(), 0);

				// Cluster cleaning. A saturated cluster that the new correction leaves
				// unchanged is left out of the NewCorr estimator; the OldCorr estimator
				// keeps it.
				bool IscleanOldCorr = sclus_clusclean2[iclu];
				bool IscleanNewCorr = sclus_clusclean2[iclu];
				if (cpt_sat > 0 && IsSameCluster) IscleanNewCorr = false;

				charge_corr_Old.push_back(SumQcorrOld);
				charge_corr_New.push_back(SumQcorrNew);
				pathlength.push_back(dedx_pathlength[iclu]);
				mustBeInside.push_back(dedx_insideTkMod[iclu]);
				bool_cleaning_Old.push_back(IscleanOldCorr);
				bool_cleaning_New.push_back(IscleanNewCorr);
			} // end of the loop on the clusters of the track

			// --- dE/dx estimators ---
			double ih_OldCorr = getdEdX(charge_corr_Old, pathlength, bool_cleaning_Old, mustBeInside);
			double ih_NewCorr = getdEdX(charge_corr_New, pathlength, bool_cleaning_New, mustBeInside);

			// --- Histograms, each track weighted by its prescale ---
			// The cut is on the number of strip clusters of the track before
			// cleaning, so it is the same for both corrections.
			int presk = track_prescale[itr];

			if (charge_corr_Old.size() >= kMinStripClusters)
			{
				if (track_p[itr]<5)
				{
					dEdX0stripVsP_lowp_OldCorr->Fill(track_p[itr],ih_OldCorr,presk);
					dEdX0stripVsP_charge_OldCorr->Fill(track_p[itr]*track_charge[itr],ih_OldCorr,presk);
				}
				dEdX0stripVsP_OldCorr->Fill(track_p[itr],ih_OldCorr,presk);
			}
			if (charge_corr_New.size() >= kMinStripClusters)
			{
				if (track_p[itr]<5)
				{
					dEdX0stripVsP_lowp_NewCorr->Fill(track_p[itr],ih_NewCorr,presk);
					dEdX0stripVsP_charge_NewCorr->Fill(track_p[itr]*track_charge[itr],ih_NewCorr,presk);
				}
				dEdX0stripVsP_NewCorr->Fill(track_p[itr],ih_NewCorr,presk);
				dEdX0stripVsP_NewCorr__large->Fill(track_p[itr],ih_NewCorr,presk);
			}

		} // end of the loop on the tracks
	} // end of the loop on the events

	// --- Write the histograms ---
	OutputHisto->cd();

	dEdX0stripVsP_lowp_OldCorr->Write();
	dEdX0stripVsP_charge_OldCorr->Write();
	dEdX0stripVsP_OldCorr->Write();

	dEdX0stripVsP_lowp_NewCorr->Write();
	dEdX0stripVsP_charge_NewCorr->Write();
	dEdX0stripVsP_NewCorr->Write();
	dEdX0stripVsP_NewCorr__large->Write();

	OutputHisto->Close();
}