#define run2analysis_cxx
#include "run2analysis.h"
#include "CorrFunctions.h"
#include <TMatrix.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>


double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside)
{
	std::vector<double> vect;
	double result = -1;

	for(int h=0;h<charge.size();h++)
	{
		if(subdetId[h]>=3 && !bool_cleaning[h]) continue;
		if(subdetId[h]>=3 && !mustBeInside[h]) continue;

		vect.push_back(3.61e-06*265 * charge[h] / pathlength[h]); //save charge
	}

	int size = vect.size();

	if(size>0)
	{
		//dEdx estimator : harmonic-2 mean
		result = 0;
		double expo = -2;
		for(int i = 0; i< size; i ++) result += pow(vect[i],expo);
		result = pow(result/size,1./expo);
	}
	else result = -1;

  return result;
}


std::vector<uint16_t> run2analysis::OldCorrection(const std::vector<uint16_t>&  Q, bool way, float thresholdSat)
{
  std::vector<uint16_t> QII;

	if(Q.size()<2 || Q.size()>8)
	{
		for (unsigned int i=0;i<Q.size();i++) QII.push_back((uint16_t) Q[i]);   // too large --> no correction
		return QII;
	}

	if(way)
	{
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
	}

 return Q; 
}


// 		LOOP FUNCTION
void run2analysis::Loop()
{
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   //nentries = 100000;

   // HISTOGRAMS OLD CORRECTION
   TH2D* dEdX0stripVsP_lowp_OldCorr = new TH2D("dEdX0stripVsP_lowp_OldCorr", "dEdX:Momentum [GeV]", 200,0,5, 400, 0.,20.);
   TH2D* dEdX0stripVsP_OldCorr = new TH2D("dEdX0stripVsP_OldCorr", "dEdX:Momentum [GeV]", 250,0,50, 400, 0.,20.);
   TH2D* dEdX0stripVsP_charge_OldCorr = new TH2D("dEdX0stripVsP_charge_OldCorr", "dEdX:Charge*Momentum [GeV]", 400,-5,5, 400, 0.,20.);
   dEdX0stripVsP_lowp_OldCorr->Sumw2();
   dEdX0stripVsP_OldCorr->Sumw2();
   dEdX0stripVsP_charge_OldCorr->Sumw2();

   // HISTOGRAMS NEW CORRECTION
   TH2D* dEdX0stripVsP_lowp_NewCorr = new TH2D("dEdX0stripVsP_lowp_NewCorr", "dEdX:Momentum [GeV]", 200,0,5, 400, 0.,20.);
   TH2D* dEdX0stripVsP_NewCorr = new TH2D("dEdX0stripVsP_NewCorr", "dEdX:Momentum [GeV]", 250,0,50, 400, 0.,20.);
   TH2D* dEdX0stripVsP_charge_NewCorr = new TH2D("dEdX0stripVsP_charge_NewCorr", "dEdX:Charge*Momentum [GeV]", 400,-5,5, 400, 0.,20.);
   TH2D* dEdX0stripVsP_NewCorr__large = new TH2D("dEdX0stripVsP_NewCorr__large", "dEdX:Charge*Momentum [GeV]", 1000, 0, 4000, 500, 0.,50.);
   dEdX0stripVsP_lowp_NewCorr->Sumw2();
   dEdX0stripVsP_NewCorr->Sumw2();
   dEdX0stripVsP_charge_NewCorr->Sumw2();
   dEdX0stripVsP_NewCorr__large->Sumw2();

   TString outputfilename="ROOT_histograms/test.root";
   TFile* OutputHisto = new TFile(outputfilename,"RECREATE");

   cout << endl;
   cout << "    Output file: " << outputfilename << endl;
   cout << "--- Total event: " << nentries << endl;

   //nentries = 10000;
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
		nb = fChain->GetEntry(jentry);
		nbytes += nb;

		if(jentry%100000 == 0 && jentry!=0) cout << jentry <<  " (" << (100.*jentry)/(1.*nentries) << " %)" <<endl;

		if (npv<1) continue;

		for (int itr=0; itr<ntracks; itr++)
		{
			bool selection=true;
			if (!track_qual2[itr]) selection=false;
			if (track_nvalidhits[itr] < 8) selection=false;
			if (track_npixhits[itr] < 2) selection=false;
			if (track_validfraction[itr] < 0.8) selection=false;
			if (track_pt[itr] < 0.5) selection = false;
			//if (track_p[itr] > 50) selection = false;
			if (abs(track_dxy[itr]) > 0.5) selection=false;
			if (abs(track_dz[itr]) > 0.5) selection=false;
			if (track_pterr[itr]/track_pt[itr] > 0.25) selection=false;
			if (track_chi2[itr] > 5) selection=false;


			if (!selection) continue;

			int presk = track_prescale[itr];
			std::vector <float> charge_corr_Old, charge_corr_New;
			std::vector <float> pathlength;
			std::vector <int> subdetId;
			std::vector <int> moduleGeometry;
			std::vector <bool> bool_cleaning_Old, bool_cleaning_New;
			std::vector <bool> mustBeInside;

			for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++)
			{
				bool IscleanOldCorr = sclus_clusclean2[iclu], IscleanNewCorr = sclus_clusclean2[iclu];
				bool IsSameCluster = false;

				if (dedx_subdetid[iclu] >= 3)	// IS STRIP
				{
					// SETUP
					vector<uint16_t> Quncor;
					int cpt_sat = 0;
					for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++) {
						Quncor.push_back(strip_ampl[istrip]);
						if (strip_ampl[istrip]>253) cpt_sat++;
					}
					
					// Old saturation correction
					vector<uint16_t> QcorOld = OldCorrection(Quncor, true, 25);		
					float SumQcorrOld = accumulate(QcorOld.begin(), QcorOld.end(), 0);

					// New saturation correction
					vector<uint16_t> QcorNew = ReturnCorr(Quncor, dedx_subdetid[iclu], dedx_detid[iclu], IsSameCluster);
					if (cpt_sat > 0 && IsSameCluster) IscleanNewCorr = false;	// if saturated but corrected bellow, no taken into account
					float SumQcorrNew = accumulate(QcorNew.begin(), QcorNew.end(), 0);
					
					charge_corr_Old.push_back(SumQcorrOld); charge_corr_New.push_back(SumQcorrNew);
					pathlength.push_back(dedx_pathlength[iclu]);
					subdetId.push_back(dedx_subdetid[iclu]);
					moduleGeometry.push_back(dedx_modulgeom[iclu]);
					mustBeInside.push_back(dedx_insideTkMod[iclu]);
					bool_cleaning_Old.push_back(IscleanOldCorr); bool_cleaning_New.push_back(IscleanNewCorr);
					
				}  //END IS STRIP
			} //END LOOP ON CLUSTERS OF A TRACK

			double ih_OldCorr = getdEdX(charge_corr_Old, pathlength, subdetId, moduleGeometry, bool_cleaning_Old, mustBeInside);
    		double ih_NewCorr = getdEdX(charge_corr_New, pathlength, subdetId, moduleGeometry, bool_cleaning_New, mustBeInside);

			if (charge_corr_Old.size() >= 10)	// at least 10 clusters
			{
				if (track_p[itr]<5)
				{
					dEdX0stripVsP_lowp_OldCorr->Fill(track_p[itr],ih_OldCorr,presk);
					dEdX0stripVsP_charge_OldCorr->Fill(track_p[itr]*track_charge[itr],ih_OldCorr,presk);
				}
				dEdX0stripVsP_OldCorr->Fill(track_p[itr],ih_OldCorr,presk);
			}
			if (charge_corr_New.size() >= 10)	// at least 10 clusters
			{
				if (track_p[itr]<5)
				{
					dEdX0stripVsP_lowp_NewCorr->Fill(track_p[itr],ih_NewCorr,presk);
					dEdX0stripVsP_charge_NewCorr->Fill(track_p[itr]*track_charge[itr],ih_NewCorr,presk);
				}
				dEdX0stripVsP_NewCorr->Fill(track_p[itr],ih_NewCorr,presk);
				dEdX0stripVsP_NewCorr__large->Fill(track_p[itr],ih_NewCorr,presk);
			}
    
		} //END LOOP ON ALL TRACKS
   }	//END LOOP ON ALL EVENTS

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
