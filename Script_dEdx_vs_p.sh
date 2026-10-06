#!/bin/bash

root -l<<EOC
.L code_dedx.C
TChain *chain_dedx = new TChain("stage/ttree");

chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunCv1_251031_101954/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunCv2_251031_102006/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunD_251031_102017/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunD_251031_102017/0001/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunE_251031_102029/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunFv1_251031_102041/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunFv2_251031_102053/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025G_fix/251209_091842/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025G_fix/251209_091842/0001/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunB_251031_102118/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunCv1_251031_102130/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunCv2_251031_102142/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunD_251031_102154/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunD_251031_102154/0001/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunE_251031_102206/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunFv1_251031_102218/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunFv2_251031_102230/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu1_Run2025G_fix/251209_091910/0000/*.root");
chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu1_Run2025G_fix/251209_091910/0001/*.root");

run2analysis t(chain_dedx);
t.Loop();
EOC


# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunCv1_251031_101954/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunCv2_251031_102006/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunD_251031_102017/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunD_251031_102017/0001/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunE_251031_102029/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunFv1_251031_102041/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025B/RunFv2_251031_102053/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025G_fix/251209_091842/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon0/Skim_ZMu0_Run2025G_fix/251209_091842/0001/*.root");


# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunB_251031_102118/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunCv1_251031_102130/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunCv2_251031_102142/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunD_251031_102154/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunD_251031_102154/0001/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunE_251031_102206/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunFv1_251031_102218/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu0_Run2025B/RunFv2_251031_102230/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu1_Run2025G_fix/251209_091910/0000/*.root");
# chain_dedx->Add("/scratch/ui14_2/ccollard/HSCP_prod/prodOct2025_CMSSW_15_0_15p4/Muon1/Skim_ZMu1_Run2025G_fix/251209_091910/0001/*.root");