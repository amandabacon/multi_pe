#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DU/PMTCalStatus.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/MC.hh> // TTree element
#include <RAT/DS/MCParticle.hh> // sub-branch of MC branch
#include <RAT/DS/EV.hh> // TTree element
#include <RAT/DS/MCPMT.hh> // sub-branch of MC branch
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
using namespace RAT;

#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TAttMarker.h>
using namespace ROOT;

#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

string filename = "snoplus_electrons.root";

void electrons2(){
    gStyle->SetOptStat(0); // no stat box
    gStyle->SetOptFit(0); // no fit info box
    gStyle->GetAttDate()->SetTextColor(1);
    gStyle->SetOptTitle(0); // no title
    gStyle->SetLabelFont(132,"XYZ");
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132,"XYZ");
    gStyle->SetLegendFont(132);

    gROOT->ForceStyle();

    TH2D* nhit_energy = new TH2D("","",100,0,10,100,0,6000);
    TH1D* ratio_1_norm = new TH1D("","",55,0.,2);
    TH1D* ratio_2_norm = new TH1D("","",55,0.,2);
    TH1D* ratio_3_norm = new TH1D("","",55,0.,2);
    TH1D* log_likelihood_1 = new TH1D("","",55,0.,2);
    TH1D* log_likelihood_2 = new TH1D("","",55,0.,2);
    TH1D* log_likelihood_3 = new TH1D("","",55,0.,2);

    RAT::DU::DSReader dsreader(filename);

    const RAT::DU::PMTCalStatus& pmtCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
    const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
//    cout << "There are " << pmtInfo.GetCount() << " PMTs" << endl; // 9728
    for(unsigned int iPMT = 0; iPMT < pmtInfo.GetCount(); iPMT++){
        TVector3 pmtPos = pmtInfo.GetPosition(iPMT);
//        pmtPos.Print();
    }

    // empty arrays to fill later
    double charges[2][9728]; // 9728 total, 9394 inward-facing
    int pmts[9728];
    int mcpmts[9728];
    int NHits[9728];

    // get the first electron only
    const RAT::DS::Entry& rds = dsreader.GetEntry(0);
    const RAT::DS::MC& rmc = rds.GetMC();
    const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
    const RAT::DS::EV& rev = rds.GetEV(0);
    int NHit = rev.GetNhits();
    const RAT::DS::MCPMT& mcPMT = rmc.GetMCPMT(0);
    double energy = rmc.GetScintEnergyDeposit();

    for(unsigned int i = 0; i < dsreader.GetEntryCount(); i++){
        const RAT::DS::Entry& rds = dsreader.GetEntry(i);
        const RAT::DS::MC& rmc = rds.GetMC();
        const RAT::DS::MCParticle& rmcparticle = rmc.GetMCParticle(0);
        TVector3 mcPos = rmcparticle.GetPosition(); // ini. pos. of 1st particle
        double mc_r = mcPos.Mag(); // gets radius

        // reset charges array
	for(unsigned int j = 0; j < 2; j++){
	    for(unsigned int k = 0; k < 9728; k++){
                charges[j][k] = 0;
            }
        }

        // reset pmt ID array
        for(unsigned int m = 0; m < 9728; m++){
            pmts[m] = 0;
        }

        // reset MC pmt ID array
        for(unsigned int n = 0; n < 9728; n++){
            mcpmts[n] = 0;
        }

        // reset NHits array
        for(unsigned int l = 0; l < 9728; l++){
            NHits[l] = 0;
        }

        unsigned int nevC = rds.GetEVCount(); // 2000 sim. evs--1191+6+803=2000
        if(nevC != 1) continue; // no retriggers--2 (and none of the non-triggered events--0)
        const RAT::DS::EV& rev = rds.GetEV(0); // get the first trigger
        const RAT::DS::CalPMTs& calpmts = rev.GetCalPMTs();
        NHits[rev.GetNhits()] = rev.GetNhits();
        int NHit = rev.GetNhits();
        for(unsigned int ipmt = 0; ipmt < calpmts.GetCount(); ipmt++){
            const RAT::DS::PMTCal& pmt = calpmts.GetPMT(ipmt);
            unsigned int status = pmtCalStatus.GetHitStatus(pmt);
            if(status != 0) continue;
            if(pmtInfo.GetType(pmt.GetID()) != RAT::DU::PMTInfo::NORMAL && pmtInfo.GetType(pmt.GetID()) != RAT::DU::PMTInfo::HQE) continue;
            if(pmtCalStatus.GetChannelStatus(pmt.GetID()) == 0 || pmt.GetID() != pmt.GetID()) continue;
            charges[0][pmt.GetID()] = pmt.GetQHS();
            charges[1][pmt.GetID()] = pmt.GetQHL();
            pmts[pmt.GetID()] = pmt.GetID();
        } // ipmt

        for(unsigned int iMCpmt = 0; iMCpmt < rmc.GetMCPMTCount(); iMCpmt++){
            const RAT::DS::MCPMT& mcPMT = rmc.GetMCPMT(iMCpmt);
            double energy = rmc.GetScintEnergyDeposit();
            if(pmtInfo.GetType(mcPMT.GetID()) != RAT::DU::PMTInfo::NORMAL && pmtInfo.GetType(mcPMT.GetID()) != RAT::DU::PMTInfo::HQE) continue;
            if(pmtCalStatus.GetChannelStatus(mcPMT.GetID()) == 0 || mcPMT.GetID() != mcPMT.GetID()) continue; // pmt.GetID() should be the same as mcPMT.GetID()
            mcpmts[mcPMT.GetID()] = mcPMT.GetID();
            if(mcPMT.GetMCPECount() == 1){
                ratio_1_norm->Fill(charges[0][mcpmts[mcPMT.GetID()]]/charges[1][mcpmts[mcPMT.GetID()]]);
                for(unsigned int i = 0; i < 55; i++){
                    ratio_1_norm->GetBinContent(i);
                    unsigned int ratio_1 = ratio_1_norm->GetBinContent(i);
                }
            } else if(mcPMT.GetMCPECount() == 2){
                ratio_2_norm->Fill(charges[0][mcpmts[mcPMT.GetID()]]/charges[1][mcpmts[mcPMT.GetID()]]);
                for(unsigned int j = 0; j < 55; j++){
                    ratio_2_norm->GetBinContent(j);
                    unsigned int ratio_2 = ratio_2_norm->GetBinContent(j);
                }
            } else if(mcPMT.GetMCPECount() == 3){
                ratio_3_norm->Fill(charges[0][mcpmts[mcPMT.GetID()]]/charges[1][mcpmts[mcPMT.GetID()]]);
                for(unsigned int k = 0; k < 55; k++){
                    ratio_3_norm->GetBinContent(k);
                    unsigned int ratio_3 = ratio_3_norm->GetBinContent(k);
                }
            }                
        } // for(unsigned int iMCpmt...

        for(unsigned int iMCpmt2 = 0; iMCpmt2 < rmc.GetMCPMTCount(); iMCpmt2++){
            const RAT::DS::MCPMT& mcPMT = rmc.GetMCPMT(iMCpmt2);
            double energy = rmc.GetScintEnergyDeposit();
            if(pmtInfo.GetType(mcPMT.GetID()) != RAT::DU::PMTInfo::NORMAL && pmtInfo.GetType(mcPMT.GetID()) != RAT::DU::PMTInfo::HQE) continue;
            if(pmtCalStatus.GetChannelStatus(mcPMT.GetID()) == 0 || mcPMT.GetID() != mcPMT.GetID()) continue; // pmt.GetID() should be the same as mcPMT.GetID()
            mcpmts[mcPMT.GetID()] = mcPMT.GetID();
            if(mcPMT.GetMCPECount() == 1){
//                log_likelihood_1->Fill();
            } else if(mcPMT.GetMCPECount() == 2){
//                log_likelihood_2->Fill();
            } else if(mcPMT.GetMCPECount() == 3){
//                log_likelihood_3->Fill();
            }
        } // for(unsigned int iMCpmt2...
    } // for(unsigned int GetEntry...
    nhit_energy->Fill(energy,NHit);

    // lowered entry count by 1,095,254, so 1 mil for each hypothesis 
    // removed the doubling due to retirggers and loops
    TCanvas *C6 = new TCanvas("C6");
    C6->UseCurrentStyle();
    C6->Update();
    ratio_1_norm->Scale(1/ratio_1_norm->GetEntries());
    ratio_1_norm->SetXTitle("Prompt Charge/Late Charge");
    ratio_1_norm->Draw();
    C6->Update();
    TLegend *legend6 = new TLegend(0.7,0.7,0.89,0.89); // R,B,L,T
    legend6->SetBorderSize(0);
    legend6->AddEntry(ratio_1_norm, "One p.e. 25 ns");
    legend6->Draw();
    C6->Update();
//    ratio_1_norm->SaveAs("ratio_1_25ns_norm.root");
    unsigned int nbins = ratio_1_norm->GetNbinsX(); // 55, checks out with bin assignment in histogram at top
    cout << nbins << endl;
    unsigned int binmax = 0; // works, but only at 0...
    cout << ratio_1_norm->GetBinContent(25) << endl; // have to put in bin to get bin content, which is y-axis value (normalized, so probability)
    for(unsigned int i = 0; i < nbins; i++){
        ratio_1_norm->GetBinContent(i);
        cout << "Bin " << i << " Y-Value = " << ratio_1_norm->GetBinContent(i) << endl; // checked using TBrowser to be sure it grabs the correct bin content for each bin
        if(ratio_1_norm->GetBinContent(i) == binmax){ 
            cout << "Max X = " << ratio_1_norm->GetXaxis()->GetBinCenter(i) << endl;
        }
    }

    TCanvas *C7 = new TCanvas("C7");
    C7->UseCurrentStyle();
    C7->Update();
    ratio_2_norm->Scale(1/ratio_2_norm->GetEntries());
    ratio_2_norm->SetXTitle("Prompt Charge/Late Charge");
    ratio_2_norm->Draw();
    C7->Update();
    TLegend *legend7 = new TLegend(0.7,0.7,0.89,0.89); // R,B,L,T
    legend7->SetBorderSize(0);
    legend7->AddEntry(ratio_2_norm, "Two p.e. 25 ns");
    legend7->Draw();
    C7->Update();
//    ratio_2_norm->SaveAs("ratio_2_25ns_norm.root");
    for(unsigned int j = 0; j < nbins; j++){
        ratio_2_norm->GetBinContent(j);
        cout << "Bin " << j << "Y-Value = " << ratio_2_norm->GetBinContent(j) << endl;
        if(ratio_2_norm->GetBinContent(j) == binmax) 
            cout << "Max X = " << ratio_2_norm->GetXaxis()->GetBinCenter(j) << endl;
    }

    TCanvas *C8 = new TCanvas("C8");
    C8->UseCurrentStyle();
    C8->Update();
    ratio_3_norm->Scale(1/ratio_3_norm->GetEntries());
    ratio_3_norm->SetXTitle("Prompt Charge/Late Charge");
    ratio_3_norm->Draw();
    C8->Update();
    TLegend *legend8 = new TLegend(0.7,0.7,0.89,0.89); // R,B,L,T
    legend8->SetBorderSize(0);
    legend8->AddEntry(ratio_3_norm, "Three p.e. 25 ns");
    legend8->Draw();
    C8->Update();
//    ratio_3_norm->SaveAs("ratio_3_25ns_norm.root");
    for(unsigned int k = 0; k < nbins; k++){
        ratio_3_norm->GetBinContent(k);
        cout << "Bin " << k << "Y-Value = " << ratio_3_norm->GetBinContent(k) << endl;
        if(ratio_3_norm->GetBinContent(k) == binmax)  
            cout << "Max X = " << ratio_3_norm->GetXaxis()->GetBinCenter(k) << endl;
    }

    TCanvas *C13 = new TCanvas("C13");
    C13->UseCurrentStyle();
    C13->Update();
    nhit_energy->SetTitleOffset(1.3,"y");
    nhit_energy->SetYTitle("Nhit");
    nhit_energy->SetXTitle("Energy [MeV]");
    nhit_energy->SetMarkerStyle(20);
    nhit_energy->Draw();
    C13->Update();
//    nhit_energy->SaveAs("nhit_energy_25ns.root");
} // electrons()
