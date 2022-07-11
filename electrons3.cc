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
#include <cmath>
#include <algorithm>
using namespace std;

string filename = "snoplus_electrons_testing.root";

void electrons3(){
    unsigned int count12 = 0;
    unsigned int count13 = 0;
    unsigned int count23 = 0;
    gStyle->SetOptStat(0); // no stat box
    gStyle->SetOptFit(0); // no fit info box
    gStyle->GetAttDate()->SetTextColor(1);
    gStyle->SetOptTitle(0); // no title
    gStyle->SetLabelFont(132,"XYZ");
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132,"XYZ");
    gStyle->SetLegendFont(132);

    gROOT->ForceStyle();

    TH1D* log_likelihood_1 = new TH1D("log_likelihood_1","",55,0.,2);
    TH1D* log_likelihood_2 = new TH1D("log_likelihood_2","",55,0.,2);
    TH1D* log_likelihood_3 = new TH1D("log_likelihood_3","",55,0.,2);
    TH1D* log_likelihood_true_12 = new TH1D("log_likelihood_true_12","",55,0.,2);
    TH1D* log_likelihood_true_13 = new TH1D("log_likelihood_true_13","",55,0.,2);
    TH1D* log_likelihood_true_23 = new TH1D("log_likelihood_true_23","",55,0.,2);

    TFile* file1 = new TFile("ratio_1_25ns_norm.root");
    TH1D* ratio_1_norm = (TH1D*)file1->Get(""); // didn't give the histo a name like above
    ratio_1_norm->SetDirectory(0);

    TFile* file2 = new TFile("ratio_2_25ns_norm.root");
    TH1D* ratio_2_norm = (TH1D*)file2->Get(""); // didn't give the histo a name like above
    ratio_2_norm->SetDirectory(0);

    TFile* file3 = new TFile("ratio_3_25ns_norm.root");
    TH1D* ratio_3_norm = (TH1D*)file3->Get(""); // didn't give the histo a name like above
    unsigned int nbins = ratio_1_norm->GetNbinsX();
    double ll12[55];
    double ll13[55];
    double ll23[55];

    // reset ll12 array
    for(unsigned int z = 1; z < 55; z++){
        ll12[z] = 0;
    }

    // reset ll13 array
    for(unsigned int y = 1; y < 55; y++){
        ll13[y] = 0;
    }

    // reset ll23 array
    for(unsigned int x = 1; x < 55; x++){
        ll23[x] = 0;
    }

    for(unsigned int j = 1; j < nbins; j++){ // THIS SHOULD BE FROM 1 TO < NBINS TO REMOVE OVERFLOW BIN
        double ratio1 = ratio_1_norm->GetBinContent(j);
        double ratio2 = ratio_2_norm->GetBinContent(j);
        double ratio3 = ratio_3_norm->GetBinContent(j);
//        cout << "1 pe Bin " << j << " Y-Value = " << ratio_1_norm->GetBinContent(j) << endl;
//        cout << "2 pe Bin " << j << " Y-Value = " << ratio_2_norm->GetBinContent(j) << endl;
//        cout << "3 pe Bin " << j << " Y-Value = " << ratio_3_norm->GetBinContent(j) << endl;
        double div12 = ratio1/ratio2;
        double div13 = ratio1/ratio3;
        double div23 = ratio2/ratio3;
        double log_likelihood12 = log(div12); // natural base e logarithm
        double log_likelihood13 = log(div13); // natural base e logarithm
        double log_likelihood23 = log(div23); // natural base e logarithm
        for(unsigned int m = 0; m < log_likelihood12; m++){
            count12++;
        }
        for(unsigned int p = 0; p < log_likelihood13; p++){
            count13++;
        }
        for(unsigned int i = 0; i < log_likelihood23; i++){
            count23++;
        }
//        cout << "Division12: " << div12 << endl;
//        cout << "Division13: " << div13 << endl;
//        cout << "Division23: " << div23 << endl;
        cout << "Log Likelihood12 (not edited): " << log_likelihood12 << endl;
        cout << "Log Likelihood13 (not edited): " << log_likelihood13 << endl;
        cout << "Log Likelihood23 (not edited): " << log_likelihood23 << endl;
        if(isinf(log_likelihood12)){
            log_likelihood12 = 0;
            ll12[j] = log_likelihood12;
            cout << "Log Likelihood12 (remove inf): " << log_likelihood12 << endl;
            cout << "Check12 (remove inf): " << ll12[j] << endl;
        }
        if(isinf(log_likelihood13)){
            log_likelihood13 = 0;
            ll13[j] = log_likelihood13;
            cout << "Log Likelihood13 (remove inf): " << log_likelihood13 << endl;
            cout << "Check13 (remove inf): " << ll13[j] << endl; 
        }
        if(isinf(log_likelihood23)){
            log_likelihood23 = 0;
            ll23[j] = log_likelihood23;
            cout << "Log Likelihood23 (remove inf): " << log_likelihood23 << endl;
            cout << "Check23 (remove inf): " << ll23[j] << endl;
        }
        ll12[j] = log_likelihood12;
        cout << "LL12 (after removal): " << ll12[j] << endl;
        ll13[j] = log_likelihood13;
        cout << "LL13 (after removal): " << ll13[j] << endl;
        ll23[j] = log_likelihood23;
        cout << "LL23 (after removal): " << ll23[j] << endl;
        log_likelihood_1->SetBinContent(j,log_likelihood12);
        log_likelihood_2->SetBinContent(j,log_likelihood23);
        log_likelihood_3->SetBinContent(j,log_likelihood13);
    }
    ratio_3_norm->SetDirectory(0);

    RAT::DU::DSReader dsreader(filename);

    const RAT::DU::PMTCalStatus& pmtCalStatus = RAT::DU::Utility::Get()->GetPMTCalStatus();
    const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
    for(unsigned int iPMT = 0; iPMT < pmtInfo.GetCount(); iPMT++){
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
            double fprompt = charges[0][mcpmts[mcPMT.GetID()]]/charges[1][mcpmts[mcPMT.GetID()]];
            if(isnan(fprompt)){
                fprompt = -99999;
            }
            if(isinf(fprompt)){
                fprompt = -99999;
            }
            if(fprompt == -99999) continue;
//            cout << "This is fprompt: " << fprompt << endl;
//            double min = fprompt;
//            double max = fprompt;
//            if(max < fprompt){
//                max = fprompt;
//            } else if(min > fprompt){
//                min = fprompt;
//            }
//            cout << "min: " << min << endl;
//            cout << "max: " << max << endl;
            int j = floor((fprompt/((2.0-0.0)/55)));
            cout << "This is j: " << j << endl;
            double delta_likelihood_12 = ll12[j];
            cout << "This is delta l12: " << delta_likelihood_12 << endl;
            double delta_likelihood_13 = ll13[j];
            cout << "This is delta l13: " << delta_likelihood_13 << endl;
            double delta_likelihood_23 = ll23[j];
            cout << "This is delta l23: " << delta_likelihood_23 << endl;
            log_likelihood_true_12->Fill(delta_likelihood_12);
            log_likelihood_true_13->Fill(delta_likelihood_13);
            log_likelihood_true_23->Fill(delta_likelihood_23);
            if(mcPMT.GetMCPECount() == 1){
//                log_likelihood_true_12->Fill(delta_likelihood_12);
            } else if(mcPMT.GetMCPECount() == 2){
//                log_likelihood_true_12->Fill(delta_likelihood_12);
//                log_likelihood_true_23->Fill(delta_likelihood_23);
            } else if(mcPMT.GetMCPECount() == 3){
//                log_likelihood_true_13->Fill(delta_likelihood_13);
//                log_likelihood_true_23->Fill(delta_likelihood_23);
            }
        } // for(unsigned int iMCpmt...
    } // for(unsigned int GetEntry...
    cout << "length12 " << count12 << endl; // 39
    cout << "length13 " << count13 << endl; // 74
    cout << "length23 " << count23 << endl; // 51

    TCanvas *C14 = new TCanvas("C14");
    C14->UseCurrentStyle();
    C14->Update();
    log_likelihood_1->SetXTitle("Fprompt");
    log_likelihood_1->SetYTitle("Log Likelihood Ratio");
    log_likelihood_1->Draw();
    C14->Update();
    TLegend *legend14 = new TLegend(0.1,0.7,0.29,0.89); // R,B,L,T
    legend14->SetBorderSize(0);
    legend14->SetFillStyle(0);
    legend14->AddEntry(log_likelihood_1, "1 p.e./2 p.e. 25 ns");
    legend14->Draw();
    C14->Update();
    C14->SaveAs("log_likelihood_12_25ns.png");
    log_likelihood_1->SaveAs("log_likelihood_12_25ns.root");

    TCanvas *C15 = new TCanvas("C15");
    C15->UseCurrentStyle();
    C15->Update();
    log_likelihood_2->SetXTitle("Fprompt");
    log_likelihood_2->SetYTitle("Log Likelihood Ratio");
    log_likelihood_2->Draw();
    C15->Update();
    TLegend *legend15 = new TLegend(0.1,0.7,0.29,0.89); // R,B,L,T
    legend15->SetBorderSize(0);
    legend15->SetFillStyle(0);
    legend15->AddEntry(log_likelihood_2, "2 p.e./3 p.e. 25 ns");
    legend15->Draw();
    C15->Update();
    C15->SaveAs("log_likelihood_23_25ns.png");
    log_likelihood_2->SaveAs("log_likelihood_23_25ns.root");

    TCanvas *C16 = new TCanvas("C16");
    C16->UseCurrentStyle();
    C16->Update();
    log_likelihood_3->SetXTitle("Fprompt");
    log_likelihood_3->SetYTitle("Log Likelihood Ratio");
    log_likelihood_3->Draw();
    C16->Update();
    TLegend *legend16 = new TLegend(0.1,0.7,0.29,0.89); // R,B,L,T
    legend16->SetBorderSize(0);
    legend16->SetFillStyle(0);
    legend16->AddEntry(log_likelihood_3, "1 p.e./3 p.e. 25 ns");
    legend16->Draw();
    C16->Update();
    C16->SaveAs("log_likelihood_13_25ns.png");
    log_likelihood_3->SaveAs("log_likelihood_13_25ns.root");

    TCanvas *C17 = new TCanvas("C17");
    C17->UseCurrentStyle();
    C17->Update();
    log_likelihood_true_12->Scale(1/log_likelihood_true_12->GetEntries());
    log_likelihood_true_12->SetXTitle("LL12");
    log_likelihood_true_12->Draw();
    C17->Update();
    TLegend *legend17 = new TLegend(0.7,0.7,0.89,0.89); // R,B,L,T
    legend17->SetBorderSize(0);
    legend17->SetFillStyle(0);
    legend17->AddEntry(log_likelihood_true_12, "LL12 25 ns");
    legend17->Draw();
    C17->Update();
    C17->SaveAs("log_likelihood_true_12_25ns.png");
    log_likelihood_true_12->SaveAs("log_likelihood_true_12_25ns.root");

    TCanvas *C18 = new TCanvas("C18");
    C18->UseCurrentStyle();
    C18->Update();
    log_likelihood_true_23->Scale(1/log_likelihood_true_23->GetEntries());
    log_likelihood_true_23->SetXTitle("LL23");
    log_likelihood_true_23->Draw();
    C18->Update();
    TLegend *legend18 = new TLegend(0.7,0.7,0.89,0.89); // R,B,L,T
    legend18->SetBorderSize(0);
    legend18->SetFillStyle(0);
    legend18->AddEntry(log_likelihood_true_23, "LL23  25 ns");
    legend18->Draw();
    C18->Update();
    C18->SaveAs("log_likelihood_true_23_25ns.png");
    log_likelihood_true_23->SaveAs("log_likelihood_true_23_25ns.root");

    TCanvas *C19 = new TCanvas("C19");
    C19->UseCurrentStyle();
    C19->Update();
    log_likelihood_true_13->Scale(1/log_likelihood_true_13->GetEntries());
    log_likelihood_true_13->SetXTitle("LL13");
    log_likelihood_true_13->Draw();
    C19->Update();
    TLegend *legend19 = new TLegend(0.7,0.7,0.89,0.89); // R,B,L,T
    legend19->SetBorderSize(0);
    legend19->SetFillStyle(0);
    legend19->AddEntry(log_likelihood_true_13, "LL13 25 ns");
    legend19->Draw();
    C19->Update();
    C19->SaveAs("log_likelihood_true_13_25ns.png");
    log_likelihood_true_13->SaveAs("log_likelihood_true_13_25ns.root");
} // electrons3()
