

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TTree.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
//#include "reader.h"

// Header file for the classes stored in the TTree if any.
#include "vector"

void compare_neutrons()
{
    // Load the HIPO library
    gSystem->Load("libhipo4");

    // Create a chain for each file
    // Create a chain for each file
    TChain *chain1 = new TChain("chain1", "");
    TChain *chain2 = new TChain("chain2", "");
    TChain *chain3 = new TChain("chain3", "");


   chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output0/dst-0.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output1/dst-1.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output2/dst-2.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output3/dst-3.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output4/dst-4.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output5/dst-5.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output6/dst-6.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output7/dst-7.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output8/dst-8.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output9/dst-9.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output10/dst-10.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output11/dst-11.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output12/dst-12.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output13/dst-13.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output14/dst-14.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output15/dst-15.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output16/dst-16.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output17/dst-17.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output18/dst-18.hipo");
    chain1->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testOSG/output19/dst-19.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output0/dst-0.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output1/dst-1.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output2/dst-2.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output3/dst-3.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output4/dst-4.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output5/dst-5.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output6/dst-6.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output7/dst-7.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output8/dst-8.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output9/dst-9.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output10/dst-10.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output11/dst-11.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output12/dst-12.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output13/dst-13.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output14/dst-14.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output15/dst-15.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output16/dst-16.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output17/dst-17.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output18/dst-18.hipo");
    chain2->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ0/output19/dst-19.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output0/dst-0.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output1/dst-1.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output2/dst-2.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output3/dst-3.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output4/dst-4.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output5/dst-5.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output6/dst-6.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output7/dst-7.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output8/dst-8.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output9/dst-9.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output10/dst-10.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output11/dst-11.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output12/dst-12.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output13/dst-13.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output14/dst-14.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output15/dst-15.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output16/dst-16.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output17/dst-17.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output18/dst-18.hipo");
    chain3->AddFile("/w/hallb-scshelf2102/clas12/lixu/singleParticle/testCJ1/output19/dst-19.hipo");


    // Create readers for each chain
    hipo::reader reader1(chain1);
    hipo::reader reader2(chain2);
    hipo::reader reader3(chain3);

     hipo::dictionary factory1;
    reader1.readDictionary(factory1);
    hipo::dictionary factory2;
    reader2.readDictionary(factory2);
    hipo::dictionary factory3;
    reader3.readDictionary(factory3);

    // Get the REC::Particle bank
    hipo::bank recParticle1(factory1.getSchema("REC::Particle"));
    hipo::bank recParticle2(factory2.getSchema("REC::Particle"));
    hipo::bank recParticle3(factory3.getSchema("REC::Particle"));
 
 
    // 1D Histograms
    TH1F *hP1 = new TH1F("hP1", "Neutron Momentum - File 1; p [GeV]; Counts", 100, 0, 5);
    TH1F *hP2 = new TH1F("hP2", "Neutron Momentum - File 2; p [GeV]; Counts", 100, 0, 5);
    TH1F *hP3 = new TH1F("hP3", "Neutron Momentum - File 3; p [GeV]; Counts", 100, 0, 5);

    TH1F *hTheta1 = new TH1F("hTheta1", "Neutron Theta - File 1; #theta [deg]; Counts", 100, 0, 60);
    TH1F *hTheta2 = new TH1F("hTheta2", "Neutron Theta - File 2; #theta [deg]; Counts", 100, 0, 60);
    TH1F *hTheta3 = new TH1F("hTheta3", "Neutron Theta - File 3; #theta [deg]; Counts", 100, 0, 60);

    TH1F *hPhi1 = new TH1F("hPhi1", "Neutron Phi - File 1; #phi [deg]; Counts", 100, -180, 180);
    TH1F *hPhi2 = new TH1F("hPhi2", "Neutron Phi - File 2; #phi [deg]; Counts", 100, -180, 180);
    TH1F *hPhi3 = new TH1F("hPhi3", "Neutron Phi - File 3; #phi [deg]; Counts", 100, -180, 180);

    // 2D Histograms
    TH2F *hPTheta1 = new TH2F("hPTheta1", "Momentum vs Theta - File 1; p [GeV]; #theta [deg]",
                              100, 0, 5, 100, 0, 60);
    TH2F *hPTheta2 = new TH2F("hPTheta2", "Momentum vs Theta - File 2; p [GeV]; #theta [deg]",
                              100, 0, 5, 100, 0, 60);
    TH2F *hPTheta3 = new TH2F("hPTheta3", "Momentum vs Theta - File 3; p [GeV]; #theta [deg]",
                              100, 0, 5, 100, 0, 60);

    TH2F *hPPhi1 = new TH2F("hPPhi1", "Momentum vs Phi - File 1; p [GeV]; #phi [deg]",
                            100, 0, 5, 100, -180, 180);
    TH2F *hPPhi2 = new TH2F("hPPhi2", "Momentum vs Phi - File 2; p [GeV]; #phi [deg]",
                            100, 0, 5, 100, -180, 180);
    TH2F *hPPhi3 = new TH2F("hPPhi3", "Momentum vs Phi - File 3; p [GeV]; #phi [deg]",
                            100, 0, 5, 100, -180, 180);

    TH2F *hThetaPhi1 = new TH2F("hThetaPhi1", "Theta vs Phi - File 1; #theta [deg]; #phi [deg]",
                                100, 0, 60, 100, -180, 180);
    TH2F *hThetaPhi2 = new TH2F("hThetaPhi2", "Theta vs Phi - File 2; #theta [deg]; #phi [deg]",
                                100, 0, 60, 100, -180, 180);
    TH2F *hThetaPhi3 = new TH2F("hThetaPhi3", "Theta vs Phi - File 3; #theta [deg]; #phi [deg]",
                                100, 0, 60, 100, -180, 180);

    // Process file 1
    int entry1 = 0;
    while (reader1.next() && entry1 < 100000)
    {
        entry1++;
        reader1.read(recParticle1);

        int nrows = recParticle1.getRows();
        for (int i = 0; i < nrows; i++)
        {
            int pid = recParticle1.getInt("pid", i);
            if (pid == 2112)
            { // Neutron PDG code
                float px = recParticle1.getFloat("px", i);
                float py = recParticle1.getFloat("py", i);
                float pz = recParticle1.getFloat("pz", i);

                TVector3 p(px, py, pz);
                float momentum = p.Mag();
                float theta = p.Theta() * TMath::RadToDeg();
                float phi = p.Phi() * TMath::RadToDeg();

                // Fill 1D histograms
                hP1->Fill(momentum);
                hTheta1->Fill(theta);
                hPhi1->Fill(phi);

                // Fill 2D histograms
                hPTheta1->Fill(momentum, theta);
                hPPhi1->Fill(momentum, phi);
                hThetaPhi1->Fill(theta, phi);
            }
        }
    }

    // Process file 2
    int entry2 = 0;
    while (reader2.next() && entry2 < 100000)
    {
        entry2++;
        reader2.read(recParticle2);

        int nrows = recParticle2.getRows();
        for (int i = 0; i < nrows; i++)
        {
            int pid = recParticle2.getInt("pid", i);
            if (pid == 2112)
            {
                float px = recParticle2.getFloat("px", i);
                float py = recParticle2.getFloat("py", i);
                float pz = recParticle2.getFloat("pz", i);

                TVector3 p(px, py, pz);
                float momentum = p.Mag();
                float theta = p.Theta() * TMath::RadToDeg();
                float phi = p.Phi() * TMath::RadToDeg();

                // Fill 1D histograms
                hP2->Fill(momentum);
                hTheta2->Fill(theta);
                hPhi2->Fill(phi);

                // Fill 2D histograms
                hPTheta2->Fill(momentum, theta);
                hPPhi2->Fill(momentum, phi);
                hThetaPhi2->Fill(theta, phi);
            }
        }
    }

    // Process file 3
    int entry3 = 0;
    while (reader3.next() && entry3 < 100000)
    {
        entry3++;
        reader3.read(recParticle3);

        int nrows = recParticle3.getRows();
        for (int i = 0; i < nrows; i++)
        {
            int pid = recParticle3.getInt("pid", i);
            if (pid == 2112)
            {
                float px = recParticle3.getFloat("px", i);
                float py = recParticle3.getFloat("py", i);
                float pz = recParticle3.getFloat("pz", i);

                TVector3 p(px, py, pz);
                float momentum = p.Mag();
                float theta = p.Theta() * TMath::RadToDeg();
                float phi = p.Phi() * TMath::RadToDeg();

                // Fill 1D histograms
                hP3->Fill(momentum);
                hTheta3->Fill(theta);
                hPhi3->Fill(phi);

                // Fill 2D histograms
                hPTheta3->Fill(momentum, theta);
                hPPhi3->Fill(momentum, phi);
                hThetaPhi3->Fill(theta, phi);
            }
        }
    }

    // Create canvases and draw comparisons

    // Canvas 1: 1D distributions
    TCanvas *c1 = new TCanvas("c1", "1D Distributions", 1200, 800);
    c1->Divide(3, 1);

    // Momentum comparison
    c1->cd(1);
    hP1->SetLineColor(kRed);
    hP2->SetLineColor(kBlue);
    hP3->SetLineColor(kGreen);
    hP1->Draw();
    hP2->Draw("same");
    hP3->Draw("same");

    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(hP1, "File 1", "l");
    leg1->AddEntry(hP2, "File 2", "l");
    leg1->AddEntry(hP3, "File 3", "l");
    leg1->Draw();

    // Theta comparison
    c1->cd(2);
    hTheta1->SetLineColor(kRed);
    hTheta2->SetLineColor(kBlue);
    hTheta3->SetLineColor(kGreen);
    hTheta1->Draw();
    hTheta2->Draw("same");
    hTheta3->Draw("same");

    // Phi comparison
    c1->cd(3);
    hPhi1->SetLineColor(kRed);
    hPhi2->SetLineColor(kBlue);
    hPhi3->SetLineColor(kGreen);
    hPhi1->Draw();
    hPhi2->Draw("same");
    hPhi3->Draw("same");

    c1->SaveAs("neutron_comparison_1D.png");

    // Canvas 2: 2D momentum vs theta
    TCanvas *c2 = new TCanvas("c2", "Momentum vs Theta", 1200, 400);
    c2->Divide(3, 1);

    c2->cd(1);
    hPTheta1->Draw("colz");
    c2->cd(2);
    hPTheta2->Draw("colz");
    c2->cd(3);
    hPTheta3->Draw("colz");

    c2->SaveAs("neutron_comparison_p_theta.png");

    // Canvas 3: 2D momentum vs phi
    TCanvas *c3 = new TCanvas("c3", "Momentum vs Phi", 1200, 400);
    c3->Divide(3, 1);

    c3->cd(1);
    hPPhi1->Draw("colz");
    c3->cd(2);
    hPPhi2->Draw("colz");
    c3->cd(3);
    hPPhi3->Draw("colz");

    c3->SaveAs("neutron_comparison_p_phi.png");

    // Canvas 4: 2D theta vs phi
    TCanvas *c4 = new TCanvas("c4", "Theta vs Phi", 1200, 400);
    c4->Divide(3, 1);

    c4->cd(1);
    hThetaPhi1->Draw("colz");
    c4->cd(2);
    hThetaPhi2->Draw("colz");
    c4->cd(3);
    hThetaPhi3->Draw("colz");

    c4->SaveAs("neutron_comparison_theta_phi.png");

    // Print statistics
    cout << "Processed " << entry1 << " events from file 1" << endl;
    cout << "Processed " << entry2 << " events from file 2" << endl;
    cout << "Processed " << entry3 << " events from file 3" << endl;
}
