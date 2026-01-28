// To run: clas12root -l -b -q 'compare_cnd_versions.C+(100000)'

#include <TROOT.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <string>
#include <cmath>
#include <TVector2.h>
#include <limits>
#include <algorithm>

#include <TFile.h>
#include <cstdlib>
#include <chrono>
#include <TTree.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TBenchmark.h>
// #include "reader.h"

#include "hipo4/reader.h"
#include "hipo4/dictionary.h"
#include "hipo4/bank.h"
#include "hipo4/event.h"

#include <sys/stat.h>
#include <sys/types.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

void ensure_dir(const char *dirname)
{
    struct stat st;
    if (stat(dirname, &st) != 0)
    {
        std::cout << "Creating directory: " << dirname << std::endl;
        mkdir(dirname, 0755);
    }
}

TChain *make_chain(const char *baseDir,
                   const char *testFolder,
                   int first = 0,
                   int last = 19)
{
    auto *ch = new TChain("hipoFiles");
    for (int i = first; i <= last; ++i)
    {
        TString path = TString::Format("%s/%s/output%d/dst-%d.hipo", baseDir, testFolder, i, i);
        ch->Add(path);
    }
    return ch;
}

struct Hists
{
    TH1F *hP;
    TH1F *hTheta;
    TH1F *hPhi;

    TH1F *hEnergy_CND;
    TH1F *hTheta_CND;
    TH1F *hPhi_CND;

    TH1F *hDTheta_CND;
    TH1F *hDPhi_CND;

    TH2F *hPTheta;
    TH2F *hPPhi;
    TH2F *hThetaPhi;
};

Hists book_hists(const char *tag)
{
    Hists h;
    h.hP = new TH1F(TString::Format("hP_%s", tag), TString::Format("Neutron p (%s) from REC::Particle;p [GeV];Counts", tag), 100, 0, 5);
    h.hTheta = new TH1F(TString::Format("hTh_%s", tag), TString::Format("Neutron #theta (%s) from REC::Particle;#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhi = new TH1F(TString::Format("hPhi_%s", tag), TString::Format("Neutron #phi (%s) from REC::Particle;#phi [deg];Counts", tag), 100, -180, 180);

    h.hEnergy_CND = new TH1F(TString::Format("hECND_%s", tag), TString::Format("Neutron E (%s) from CND-hit;E [GeV];Counts", tag), 100, 0, 5);
    h.hTheta_CND = new TH1F(TString::Format("hThCND_%s", tag), TString::Format("Neutron #theta (%s) from CND-hit;#theta [deg];Counts", tag), 100, 0, 60);
    h.hPhi_CND = new TH1F(TString::Format("hPhiCND_%s", tag), TString::Format("Neutron #phi (%s) from CND-hit;#phi [deg];Counts", tag), 100, -180, 180);

    h.hDTheta_CND = new TH1F(TString::Format("hDThCND_%s", tag), "#Delta#theta(CND-hit - particle);#Delta#theta [deg];Counts", 120, -30, 30);
    h.hDPhi_CND = new TH1F(TString::Format("hDPhCND_%s", tag), "#Delta#phi(CND-hit - particle);#Delta#phi [deg];Counts", 180, -180, 180);

    h.hPTheta = new TH2F(TString::Format("hPTh_%s", tag), TString::Format("p vs #theta (%s) from REC::Particle;p [GeV];#theta [deg]", tag),
                         100, 0, 5, 100, 0, 180);
    h.hPPhi = new TH2F(TString::Format("hPPhi_%s", tag), TString::Format("p vs #phi (%s) from REC::Particle;p [GeV];#phi [deg]", tag),
                       100, 0, 5, 100, -180, 180);
    h.hThetaPhi = new TH2F(TString::Format("hThPhi_%s", tag), TString::Format("#theta vs #phi (%s) from REC::Particle;#theta [deg];#phi [deg]", tag),
                           100, 0, 180, 100, -180, 180);
    return h;
}

void process_chain(TChain *chain, Hists &h, const char *tag, int maxEvents = 300000, int cnd_id = 3)
{
    if (!chain)
    {
        std::cerr << "Null chain for " << tag << std::endl;
        return;
    }

    TObjArray *files = chain->GetListOfFiles();
    if (!files || files->GetLast() < 0)
    {
        std::cerr << "No files in chain for " << tag << std::endl;
        return;
    }

    const double Nmass = TDatabasePDG::Instance()->GetParticle(2112)->Mass();

    int processed = 0;

    hipo::reader reader;
    hipo::dictionary factory;
    hipo::event event;

    TLorentzVector N_Vec_temp;
    // vector<TLorentzVector> N_Vec_;
    // vector<double> N_info_temp;
    // vector<vector<double>> N_info_;

    for (int fi = 0; fi <= files->GetLast(); ++fi)
    {
        const char *fname = files->At(fi)->GetTitle();

        reader.open(fname);
        reader.readDictionary(factory);

        // Check required banks
        if (!factory.hasSchema("REC::Particle") ||
            !factory.hasSchema("REC::Scintillator") ||
            !factory.hasSchema("RUN::config"))
        {
            std::cerr << "Missing required schema in " << fname
                      << " for " << tag << std::endl;
            continue;
        }

        hipo::bank CONF(factory.getSchema("RUN::config"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank SCINT(factory.getSchema("REC::Scintillator"));

        // bool hasScintX = factory.hasSchema("REC::ScintExtras");
        // hipo::bank SCINTX(hasScintX ? factory.getSchema("REC::ScintExtras")
        //                             : factory.getSchema("REC::Scintillator"));

        while (reader.next())
        {
            reader.read(event);

            event.getStructure(CONF);
            event.getStructure(PART);
            event.getStructure(SCINT);
            // if (hasScintX) event.getStructure(SCINTX);

            int RunNumber = CONF.getInt("run", 0);
            int EventNumber = CONF.getInt("event", 0);

            const int nPart = PART.getRows();
            const int nSc = SCINT.getRows();
            // const int nScX = hasScintX ? SCINTX.getRows() : 0;

            // // Build a fast lookup for ScintExtras by "index" if possible
            // std::map<int, int> scintX_by_index;
            // if (hasScintX && nScX > 0)
            // {
            //     // Many CLAS12 banks use "index" to link related rows.
            //     // If SCINTX does not have "index", this will need adjustment.
            //     for (int ix = 0; ix < nScX; ++ix)
            //     {
            //         if (SCINTX.getSchema().hasEntry("index"))
            //         {
            //             scintX_by_index[SCINTX.getInt("index", ix)] = ix;
            //         }
            //     }
            // }

            std::vector<int> bestScRow(nPart, -1);
            std::vector<float> bestE(nPart, -1.0f);

            // std::vector<double> part_Scint_CND_E(nPart, NAN);
            // std::vector<double> part_Scint_CND_t(nPart, NAN);
            // std::vector<double> part_Scint_CND_x(nPart, NAN);
            // std::vector<double> part_Scint_CND_y(nPart, NAN);
            // std::vector<double> part_Scint_CND_z(nPart, NAN);
            // std::vector<double> part_ScintX_CND_dedx(nPart, NAN);
            // std::vector<double> part_ScintX_CND_size(nPart, NAN);
            // std::vector<double> part_ScintX_CND_layermult(nPart, NAN);

            // Fill per-particle CND info from scintillator rows (detector == CND)
            for (int sc = 0; sc < nSc; ++sc)
            {
                const int detector = SCINT.getInt("detector", sc);
                if (detector != cnd_id)
                    continue;

                const int pindex = SCINT.getInt("pindex", sc);
                if (pindex < 0 || pindex >= nPart)
                    continue;

                const float E = SCINT.getFloat("energy", sc);
                if (E > bestE[pindex])
                {
                    bestE[pindex] = E;
                    bestScRow[pindex] = sc;
                }

                // Found a CND-related scintillator hit
                // part_Scint_CND_E[pindex] = SCINT.getFloat("energy", sc);
                // part_Scint_CND_t[pindex] = SCINT.getFloat("time", sc);
                // part_Scint_CND_x[pindex] = SCINT.getFloat("x", sc);
                // part_Scint_CND_y[pindex] = SCINT.getFloat("y", sc);
                // part_Scint_CND_z[pindex] = SCINT.getFloat("z", sc);

                // Extras: try to match by "index" if both have it
                // if (hasScintX && SCINT.getSchema().hasEntry("index") && SCINTX.getSchema().hasEntry("index"))
                // {
                //     const int idx = SCINT.getInt("index", sc);
                //     auto it = scintX_by_index.find(idx);
                //     if (it != scintX_by_index.end())
                //     {
                //         const int ix = it->second;
                //         if (SCINTX.getSchema().hasEntry("dedx"))
                //             cndDedx[pindex] = SCINTX.getFloat("dedx", ix);
                //         if (SCINTX.getSchema().hasEntry("size"))
                //             cndSize[pindex] = SCINTX.getInt("size", ix);
                //         if (SCINTX.getSchema().hasEntry("layermult"))
                //             cndLayermult[pindex] = SCINTX.getInt("layermult", ix);
                //     }
                // }
            }

            // Loop particles and fill neutron hists
            for (int ip = 0; ip < nPart; ++ip)
            {

                const int pid = PART.getInt("pid", ip);
                if (pid != 2112)
                    continue;

                const int pcharge = PART.getInt("charge", ip);
                const int psatus = PART.getInt("status", ip);

                const float px = PART.getFloat("px", ip);
                const float py = PART.getFloat("py", ip);
                const float pz = PART.getFloat("pz", ip);

                const double pvx = PART.getFloat("vx", ip);
                const double pvy = PART.getFloat("vy", ip);
                const double pvz = PART.getFloat("vz", ip);

                const double pchi2pid = PART.getFloat("chi2pid", ip);
                const double pbeta = PART.getFloat("beta", ip);

                TVector3 p3(px, py, pz);
                const float p = p3.Mag();
                const float theta = p3.Theta() * TMath::RadToDeg();
                const float phi = p3.Phi() * TMath::RadToDeg();
                N_Vec_temp.SetPxPyPzE(px, py, pz, TMath::Sqrt(p * p + Nmass * Nmass));

                h.hP->Fill(p);
                h.hTheta->Fill(theta);
                h.hPhi->Fill(phi);
                h.hPTheta->Fill(p, theta);
                h.hPPhi->Fill(p, phi);
                h.hThetaPhi->Fill(theta, phi);

                // N_Vec_.push_back(N_Vec_temp);
                // N_info_temp.push_back(pvx);
                // N_info_temp.push_back(pvy);
                // N_info_temp.push_back(pvz);
                // N_info_temp.push_back(psatus);
                // N_info_temp.push_back(pchi2pid);
                // N_info_temp.push_back(pbeta);

                // N_info_temp.push_back(part_Scint_CND_E[ip]);
                // N_info_temp.push_back(part_Scint_CND_t[ip]);
                // N_info_temp.push_back(part_Scint_CND_x[ip]);
                // N_info_temp.push_back(part_Scint_CND_y[ip]);
                // N_info_temp.push_back(part_Scint_CND_z[ip]);

                // N_info_temp.push_back(part_ScintX_CND_dedx[ip]);
                // N_info_temp.push_back(part_ScintX_CND_size[ip]);
                // N_info_temp.push_back(part_ScintX_CND_layermult[ip]);

                // N_info_.push_back(N_info_temp);

                // Now fill CND-related info if we have an associated row
                const int scBest = bestScRow[ip];
                if (scBest < 0)
                    continue;

                const float E = SCINT.getFloat("energy", scBest);
                const float x = SCINT.getFloat("x", scBest);
                const float y = SCINT.getFloat("y", scBest);
                const float z = SCINT.getFloat("z", scBest);

                TVector3 r(x, y, z);
                const float th_hit = r.Theta() * TMath::RadToDeg();
                const float ph_hit = r.Phi() * TMath::RadToDeg();

                // Delta phi with wrapping
                const float dth = th_hit - theta;
                const float dph = TVector2::Phi_mpi_pi((ph_hit - phi) * TMath::DegToRad()) * TMath::RadToDeg();

                h.hEnergy_CND->Fill(E);
                h.hTheta_CND->Fill(th_hit);
                h.hPhi_CND->Fill(ph_hit);
                h.hDTheta_CND->Fill(dth);
                h.hDPhi_CND->Fill(dph);
            }

            processed++;
            if (processed >= maxEvents)
                break;
        }

        if (processed >= maxEvents)
            break;
    }

    std::cout << "Processed " << processed << " events for " << tag << std::endl;
}

// Helper: style and overlay 3 TH1s (optionally normalized)
void draw_overlay_1D(TH1 *h1_in, TH1 *h2_in, TH1 *h3_in,
                     const char *leg1, const char *leg2, const char *leg3,
                     bool normalize = true)
{
    if (!h1_in || !h2_in || !h3_in) return;

    gPad->SetLogy();

    // Clone so normalization does not permanently change original histograms
    auto *h1 = (TH1*)h1_in->Clone(Form("%s_clone1", h1_in->GetName()));
    auto *h2 = (TH1*)h2_in->Clone(Form("%s_clone2", h2_in->GetName()));
    auto *h3 = (TH1*)h3_in->Clone(Form("%s_clone3", h3_in->GetName()));

    h1->SetDirectory(nullptr);
    h2->SetDirectory(nullptr);
    h3->SetDirectory(nullptr);

    if (normalize)
    {
        if (h1->Integral() > 0) h1->Scale(1.0 / h1->Integral());
        if (h2->Integral() > 0) h2->Scale(1.0 / h2->Integral());
        if (h3->Integral() > 0) h3->Scale(1.0 / h3->Integral());
        h1->GetYaxis()->SetTitle("Normalized counts");
    }

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kGreen + 2);

    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h3->SetLineWidth(2);

    double m = std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
    h1->SetMaximum(1.15 * m);

    h1->Draw("hist");
    h2->Draw("hist same");
    h3->Draw("hist same");

    auto *leg = new TLegend(0.30, 0.72, 0.55, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h1, leg1, "l");
    leg->AddEntry(h2, leg2, "l");
    leg->AddEntry(h3, leg3, "l");
    leg->Draw();
}

// Helper: draw 3 TH2 side-by-side
void draw_triptych_2D(TH2 *h1, TH2 *h2, TH2 *h3,
                      const char *title,
                      const char *outname)
{
    if (!h1 || !h2 || !h3)
        return;

    auto *c = new TCanvas(Form("c2D_%s", outname), title, 1500, 450);
    c->Divide(3, 1);

    c->cd(1);
    gPad->SetRightMargin(0.14);
    h1->Draw("colz");
    c->cd(2);
    gPad->SetRightMargin(0.14);
    h2->Draw("colz");
    c->cd(3);
    gPad->SetRightMargin(0.14);
    h3->Draw("colz");

    c->SaveAs(outname);
}

void compare_cnd_versions(int maxEvents = 300000)
{
    // gSystem->Load("libhipo4");

    const char *BASE = "/ceph24/hallb/clas12/users/lixu/singleParticle";

    const char *OUTDIR = "Histograms";
    ensure_dir(OUTDIR);

    // Build chains
    TChain *chOSG = make_chain(BASE, "testOSG");
    TChain *chCJ0 = make_chain(BASE, "testCJ0");
    TChain *chCJ1 = make_chain(BASE, "testCJ1");

    // Book hists
    Hists hOSG = book_hists("OSG");
    Hists hCJ0 = book_hists("CJ0");
    Hists hCJ1 = book_hists("CJ1");

    // Fill
    process_chain(chOSG, hOSG, "OSG", maxEvents);
    process_chain(chCJ0, hCJ0, "CJ0", maxEvents);
    process_chain(chCJ1, hCJ1, "CJ1", maxEvents);

    // ------------------------------------------------------------
    // 1) REC::Particle neutrons: p, theta, phi
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_particle_1D", "REC::Particle neutrons", 1500, 450);
        c->Divide(3, 1);
        
        c->cd(1);
        draw_overlay_1D(hOSG.hP, hCJ0.hP, hCJ1.hP, "OSG", "CJ0", "CJ1", false);
        c->cd(2);
        draw_overlay_1D(hOSG.hTheta, hCJ0.hTheta, hCJ1.hTheta, "OSG", "CJ0", "CJ1", false);
        c->cd(3);
        draw_overlay_1D(hOSG.hPhi, hCJ0.hPhi, hCJ1.hPhi, "OSG", "CJ0", "CJ1", false);

        c->SaveAs("Histograms/cmp_RECParticle_neutrons_1D.png");
    }

    // ------------------------------------------------------------
    // 2) CND-matched scintillator info: energy, theta_hit, phi_hit
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_cndhit_1D", "CND matched scintillator (REC::Scintillator, detector==CND)", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D(hOSG.hEnergy_CND, hCJ0.hEnergy_CND, hCJ1.hEnergy_CND, "OSG", "CJ0", "CJ1", false);
        c->cd(2);
        draw_overlay_1D(hOSG.hTheta_CND, hCJ0.hTheta_CND, hCJ1.hTheta_CND, "OSG", "CJ0", "CJ1", false);
        c->cd(3);
        draw_overlay_1D(hOSG.hPhi_CND, hCJ0.hPhi_CND, hCJ1.hPhi_CND, "OSG", "CJ0", "CJ1", false);
        c->SaveAs("Histograms/cmp_CND_scint_1D.png");
    }

    // ------------------------------------------------------------
    // 3) Residuals: Δtheta, Δphi  (hit - particle)
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_residuals_1D", "CND residuals", 1200, 450);
        c->Divide(2, 1);

        c->cd(1);
        draw_overlay_1D(hOSG.hDTheta_CND, hCJ0.hDTheta_CND, hCJ1.hDTheta_CND, "OSG", "CJ0", "CJ1", false);
        c->cd(2);
        draw_overlay_1D(hOSG.hDPhi_CND, hCJ0.hDPhi_CND, hCJ1.hDPhi_CND, "OSG", "CJ0", "CJ1", false);

        c->SaveAs("Histograms/cmp_CND_residuals_1D.png");
    }

    // ------------------------------------------------------------
    // 4) 2D histograms (triptychs)
    // ------------------------------------------------------------
    draw_triptych_2D(hOSG.hPTheta, hCJ0.hPTheta, hCJ1.hPTheta,
                     "p vs theta (REC::Particle)", "Histograms/cmp_2D_p_vs_theta.png");

    draw_triptych_2D(hOSG.hPPhi, hCJ0.hPPhi, hCJ1.hPPhi,
                     "p vs phi (REC::Particle)", "Histograms/cmp_2D_p_vs_phi.png");

    draw_triptych_2D(hOSG.hThetaPhi, hCJ0.hThetaPhi, hCJ1.hThetaPhi,
                     "theta vs phi (REC::Particle)", "Histograms/cmp_2D_theta_vs_phi.png");

    // ------------------------------------------------------------
    // 5) Save all histograms into a ROOT file
    // ------------------------------------------------------------
    {
        TFile *fout = TFile::Open("Histograms/cnd_comparison.root", "RECREATE");
        if (!fout || fout->IsZombie())
        {
            std::cerr << "Error opening output ROOT file" << std::endl;
            return;
        }

        // OSG
        hOSG.hP->Write("", TObject::kOverwrite);
        hOSG.hTheta->Write("", TObject::kOverwrite);
        hOSG.hPhi->Write("", TObject::kOverwrite);
        hOSG.hEnergy_CND->Write("", TObject::kOverwrite);
        hOSG.hTheta_CND->Write("", TObject::kOverwrite);
        hOSG.hPhi_CND->Write("", TObject::kOverwrite);
        hOSG.hDTheta_CND->Write("", TObject::kOverwrite);
        hOSG.hDPhi_CND->Write("", TObject::kOverwrite);
        hOSG.hPTheta->Write("", TObject::kOverwrite);
        hOSG.hPPhi->Write("", TObject::kOverwrite);
        hOSG.hThetaPhi->Write("", TObject::kOverwrite);

        // CJ0
        hCJ0.hP->Write("", TObject::kOverwrite);
        hCJ0.hTheta->Write("", TObject::kOverwrite);
        hCJ0.hPhi->Write("", TObject::kOverwrite);
        hCJ0.hEnergy_CND->Write("", TObject::kOverwrite);
        hCJ0.hTheta_CND->Write("", TObject::kOverwrite);
        hCJ0.hPhi_CND->Write("", TObject::kOverwrite);
        hCJ0.hDTheta_CND->Write("", TObject::kOverwrite);
        hCJ0.hDPhi_CND->Write("", TObject::kOverwrite);
        hCJ0.hPTheta->Write("", TObject::kOverwrite);
        hCJ0.hPPhi->Write("", TObject::kOverwrite);
        hCJ0.hThetaPhi->Write("", TObject::kOverwrite);

        // CJ1
        hCJ1.hP->Write("", TObject::kOverwrite);
        hCJ1.hTheta->Write("", TObject::kOverwrite);
        hCJ1.hPhi->Write("", TObject::kOverwrite);
        hCJ1.hEnergy_CND->Write("", TObject::kOverwrite);
        hCJ1.hTheta_CND->Write("", TObject::kOverwrite);
        hCJ1.hPhi_CND->Write("", TObject::kOverwrite);
        hCJ1.hDTheta_CND->Write("", TObject::kOverwrite);
        hCJ1.hDPhi_CND->Write("", TObject::kOverwrite);
        hCJ1.hPTheta->Write("", TObject::kOverwrite);
        hCJ1.hPPhi->Write("", TObject::kOverwrite);
        hCJ1.hThetaPhi->Write("", TObject::kOverwrite);

        fout->Close();
    }
    std::cout << "Done. Wrote PNGs + cmp_cnd_versions_hists.root" << std::endl;
}