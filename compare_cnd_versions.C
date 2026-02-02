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

    TH1F *hPMC;
    TH1F *hThetaMC;
    TH1F *hPhiMC;

    TH1F *hDTheta_REC_MC;
    TH1F *hDPhi_REC_MC;
    TH1F *hAngle_REC_MC; // opening angle between REC and MC neutrons (deg)
    TH1F *hDP_REC_MC;    // delta p between REC and MC neutrons (p_rec - p_mc [GeV])
};

Hists book_hists(const char *tag)
{
    Hists h;
    h.hP = new TH1F(TString::Format("hP_%s", tag), TString::Format("Neutron p (%s) from REC::Particle;p [GeV];Counts", tag), 100, 0, 5);
    h.hTheta = new TH1F(TString::Format("hTh_%s", tag), TString::Format("Neutron #theta (%s) from REC::Particle;#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhi = new TH1F(TString::Format("hPhi_%s", tag), TString::Format("Neutron #phi (%s) from REC::Particle;#phi [deg];Counts", tag), 100, -180, 180);

    h.hEnergy_CND = new TH1F(TString::Format("hECND_%s", tag), TString::Format("Neutron E (%s) from CND-hit;E [GeV];Counts", tag), 100, 0, 12);
    h.hTheta_CND = new TH1F(TString::Format("hThCND_%s", tag), TString::Format("Neutron #theta (%s) from CND-hit;#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhi_CND = new TH1F(TString::Format("hPhiCND_%s", tag), TString::Format("Neutron #phi (%s) from CND-hit;#phi [deg];Counts", tag), 100, -180, 180);

    h.hDTheta_CND = new TH1F(TString::Format("hDThCND_%s", tag), "#Delta#theta(CND-hit - particle);#Delta#theta [deg];Counts", 100, -30, 30);
    h.hDPhi_CND = new TH1F(TString::Format("hDPhCND_%s", tag), "#Delta#phi(CND-hit - particle);#Delta#phi [deg];Counts", 100, -180, 180);

    h.hPTheta = new TH2F(TString::Format("hPTh_%s", tag), TString::Format("#theta vs p (%s) from REC::Particle;p [GeV];#theta [deg]", tag),
                         100, 0, 10, 100, 0, 180);
    h.hPPhi = new TH2F(TString::Format("hPPhi_%s", tag), TString::Format("#phi vs p (%s) from REC::Particle;p [GeV];#phi [deg]", tag),
                       100, 0, 10, 100, -180, 180);
    h.hThetaPhi = new TH2F(TString::Format("hThPhi_%s", tag), TString::Format("#phi vs #theta (%s) from REC::Particle;#theta [deg];#phi [deg]", tag),
                           100, 0, 180, 100, -180, 180);
                           
    h.hPMC = new TH1F(TString::Format("hPMC_%s", tag), TString::Format("MC Neutron p (%s);p [GeV];Counts", tag), 100, 0, 5);
    h.hThetaMC = new TH1F(TString::Format("hThMC_%s", tag), TString::Format("MC Neutron #theta (%s);#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhiMC = new TH1F(TString::Format("hPhiMC_%s", tag), TString::Format("MC Neutron #phi (%s);#phi [deg];Counts", tag), 100, -180, 180);

    h.hDTheta_REC_MC = new TH1F(TString::Format("hDTh_REC_MC_%s", tag), "#Delta#theta(REC - MC);#Delta#theta [deg];Counts", 200, -12, 12);
    h.hDPhi_REC_MC = new TH1F(TString::Format("hDPh_REC_MC_%s", tag), "#Delta#phi(REC - MC);#Delta#phi [deg];Counts", 200, -15, 15);
    h.hAngle_REC_MC = new TH1F(TString::Format("hAngle_REC_MC_%s", tag), "Opening angle(REC - MC);Opening angle [deg];Counts", 200, 0, 180);
    h.hDP_REC_MC = new TH1F(TString::Format("hDP_REC_MC_%s", tag), "#Delta p(REC) - p(MC);#Delta p [GeV];Counts", 200, -1.5, 1.0);
    return h;
}

void write_hists(const Hists &h)
{
    std::vector<TObject *> objs = {
        h.hP, h.hTheta, h.hPhi,
        h.hEnergy_CND, h.hTheta_CND, h.hPhi_CND,
        h.hDTheta_CND, h.hDPhi_CND,
        h.hPTheta, h.hPPhi, h.hThetaPhi,
        h.hPMC, h.hThetaMC, h.hPhiMC,
        h.hDTheta_REC_MC, h.hDPhi_REC_MC, h.hAngle_REC_MC, h.hDP_REC_MC};

    for (auto *o : objs)
    {
        if (!o)
            continue;
        o->Write("", TObject::kOverwrite);
    }
}

struct MatchResult
{
    int mcIndex = -1;
    float angleDeg = 1e9;
};

MatchResult match_neutron_rec_to_mc(const TVector3 &pREC,
                                    const hipo::bank &MCPT,
                                    float maxAngleDeg = 10.0)
{
    MatchResult out;

    const int nMC = MCPT.getRows();
    for (int imc = 0; imc < nMC; ++imc)
    {
        const int pid = MCPT.getInt("pid", imc);
        if (pid != 2112)
            continue; // neutron

        const float px = MCPT.getFloat("px", imc);
        const float py = MCPT.getFloat("py", imc);
        const float pz = MCPT.getFloat("pz", imc);

        TVector3 pMC(px, py, pz);
        if (pMC.Mag() < 1e-6)
            continue; // skip zero-momentum

        const float ang = pREC.Angle(pMC) * TMath::RadToDeg();

        if (ang < out.angleDeg)
        {
            out.angleDeg = ang;
            out.mcIndex = imc;
        }
    }

    if (out.mcIndex >= 0 && out.angleDeg > maxAngleDeg)
    {
        // No good match found
        out.mcIndex = -1;
    }
    return out;
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
    const float pt_min = 0.05; // minimum neutron transverse momentum [GeV]

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
            !factory.hasSchema("RUN::config") ||
            !factory.hasSchema("MC::Particle"))
        {
            std::cerr << "Missing required schema in " << fname
                      << " for " << tag << std::endl;
            continue;
        }

        hipo::bank CONF(factory.getSchema("RUN::config"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank SCINT(factory.getSchema("REC::Scintillator"));
        hipo::bank MCPT(factory.getSchema("MC::Particle"));

        // bool hasScintX = factory.hasSchema("REC::ScintExtras");
        // hipo::bank SCINTX(hasScintX ? factory.getSchema("REC::ScintExtras")
        //                             : factory.getSchema("REC::Scintillator"));

        while (reader.next())
        {
            reader.read(event);

            event.getStructure(CONF);
            event.getStructure(PART);
            event.getStructure(SCINT);
            event.getStructure(MCPT);
            // if (hasScintX) event.getStructure(SCINTX);

            int RunNumber = CONF.getInt("run", 0);
            int EventNumber = CONF.getInt("event", 0);

            const int nPart = PART.getRows();
            const int nSc = SCINT.getRows();
            const int nMC = MCPT.getRows();
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
                const float t = SCINT.getFloat("time", sc);
                const float x = SCINT.getFloat("x", sc);
                const float y = SCINT.getFloat("y", sc);
                const float z = SCINT.getFloat("z", sc);

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

                TVector3 pREC(px, py, pz);
                const float pt_rec = std::sqrt(px * px + py * py);

                if (pt_rec < pt_min)
                    continue; // too low momentum

                MatchResult match = match_neutron_rec_to_mc(pREC, MCPT, /*maxAngleDeg=*/10.0);
                if (match.mcIndex < 0)
                    continue; // no good MC match

                const float p_rec = pREC.Mag();
                const float th_rec = pREC.Theta() * TMath::RadToDeg();
                const float ph_rec = pREC.Phi() * TMath::RadToDeg();

                if (std::abs(ph_rec) < 1e-6)
                {
                    std::cout << "phi=0: px=" << px << " py=" << py << " pz=" << pz << "\n";
                }

                N_Vec_temp.SetPxPyPzE(px, py, pz, TMath::Sqrt(p_rec * p_rec + Nmass * Nmass));

                h.hP->Fill(p_rec);
                h.hTheta->Fill(th_rec);
                h.hPhi->Fill(ph_rec);
                h.hPTheta->Fill(p_rec, th_rec);
                h.hPPhi->Fill(p_rec, ph_rec);
                h.hThetaPhi->Fill(th_rec, ph_rec);

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

                const float mc_px = MCPT.getFloat("px", match.mcIndex);
                const float mc_py = MCPT.getFloat("py", match.mcIndex);
                const float mc_pz = MCPT.getFloat("pz", match.mcIndex);

                TVector3 pMC(mc_px, mc_py, mc_pz);

                const float th_mc = pMC.Theta() * TMath::RadToDeg();
                const float ph_mc = pMC.Phi() * TMath::RadToDeg();

                const float dth_rec_mc = th_rec - th_mc;
                const float dph_rec_mc = TVector2::Phi_mpi_pi((ph_rec - ph_mc) * TMath::DegToRad()) * TMath::RadToDeg();

                h.hDTheta_REC_MC->Fill(dth_rec_mc);
                h.hDPhi_REC_MC->Fill(dph_rec_mc);
                h.hAngle_REC_MC->Fill(match.angleDeg);
                h.hDP_REC_MC->Fill(p_rec - pMC.Mag());

                h.hPMC->Fill(pMC.Mag());
                h.hThetaMC->Fill(th_mc);
                h.hPhiMC->Fill(ph_mc);

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
                const float dth = th_hit - th_rec;
                const float dph = TVector2::Phi_mpi_pi((ph_hit - ph_rec) * TMath::DegToRad()) * TMath::RadToDeg();

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

void draw_overlay_1D_N(const std::vector<TH1 *> &h_in,
                       const std::vector<const char *> &labels,
                       bool normalize = true,
                       bool logy = true,
                       double leg_x1 = 0.15, double leg_y1 = 0.72,
                       double leg_x2 = 0.35, double leg_y2 = 0.90,
                       int leg_ncols = 1)
{
    if (h_in.empty() || h_in.size() != labels.size())
        return;

    // Require all non-null
    for (auto *h : h_in)
        if (!h)
            return;

    if (logy)
        gPad->SetLogy();

    // Clone to avoid modifying originals
    std::vector<TH1 *> h;
    h.reserve(h_in.size());
    for (size_t i = 0; i < h_in.size(); ++i)
    {
        auto *c = (TH1 *)h_in[i]->Clone(Form("%s_clone_%zu", h_in[i]->GetName(), i));
        c->SetDirectory(nullptr);
        h.push_back(c);
    }

    // Normalize
    if (normalize)
    {
        for (auto *hi : h)
        {
            const double I = hi->Integral(0, hi->GetNbinsX() + 1);
            if (I > 0)
                hi->Scale(1.0 / I);
        }
        h[0]->GetYaxis()->SetTitle("Normalized counts");
    }

    // Colors (edit as you like)
    const std::vector<int> colors = {
        kRed, kBlue, kGreen + 2, kMagenta, kOrange - 3, kCyan,
        kViolet + 1, kGray + 2};

    for (size_t i = 0; i < h.size(); ++i)
    {
        h[i]->SetLineColor(colors[i % colors.size()]);
        h[i]->SetLineWidth(2);
    }

    // Max for axis range
    double m = 0.0;
    for (auto *hi : h)
        m = std::max(m, (double)hi->GetMaximum());

    if (logy)
        h[0]->SetMaximum(10.0 * m);
    else
        h[0]->SetMaximum(2.5 * m);

    // Draw
    h[0]->Draw("hist");
    for (size_t i = 1; i < h.size(); ++i)
        h[i]->Draw("hist same");

    // Legend
    auto *leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    leg->SetNColumns(leg_ncols);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (size_t i = 0; i < h.size(); ++i)
        leg->AddEntry(h[i], labels[i], "f"); 

    leg->Draw();
}

// Helper: draw 3 TH2 side-by-side
void draw_triptych_2D(TH2 *h1_in, TH2 *h2_in, TH2 *h3_in,
                      const char *title,
                      const char *outname,
                      bool normalize = true,
                      bool logz = true)
{
    if (!h1_in || !h2_in || !h3_in)
        return;

    // Clone to avoid modifying originals
    auto *h1 = (TH2 *)h1_in->Clone(Form("%s_clone1", h1_in->GetName()));
    auto *h2 = (TH2 *)h2_in->Clone(Form("%s_clone2", h2_in->GetName()));
    auto *h3 = (TH2 *)h3_in->Clone(Form("%s_clone3", h3_in->GetName()));

    h1->SetDirectory(nullptr);
    h2->SetDirectory(nullptr);
    h3->SetDirectory(nullptr);

    // Optional normalization (global integral)
    if (normalize)
    {
        if (h1->Integral() > 0)
            h1->Scale(1.0 / h1->Integral());
        if (h2->Integral() > 0)
            h2->Scale(1.0 / h2->Integral());
        if (h3->Integral() > 0)
            h3->Scale(1.0 / h3->Integral());
    }

    gPad->Modified();
    gPad->Update();

    auto *c = new TCanvas(Form("c2D_%s", outname), title, 1800, 450);
    c->Divide(3, 1);

    c->cd(1);
    gPad->SetRightMargin(0.14);
    if (logz)
        gPad->SetLogz();
    h1->SetStats(0);
    h1->Draw("colz");

    c->cd(2);
    gPad->SetRightMargin(0.14);
    if (logz)
        gPad->SetLogz();
    h2->SetStats(0);
    h2->Draw("colz");

    c->cd(3);
    gPad->SetRightMargin(0.14);
    if (logz)
        gPad->SetLogz();
    h3->SetStats(0);
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
        auto *c = new TCanvas("c_particle_1D", "REC::Particle neutrons and MC::Particle", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hP, hCJ0.hP, hCJ1.hP, hOSG.hPMC, hCJ0.hPMC, hCJ1.hPMC}, {"OSG REC::Particle", "CJ0 REC::Particle", "CJ1 REC::Particle", "OSG MC::Particle", "CJ0 MC::Particle", "CJ1 MC::Particle"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.12,0.75,0.83,0.90, /*ncols=*/2);
        c->cd(2);
        draw_overlay_1D_N({hOSG.hTheta, hCJ0.hTheta, hCJ1.hTheta, hOSG.hThetaMC, hCJ0.hThetaMC, hCJ1.hThetaMC}, {"OSG REC::Particle", "CJ0 REC::Particle", "CJ1 REC::Particle", "OSG MC::Particle", "CJ0 MC::Particle", "CJ1 MC::Particle"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.12,0.75,0.83,0.90, /*ncols=*/2);
        c->cd(3);
        draw_overlay_1D_N({hOSG.hPhi, hCJ0.hPhi, hCJ1.hPhi, hOSG.hPhiMC, hCJ0.hPhiMC, hCJ1.hPhiMC}, {"OSG REC::Particle", "CJ0 REC::Particle", "CJ1 REC::Particle", "OSG MC::Particle", "CJ0 MC::Particle", "CJ1 MC::Particle"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.12,0.75,0.83,0.90, /*ncols=*/2);
        c->SaveAs("Histograms/cmp_RECParticle_neutrons_1D.png");
    }

    // ------------------------------------------------------------
    // 2) CND-matched scintillator info: energy, theta_hit, phi_hit
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_cndhit_1D", "CND matched scintillator (REC::Scintillator, detector==CND)", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hEnergy_CND, hCJ0.hEnergy_CND, hCJ1.hEnergy_CND}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->cd(2);
        draw_overlay_1D_N({hOSG.hTheta_CND, hCJ0.hTheta_CND, hCJ1.hTheta_CND}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->cd(3);
        draw_overlay_1D_N({hOSG.hPhi_CND, hCJ0.hPhi_CND, hCJ1.hPhi_CND}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->SaveAs("Histograms/cmp_CND_scint_1D.png");
    }

    // ------------------------------------------------------------
    // 3) Residuals: Δtheta, Δphi  (hit - particle)
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_residuals_1D", "CND residuals", 1500, 450);
        c->Divide(2, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hDTheta_CND, hCJ0.hDTheta_CND, hCJ1.hDTheta_CND}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->cd(2);
        draw_overlay_1D_N({hOSG.hDPhi_CND, hCJ0.hDPhi_CND, hCJ1.hDPhi_CND}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);

        c->SaveAs("Histograms/cmp_CND_residuals_1D.png");
    }

    // ------------------------------------------------------------
    // 4) Residuals: Δp, Δtheta, Δphi (particle(REC) - particle(MC))
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_residuals_MC_1D", "CND residuals MC", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hDP_REC_MC, hCJ0.hDP_REC_MC, hCJ1.hDP_REC_MC}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->cd(2);
        draw_overlay_1D_N({hOSG.hDTheta_REC_MC, hCJ0.hDTheta_REC_MC, hCJ1.hDTheta_REC_MC}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->cd(3);
        draw_overlay_1D_N({hOSG.hDPhi_REC_MC, hCJ0.hDPhi_REC_MC, hCJ1.hDPhi_REC_MC}, {"OSG", "CJ0", "CJ1"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15,0.72,0.35,0.90, /*ncols=*/1);
        c->SaveAs("Histograms/cmp_residuals_MC_1D.png");
    }

    // ------------------------------------------------------------
    // 5) 2D histograms (triptychs)
    // ------------------------------------------------------------
    draw_triptych_2D(hOSG.hPTheta, hCJ0.hPTheta, hCJ1.hPTheta,
                     "p vs theta (REC::Particle)", "Histograms/cmp_2D_p_vs_theta.png",
                     /*normalize =*/false, /*logz =*/true);

    draw_triptych_2D(hOSG.hPPhi, hCJ0.hPPhi, hCJ1.hPPhi,
                     "p vs phi (REC::Particle)", "Histograms/cmp_2D_p_vs_phi.png",
                     /*normalize =*/false, /*logz =*/true);
    draw_triptych_2D(hOSG.hThetaPhi, hCJ0.hThetaPhi, hCJ1.hThetaPhi,
                     "theta vs phi (REC::Particle)", "Histograms/cmp_2D_theta_vs_phi.png",
                     /*normalize =*/false, /*logz =*/true);

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

        fout->mkdir("OSG");
        fout->cd("OSG");
        write_hists(hOSG);

        fout->mkdir("CJ0");
        fout->cd("CJ0");
        write_hists(hCJ0);

        fout->mkdir("CJ1");
        fout->cd("CJ1");
        write_hists(hCJ1);

        fout->Write();

        std::cout << "hPhi_OSG entries = " << hOSG.hPhi->GetEntries()
                  << "  bin(0) = " << hOSG.hPhi->GetBinContent(hOSG.hPhi->FindBin(0.0))
                  << std::endl;

        int bx = hOSG.hPPhi->GetXaxis()->FindBin(1e-6); // near p=0
        int by = hOSG.hPPhi->GetYaxis()->FindBin(0.0);  // phi=0
        std::cout << "2D bin(p~0,phi=0) = " << hOSG.hPPhi->GetBinContent(bx, by) << "\n";
        std::cout << "2D integral = " << hOSG.hPPhi->Integral() << "\n";
        fout->Close();
    }
    std::cout << "Done. Wrote PNGs + cmp_cnd_versions_hists.root" << std::endl;
}