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
#include <cctype>

#include <TFile.h>
#include <cstdlib>
#include <chrono>
#include <TTree.h>
#include <TApplication.h>
#include <TDatabasePDG.h>
#include <TBenchmark.h>
#include <TF1.h>
// #include "reader.h"

#include "hipo4/reader.h"
#include "hipo4/dictionary.h"
#include "hipo4/bank.h"
#include "hipo4/event.h"

#include <sys/stat.h>
#include <sys/types.h>

#include <unordered_set>
#include <unordered_map>
#include <cstdint>

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

struct GaussCut
{
    double mean = 0.0;
    double sigma = 0.0;
    double nsigma = 3.0;
    bool valid = false;
};

GaussCut make_gauss_cut(TH1 *h, double nsigma = 3.0, double fitHalfRange = 4.0)
{
    GaussCut c;
    c.nsigma = nsigma;
    if (!h || h->GetEntries() < 200)
        return c; // protect against low stats

    // Seed around peak
    const int ib = h->GetMaximumBin();
    const double mu0 = h->GetXaxis()->GetBinCenter(ib);
    double s0 = h->GetRMS();
    if (s0 <= 0)
        s0 = 1.0;

    const double lo = mu0 - std::max(fitHalfRange, 2.5 * s0);
    const double hi = mu0 + std::max(fitHalfRange, 2.5 * s0);

    TF1 f("fgaus", "gaus", lo, hi);
    f.SetParameters(h->GetMaximum(), mu0, std::min(s0, fitHalfRange / 2.0));

    // Q: quiet, S: return fit result, R: use range, N: do not draw
    auto r = h->Fit(&f, "QSRN");
    if (int(r) != 0)
        return c;

    c.mean = f.GetParameter(1);
    c.sigma = std::fabs(f.GetParameter(2));
    c.valid = (c.sigma > 0.0);
    return c;
}

inline bool pass_gauss_cut(double x, const GaussCut &c)
{
    if (!c.valid)
        return true;
    return (std::fabs(x - c.mean) < c.nsigma * c.sigma);
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

int pick_cnd_row_angle_then_energy(
    int ip,
    const std::vector<std::vector<int>> &cndRows,
    const hipo::bank &SCINT,
    const TVector3 &pREC,
    double pvx, double pvy, double pvz,
    double maxCndAngDeg,
    bool fallbackToBestE,
    int bestScRowIp)
{
    int bestRow = -1;
    float bestE = -1.0f;
    float bestAng = 1e9f;

    if (ip < 0 || ip >= (int)cndRows.size())
        return -1;

    for (int sc : cndRows[ip])
    {
        float x = SCINT.getFloat("x", sc);
        float y = SCINT.getFloat("y", sc);
        float z = SCINT.getFloat("z", sc);

        TVector3 dir(x - pvx, y - pvy, z - pvz); // vertex-corrected hit direction
        if (dir.Mag2() < 1e-12)
            continue;

        float ang = pREC.Angle(dir) * TMath::RadToDeg();
        // if (ang > maxCndAngDeg)
        //     continue;

        float E = SCINT.getFloat("energy", sc);

        // Pick highest energy; tie-breaker: smaller angle
        if (E > bestE || (E == bestE && ang < bestAng))
        {
            bestE = E;
            bestAng = ang;
            bestRow = sc;
        }
    }

    if (bestRow >= 0)
        return bestRow;

    // If no candidate passed the angle window:
    if (fallbackToBestE)
        return bestScRowIp; // may still be -1
    return -1;              // strict mode
}

inline uint32_t parse_dst_index(const char *fname)
{
    // find the last occurrence of "dst-" and parse the number after it until the next non-digit character
    std::string s(fname);
    auto p = s.rfind("dst-");
    if (p == std::string::npos)
        return 0;
    p += 4;
    size_t q = p;
    while (q < s.size() && std::isdigit((unsigned char)s[q]))
        q++;

    return (uint32_t)std::stoul(s.substr(p, q - p));
}

using EventKey = uint64_t;

inline EventKey make_event_key_file(int dst, int event)
{
    // pack (dst,event) into 64-bit key
    return ((uint64_t)(uint32_t)dst << 32) | (uint32_t)event;
}

std::unordered_set<EventKey> collect_keys(TChain *chain, int maxEvents)
{
    std::unordered_set<EventKey> keys;
    if (!chain)
        return keys;

    TObjArray *files = chain->GetListOfFiles();
    if (!files || files->GetLast() < 0)
        return keys;

    hipo::reader reader;
    hipo::dictionary factory;
    hipo::event event;

    int processed = 0;
    for (int fi = 0; fi <= files->GetLast(); ++fi)
    {
        const char *fname = files->At(fi)->GetTitle();
        uint32_t dst = parse_dst_index(fname);

        reader.open(fname);
        reader.readDictionary(factory);
        if (!factory.hasSchema("RUN::config"))
            continue;
        hipo::bank CONF(factory.getSchema("RUN::config"));

        while (reader.next())
        {
            reader.read(event);
            event.getStructure(CONF);
            uint32_t evn = (uint32_t)CONF.getInt("event", 0);

            keys.insert(make_event_key_file(dst, evn));

            processed++;
            if (maxEvents > 0 && processed >= maxEvents)
                break;
        }
        if (maxEvents > 0 && processed >= maxEvents)
            break;
    }
    return keys;
}

struct Hists
{
    TH1F *hP;
    TH1F *hTheta;
    TH1F *hPhi;

    TH1F *hEnergy_CND;
    TH1F *hTheta_CND;
    TH1F *hPhi_CND;

    TH1F *hCNDLayer;           // Which CND layer contains the highest-energy CND hit for this neutron?
    TH1F *hCNDLayer_occupancy; // Among selected neutrons, how often is layer 1/2/3 present? Neutrons hitting 2 layers contribute twice.
    TH1F *hCND_NLayers;        // number of layers hit per neutron (1, 2, or 3)

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

inline float phi_0_360_deg(float phi_rad)
{
    // phi_rad typically in [-pi, pi], convert to [0, 2pi) then to degrees
    return TVector2::Phi_0_2pi(phi_rad) * TMath::RadToDeg();
}

Hists book_hists(const char *tag)
{
    Hists h;
    h.hP = new TH1F(TString::Format("hP_%s", tag), TString::Format("Neutron p (%s) from REC::Particle;p [GeV];Counts", tag), 100, 0, 5);
    h.hTheta = new TH1F(TString::Format("hTh_%s", tag), TString::Format("Neutron #theta (%s) from REC::Particle;#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhi = new TH1F(TString::Format("hPhi_%s", tag), TString::Format("Neutron #phi (%s) from REC::Particle;#phi [deg];Counts", tag), 100, 0, 361);
    // h.hP = new TH1F(TString::Format("hP_%s", tag), TString::Format("Neutron p (%s) from REC::Particle;p [GeV];Counts", tag), 100, 0, 1.3);
    // h.hTheta = new TH1F(TString::Format("hTh_%s", tag), TString::Format("Neutron #theta (%s) from REC::Particle;#theta [deg];Counts", tag), 100, 35, 115);
    // h.hPhi = new TH1F(TString::Format("hPhi_%s", tag), TString::Format("Neutron #phi (%s) from REC::Particle;#phi [deg];Counts", tag), 100, 90, 150);

    h.hEnergy_CND = new TH1F(TString::Format("hECND_%s", tag), TString::Format("Neutron E (%s) from CND-cluster;E [GeV];Counts", tag), 100, 0, 12);
    h.hTheta_CND = new TH1F(TString::Format("hThCND_%s", tag), TString::Format("Neutron #theta (%s) from CND-cluster;#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhi_CND = new TH1F(TString::Format("hPhiCND_%s", tag), TString::Format("Neutron #phi (%s) from CND-cluster;#phi [deg];Counts", tag), 100, 0, 361);
    // h.hEnergy_CND = new TH1F(TString::Format("hECND_%s", tag), TString::Format("Neutron E (%s) from CND-cluster;E [GeV];Counts", tag), 100, 0, 12);
    // h.hTheta_CND = new TH1F(TString::Format("hThCND_%s", tag), TString::Format("Neutron #theta (%s) from CND-cluster;#theta [deg];Counts", tag), 100, 0, 180);
    // h.hPhi_CND = new TH1F(TString::Format("hPhiCND_%s", tag), TString::Format("Neutron #phi (%s) from CND-cluster;#phi [deg];Counts", tag), 100, 90, 150);

    h.hCNDLayer = new TH1F(TString::Format("hLayerCND_%s", tag), TString::Format("CND layer (%s) for neutron clusters;Layer;Counts", tag), 3, 0.5, 3.5);
    h.hCNDLayer_occupancy = new TH1F(TString::Format("hLayerCND_mult_%s", tag), TString::Format("CND layer occupancy (%s) for neutron clusters;Layer;Counts", tag), 3, 0.5, 3.5);
    h.hCND_NLayers = new TH1F(TString::Format("hNLayerCND_%s", tag), TString::Format("Number of CND layers hit (%s) for neutron clusters;N layers;Counts", tag), 3, 0.5, 3.5);

    h.hDTheta_CND = new TH1F(TString::Format("hDThCND_%s", tag), "#Delta#theta(CND-cluster - particle);#Delta#theta [deg];Counts", 50, -5, 5);
    h.hDPhi_CND = new TH1F(TString::Format("hDPhCND_%s", tag), "#Delta#phi(CND-cluster - particle);#Delta#phi [deg];Counts", 50, -5, 5);

    h.hPTheta = new TH2F(TString::Format("hPTh_%s", tag), TString::Format("#theta vs p (%s) from REC::Particle;p [GeV];#theta [deg]", tag),
                         100, 0, 10, 100, 0, 180);
    h.hPPhi = new TH2F(TString::Format("hPPhi_%s", tag), TString::Format("#phi vs p (%s) from REC::Particle;p [GeV];#phi [deg]", tag),
                       100, 0, 10, 100, 0, 361);
    h.hThetaPhi = new TH2F(TString::Format("hThPhi_%s", tag), TString::Format("#phi vs #theta (%s) from REC::Particle;#theta [deg];#phi [deg]", tag),
                           100, 0, 180, 100, 0, 361);

    h.hPMC = new TH1F(TString::Format("hPMC_%s", tag), TString::Format("MC Neutron p (%s);p [GeV];Counts", tag), 100, 0, 5);
    h.hThetaMC = new TH1F(TString::Format("hThMC_%s", tag), TString::Format("MC Neutron #theta (%s);#theta [deg];Counts", tag), 100, 0, 180);
    h.hPhiMC = new TH1F(TString::Format("hPhiMC_%s", tag), TString::Format("MC Neutron #phi (%s);#phi [deg];Counts", tag), 100, 0, 361);
    // h.hPMC = new TH1F(TString::Format("hPMC_%s", tag), TString::Format("MC Neutron p (%s);p [GeV];Counts", tag), 100, 0, 1.3);
    // h.hThetaMC = new TH1F(TString::Format("hThMC_%s", tag), TString::Format("MC Neutron #theta (%s);#theta [deg];Counts", tag), 100, 35, 115);
    // h.hPhiMC = new TH1F(TString::Format("hPhiMC_%s", tag), TString::Format("MC Neutron #phi (%s);#phi [deg];Counts", tag), 100, 90, 150);

    h.hDTheta_REC_MC = new TH1F(TString::Format("hDTh_REC_MC_%s", tag), "#Delta#theta(REC - MC);#Delta#theta [deg];Counts", 200, -180, 180);
    h.hDPhi_REC_MC = new TH1F(TString::Format("hDPh_REC_MC_%s", tag), "#Delta#phi(REC - MC);#Delta#phi [deg];Counts", 200, -180, 180);
    h.hAngle_REC_MC = new TH1F(TString::Format("hAngle_REC_MC_%s", tag), "Opening angle(REC - MC);Opening angle [deg];Counts", 200, 0, 180);
    h.hDP_REC_MC = new TH1F(TString::Format("hDP_REC_MC_%s", tag), "#Delta p(REC) - p(MC);#Delta p [GeV];Counts", 200, -1.5, 1.0);

    return h;
}

void write_hists(const Hists &h)
{
    std::vector<TObject *> objs = {
        h.hP, h.hTheta, h.hPhi,
        h.hEnergy_CND, h.hTheta_CND, h.hPhi_CND,
        h.hCNDLayer, h.hCNDLayer_occupancy, h.hCND_NLayers,
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

void process_chain_allow_keys(
    TChain *chain, Hists &h, const char *tag,
    int maxEvents, int cnd_id, double maxAngleDeg,
    const GaussCut *dThetaCut,
    const std::unordered_set<EventKey> &allowKeys)
{
    if (!chain)
        return;

    TObjArray *files = chain->GetListOfFiles();
    hipo::reader reader;
    hipo::dictionary factory;
    hipo::event event;
    TLorentzVector N_Vec_temp;
    const double Nmass = TDatabasePDG::Instance()->GetParticle(2112)->Mass();

    long long scanned = 0;
    long long kept = 0;
    for (int fi = 0; fi <= files->GetLast(); ++fi)
    {
        TObjArray *files = chain->GetListOfFiles();
        if (!files || files->GetLast() < 0)
            return;

        const char *fname = files->At(fi)->GetTitle();
        uint32_t dst = parse_dst_index(fname);

        reader.open(fname);
        reader.readDictionary(factory);

        // if (factory.hasSchema("MC::RecMatch"))
        // {
        //     std::cout << "\n--- MC::RecMatch schema ---\n";
        //     factory.getSchema("MC::RecMatch").show();
        // }

        if (!factory.hasSchema("RUN::config") ||
            !factory.hasSchema("REC::Particle") ||
            !factory.hasSchema("REC::Scintillator") ||
            !factory.hasSchema("MC::Particle"))
            continue;

        bool hasRecMatch = factory.hasSchema("MC::RecMatch");
        hipo::bank RECM(hasRecMatch ? factory.getSchema("MC::RecMatch")
                                    : factory.getSchema("RUN::config")); // dummy schema if missing

        hipo::bank CONF(factory.getSchema("RUN::config"));
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank SCINT(factory.getSchema("REC::Scintillator"));
        hipo::bank MCPT(factory.getSchema("MC::Particle"));

        while (reader.next() && (maxEvents <= 0 || scanned < maxEvents))
        {
            reader.read(event);
            event.getStructure(CONF);

            scanned++;
            if (maxEvents > 0 && scanned >= maxEvents)
                break;

            uint32_t evn = (uint32_t)CONF.getInt("event", 0);
            EventKey key = make_event_key_file(dst, evn);

            if (!allowKeys.count(key))
            {
                continue;
            }

            kept++;

            if (hasRecMatch)
                event.getStructure(RECM);
            event.getStructure(PART);
            event.getStructure(SCINT);
            event.getStructure(MCPT);

            const int nPart = PART.getRows();
            const int nSc = SCINT.getRows();
            const int nMC = MCPT.getRows();

            std::vector<int> bestScRow(nPart, -1);
            std::vector<float> bestE(nPart, -1.0f);
            std::vector<int> cndLayerMask(nPart, 0); // bitmask of layers with associated scint hits

            // Fill per-particle CND info from scintillator rows (detector == CND)
            std::vector<std::vector<int>> cndRows(nPart);

            for (int sc = 0; sc < nSc; ++sc)
            {
                const int detector = SCINT.getInt("detector", sc);
                if (detector != cnd_id)
                    continue;

                const int pindex = SCINT.getInt("pindex", sc);
                if (pindex < 0 || pindex >= nPart)
                    continue;

                cndRows[pindex].push_back(sc);

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

                // Found a CND-related scintillator cluster
                const int layer = SCINT.getInt("layer", sc);
                if (layer >= 1 && layer <= 3)
                {
                    cndLayerMask[pindex] |= (1 << (layer - 1));
                }
            }

            std::unordered_map<int, int> p2mc;
            std::unordered_map<int, float> p2q;

            if (hasRecMatch)
            {
                const int nRM = RECM.getRows();

                for (int ir = 0; ir < nRM; ++ir)
                {
                    const int pidx = RECM.getInt("pindex", ir);
                    const int mcidx = RECM.getInt("mcindex", ir);
                    const float quality = RECM.getFloat("quality", ir);

                    if (mcidx < 0)
                        continue; // ignore unmatched rows
                    // if (quality < 0.5)
                    //     continue;

                    // Keep best score if multiple rows map to same pindex
                    auto itq = p2q.find(pidx);
                    if (itq == p2q.end() || quality > itq->second)
                    {
                        p2q[pidx] = quality;
                        p2mc[pidx] = mcidx;
                    }
                }
            }

            // static bool once = true;
            // if (once && hasRecMatch)
            // {
            //     std::cout << "RecMatch rows: " << RECM.getRows() << "\n";
            //     for (int k = 0; k < 5; k++)
            //     {
            //         std::cout << "  pindex=" << RECM.getInt("pindex", k)
            //                   << " mcindex=" << RECM.getInt("mcindex", k)
            //                   << " quality=" << RECM.getFloat("quality", k) << "\n";
            //     }
            //     once = false;
            // }

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

                if (pREC.Mag() < 1e-6 && pREC.Theta() * TMath::RadToDeg() < 1e-6 && phi_0_360_deg(pREC.Phi()) < 1e-6)
                    continue; // skip zero-momentum

                MatchResult match;
                auto it = p2mc.find(ip);

                if (it != p2mc.end())
                {
                    const int mcidx = it->second;

                    // ensure mcidx is valid and points to a neutron
                    if (mcidx >= 0 && mcidx < MCPT.getRows() && MCPT.getInt("pid", mcidx) == 2112)
                    {
                        match.mcIndex = mcidx;
                        TVector3 pMC(MCPT.getFloat("px", mcidx),
                                     MCPT.getFloat("py", mcidx),
                                     MCPT.getFloat("pz", mcidx));
                        match.angleDeg = pREC.Angle(pMC) * TMath::RadToDeg();

                        // if (maxAngleDeg > 0 && match.angleDeg > maxAngleDeg)
                        //     match.mcIndex = -1;
                    }
                }

                // if (match.mcIndex < 0)
                // {
                //     match = match_neutron_rec_to_mc(pREC, MCPT, /*maxAngleDeg=*/maxAngleDeg);
                // }
                // if (match.mcIndex < 0)
                //     continue; // no good MC match

                // MatchResult match = match_neutron_rec_to_mc(pREC, MCPT, /*maxAngleDeg=*/maxAngleDeg);
                // if (match.mcIndex < 0)
                //     continue; // no good MC match

                const float p_rec = pREC.Mag();
                const float th_rec = pREC.Theta() * TMath::RadToDeg();
                const float ph_rec = phi_0_360_deg(pREC.Phi());

                // if (std::abs(ph_rec) < 1e-6)
                // {
                //     std::cout << "phi=0: px=" << px << " py=" << py << " pz=" << pz << "\n";
                // }

                N_Vec_temp.SetPxPyPzE(px, py, pz, TMath::Sqrt(p_rec * p_rec + Nmass * Nmass));

                const float mc_px = MCPT.getFloat("px", match.mcIndex);
                const float mc_py = MCPT.getFloat("py", match.mcIndex);
                const float mc_pz = MCPT.getFloat("pz", match.mcIndex);

                TVector3 pMC(mc_px, mc_py, mc_pz);

                const float th_mc = pMC.Theta() * TMath::RadToDeg();
                const float ph_mc = phi_0_360_deg(pMC.Phi());

                const float dth_rec_mc = th_rec - th_mc;
                const float dph_rec_mc = TVector2::Phi_mpi_pi((ph_rec - ph_mc) * TMath::DegToRad()) * TMath::RadToDeg();

                // Now fill CND-related info if we have an associated row
                const bool fallbackToBestE = false; // recommended: keep strict to kill tails

                int scBest = pick_cnd_row_angle_then_energy(
                    ip, cndRows, SCINT,
                    pREC, pvx, pvy, pvz,
                    /*maxCndAngDeg=*/maxAngleDeg,
                    fallbackToBestE,
                    bestScRow[ip]);

                // if (scBest < 0)
                //     continue; // no good CND cluster for this neutron

                // if (dThetaCut && !pass_gauss_cut(dth_rec_mc, *dThetaCut))
                //     continue;

                const int layerBest = SCINT.getInt("layer", scBest);
                h.hCNDLayer->Fill(layerBest);

                int nLayers = 0;
                for (int layer = 1; layer <= 3; ++layer)
                {
                    if (cndLayerMask[ip] & (1 << (layer - 1)))
                    {
                        h.hCNDLayer_occupancy->Fill(layer);
                        nLayers++;
                    }
                }

                if (nLayers > 0)
                    h.hCND_NLayers->Fill(nLayers);

                const float E = SCINT.getFloat("energy", scBest);
                const float x = SCINT.getFloat("x", scBest);
                const float y = SCINT.getFloat("y", scBest);
                const float z = SCINT.getFloat("z", scBest);

                TVector3 r(x, y, z);
                TVector3 dir = r - TVector3(pvx, pvy, pvz); // vertex-corrected
                const float th_cluster = dir.Theta() * TMath::RadToDeg();
                const float ph_cluster = phi_0_360_deg(dir.Phi());

                // Delta phi with wrapping
                const float dth = th_cluster - th_rec;
                const float dph = TVector2::Phi_mpi_pi((ph_cluster - ph_rec) * TMath::DegToRad()) * TMath::RadToDeg();

                // if (pREC.Mag() < 0.30)
                // {
                //     float th_mc_dbg = th_mc;
                //     float ph_mc_dbg = ph_mc;
                //     std::cout
                //         << "LOW-P survivor: p_rec=" << p_rec
                //         << " th_rec=" << th_rec << " ph_rec=" << ph_rec
                //         << " |  th_mc=" << th_mc_dbg << " ph_mc=" << ph_mc_dbg
                //         << " | angle=" << match.angleDeg
                //         << " | quality=" << p2q[ip]
                //         << "\n";
                // }

                h.hP->Fill(p_rec);
                h.hTheta->Fill(th_rec);
                h.hPhi->Fill(ph_rec);
                h.hPTheta->Fill(p_rec, th_rec);
                h.hPPhi->Fill(p_rec, ph_rec);
                h.hThetaPhi->Fill(th_rec, ph_rec);

                h.hDTheta_REC_MC->Fill(dth_rec_mc);
                h.hDPhi_REC_MC->Fill(dph_rec_mc);
                h.hAngle_REC_MC->Fill(match.angleDeg);
                h.hDP_REC_MC->Fill(p_rec - pMC.Mag());

                h.hPMC->Fill(pMC.Mag());
                h.hThetaMC->Fill(th_mc);
                h.hPhiMC->Fill(ph_mc);

                h.hEnergy_CND->Fill(E);
                h.hTheta_CND->Fill(th_cluster);
                h.hPhi_CND->Fill(ph_cluster);
                h.hDTheta_CND->Fill(dth);
                h.hDPhi_CND->Fill(dph);
            }

            if (maxEvents > 0 && kept > maxEvents)
                break;
        }
        if (maxEvents > 0 && kept > maxEvents)
            break;
    }

    std::cout << "Processed " << kept << " events for " << tag << "\n";
    std::cout << tag << ": scanned=" << scanned << "\n";
}

void process_chain(TChain *chain,
                   Hists &h,
                   const char *tag,
                   int maxEvents = 300000,
                   int cnd_id = 3,
                   double maxAngleDeg = 10.0,
                   const GaussCut *dThetaCut = nullptr,
                   const std::unordered_set<EventKey> *eventVeto = nullptr)
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

        // static bool printedRM = false;
        // if (factory.hasSchema("MC::RecMatch") && !printedRM)
        // {
        //     std::cout << "\n--- MC::RecMatch schema ---\n";
        //     factory.getSchema("MC::RecMatch").show();
        //     printedRM = true;
        // }

        uint32_t dst = parse_dst_index(fname);

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

        bool hasRecMatch = factory.hasSchema("MC::RecMatch");
        hipo::bank RECM(hasRecMatch ? factory.getSchema("MC::RecMatch")
                                    : factory.getSchema("RUN::config")); // dummy schema if missing

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
            int RunNumber = CONF.getInt("run", 0);
            int EventNumber = CONF.getInt("event", 0);

            EventKey key = make_event_key_file(dst, EventNumber);
            if (eventVeto && eventVeto->count(key))
            {
                // This OSG event exists in the other dataset, so skip it (we only want missing)
                continue;
            }

            event.getStructure(PART);
            event.getStructure(SCINT);
            event.getStructure(MCPT);
            if (hasRecMatch)
                event.getStructure(RECM);
            // if (hasScintX) event.getStructure(SCINTX);

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
            std::vector<int> cndLayerMask(nPart, 0); // bitmask of layers with associated scint hits

            // std::vector<double> part_Scint_CND_E(nPart, NAN);
            // std::vector<double> part_Scint_CND_t(nPart, NAN);
            // std::vector<double> part_Scint_CND_x(nPart, NAN);
            // std::vector<double> part_Scint_CND_y(nPart, NAN);
            // std::vector<double> part_Scint_CND_z(nPart, NAN);
            // std::vector<double> part_ScintX_CND_dedx(nPart, NAN);
            // std::vector<double> part_ScintX_CND_size(nPart, NAN);
            // std::vector<double> part_ScintX_CND_layermult(nPart, NAN);

            // Fill per-particle CND info from scintillator rows (detector == CND)
            std::vector<std::vector<int>> cndRows(nPart);

            for (int sc = 0; sc < nSc; ++sc)
            {
                const int detector = SCINT.getInt("detector", sc);
                if (detector != cnd_id)
                    continue;

                const int pindex = SCINT.getInt("pindex", sc);
                if (pindex < 0 || pindex >= nPart)
                    continue;

                cndRows[pindex].push_back(sc);

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

                // Found a CND-related scintillator cluster
                const int layer = SCINT.getInt("layer", sc);
                if (layer >= 1 && layer <= 3)
                {
                    cndLayerMask[pindex] |= (1 << (layer - 1));
                }
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

            std::unordered_map<int, int> p2mc;
            std::unordered_map<int, float> p2q;

            if (hasRecMatch)
            {
                const int nRM = RECM.getRows();

                for (int ir = 0; ir < nRM; ++ir)
                {
                    const int pidx = RECM.getInt("pindex", ir);
                    const int mcidx = RECM.getInt("mcindex", ir);
                    const float quality = RECM.getFloat("quality", ir);

                    if (mcidx < 0)
                        continue; // ignore unmatched rows

                    // if (quality < 0.5)
                    //     continue;

                    // Keep best score if multiple rows map to same pindex
                    auto itq = p2q.find(pidx);
                    if (itq == p2q.end() || quality > itq->second)
                    {
                        p2q[pidx] = quality;
                        p2mc[pidx] = mcidx;
                    }
                }
            }

            // static bool once = true;
            // if (once && hasRecMatch)
            // {
            //     std::cout << "RecMatch rows: " << RECM.getRows() << "\n";
            //     for (int k = 0; k < 5; k++)
            //     {
            //         std::cout << "  pindex=" << RECM.getInt("pindex", k)
            //                   << " mcindex=" << RECM.getInt("mcindex", k)
            //                   << " quality=" << RECM.getFloat("quality", k) << "\n";
            //     }
            //     once = false;
            // }

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

                if (pREC.Mag() < 1e-6 && pREC.Theta() * TMath::RadToDeg() < 1e-6 && phi_0_360_deg(pREC.Phi()) < 1e-6)
                    continue; // skip zero-momentum

                // First try truth link
                MatchResult match;
                auto it = p2mc.find(ip);

                if (it != p2mc.end())
                {
                    const int mcidx = it->second;

                    // ensure mcidx is valid and points to a neutron
                    if (mcidx >= 0 && mcidx < MCPT.getRows() && MCPT.getInt("pid", mcidx) == 2112)
                    {
                        match.mcIndex = mcidx;
                        TVector3 pMC(MCPT.getFloat("px", mcidx),
                                     MCPT.getFloat("py", mcidx),
                                     MCPT.getFloat("pz", mcidx));
                        match.angleDeg = pREC.Angle(pMC) * TMath::RadToDeg();

                        // if (maxAngleDeg > 0 && match.angleDeg > maxAngleDeg)
                        //     match.mcIndex = -1;
                    }
                }

                // if (match.mcIndex < 0)
                // {
                //     match = match_neutron_rec_to_mc(pREC, MCPT, /*maxAngleDeg=*/maxAngleDeg);
                // }
                // if (match.mcIndex < 0)
                //     continue; // no good MC match

                const float p_rec = pREC.Mag();
                const float th_rec = pREC.Theta() * TMath::RadToDeg();
                const float ph_rec = phi_0_360_deg(pREC.Phi());

                // if (std::abs(ph_rec) < 1e-6)
                // {
                //     std::cout << "phi=0: px=" << px << " py=" << py << " pz=" << pz << "\n";
                // }

                N_Vec_temp.SetPxPyPzE(px, py, pz, TMath::Sqrt(p_rec * p_rec + Nmass * Nmass));

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
                const float ph_mc = phi_0_360_deg(pMC.Phi());

                const float dth_rec_mc = th_rec - th_mc;
                const float dph_rec_mc = TVector2::Phi_mpi_pi((ph_rec - ph_mc) * TMath::DegToRad()) * TMath::RadToDeg();

                // Now fill CND-related info if we have an associated row
                // Now fill CND-related info if we have an associated row
                const bool fallbackToBestE = false; // recommended: keep strict to kill tails
                // const double maxCndAngDeg = 15.0;

                int scBest = pick_cnd_row_angle_then_energy(
                    ip, cndRows, SCINT,
                    pREC, pvx, pvy, pvz,
                    /*maxCndAngDeg=*/maxAngleDeg,
                    fallbackToBestE,
                    bestScRow[ip]);

                // if (scBest < 0)
                //     continue; // no good CND cluster for this neutron

                // if (dThetaCut && !pass_gauss_cut(dth_rec_mc, *dThetaCut))
                //     continue;

                const int layerBest = SCINT.getInt("layer", scBest);
                h.hCNDLayer->Fill(layerBest);

                int nLayers = 0;
                for (int layer = 1; layer <= 3; ++layer)
                {
                    if (cndLayerMask[ip] & (1 << (layer - 1)))
                    {
                        h.hCNDLayer_occupancy->Fill(layer);
                        nLayers++;
                    }
                }

                if (nLayers > 0)
                    h.hCND_NLayers->Fill(nLayers);

                const float E = SCINT.getFloat("energy", scBest);
                const float x = SCINT.getFloat("x", scBest);
                const float y = SCINT.getFloat("y", scBest);
                const float z = SCINT.getFloat("z", scBest);

                TVector3 r(x, y, z);
                TVector3 dir = r - TVector3(pvx, pvy, pvz); // vertex-corrected
                const float th_cluster = dir.Theta() * TMath::RadToDeg();
                const float ph_cluster = phi_0_360_deg(dir.Phi());

                // Delta phi with wrapping
                const float dth = th_cluster - th_rec;
                const float dph = TVector2::Phi_mpi_pi((ph_cluster - ph_rec) * TMath::DegToRad()) * TMath::RadToDeg();

                if ((tag == nullptr || std::string(tag) != "OSG_cal") && (pREC.Mag() < 1e-6 || pREC.Theta() * TMath::RadToDeg() < 1e-6 || phi_0_360_deg(pREC.Phi()) < 1e-6))
                {
                    float th_mc_dbg = th_mc;
                    float ph_mc_dbg = ph_mc;
                    std::cout
                        << "LOW-P survivor: p_rec=" << p_rec
                        << " th_rec=" << th_rec << " ph_rec=" << ph_rec
                        << " |  th_mc=" << th_mc_dbg << " ph_mc=" << ph_mc_dbg
                        << " | angle=" << match.angleDeg
                        << " | quality=" << p2q[ip]
                        << "\n";
                }

                h.hP->Fill(p_rec);
                h.hTheta->Fill(th_rec);
                h.hPhi->Fill(ph_rec);
                h.hPTheta->Fill(p_rec, th_rec);
                h.hPPhi->Fill(p_rec, ph_rec);
                h.hThetaPhi->Fill(th_rec, ph_rec);

                h.hDTheta_REC_MC->Fill(dth_rec_mc);
                h.hDPhi_REC_MC->Fill(dph_rec_mc);
                h.hAngle_REC_MC->Fill(match.angleDeg);
                h.hDP_REC_MC->Fill(p_rec - pMC.Mag());

                h.hPMC->Fill(pMC.Mag());
                h.hThetaMC->Fill(th_mc);
                h.hPhiMC->Fill(ph_mc);

                h.hEnergy_CND->Fill(E);
                h.hTheta_CND->Fill(th_cluster);
                h.hPhi_CND->Fill(ph_cluster);
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

void addEntriesColumnMajor(TLegend *leg,
                           const std::vector<TH1 *> &h,
                           const std::vector<const char *> &labels,
                           int ncols,
                           const char *opt = "f")
{
    const int n = (int)h.size();
    const int nrows = (int)std::ceil((double)n / ncols);

    // Column-major order: (row, col) -> i = col*nrows + row
    // Add in the order: row=0..nrows-1, col=0..ncols-1
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            int i = col * nrows + row;
            if (i >= n)
                continue;
            leg->AddEntry(h[i], labels[i], opt);
        }
    }
}

void draw_overlay_1D_N(const std::vector<TH1 *> &h_in,
                       const std::vector<const char *> &labels,
                       bool normalize = true,
                       bool logy = true,
                       double leg_x1 = 0.15, double leg_y1 = 0.72,
                       double leg_x2 = 0.35, double leg_y2 = 0.90,
                       int leg_ncols = 1,
                       const char *plot_title = nullptr)
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

    // --- Override title ONLY for the plot (PNG), originals remain unchanged
    if (plot_title && plot_title[0] != '\0')
        h[0]->SetTitle(plot_title);

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
        kRed,
        kBlue,
        kGreen + 2,
        kMagenta - 8,
        kOrange - 3,
        kMagenta + 1,
        kViolet + 1,
        kGray + 2};

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
        h[0]->SetMaximum(1.2 * m);

    // Draw
    h[0]->Draw("hist");
    for (size_t i = 1; i < h.size(); ++i)
        h[i]->Draw("hist same");

    // Legend
    auto *leg = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
    leg->SetNColumns(leg_ncols);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    // for (size_t i = 0; i < h.size(); ++i)
    //     leg->AddEntry(h[i], labels[i], "f");

    addEntriesColumnMajor(leg, h, labels, leg_ncols, "f");
    leg->Draw();
}

// Helper: draw 3 TH2 side-by-side
void draw_triptych_2D(TH2 *h1_in, TH2 *h2_in, TH2 *h3_in, TH2 *h4_in,
                      const char *title,
                      const char *outname,
                      bool normalize = true,
                      bool logz = true)
{
    if (!h1_in || !h2_in || !h3_in || !h4_in)
        return;

    // Clone to avoid modifying originals
    auto *h1 = (TH2 *)h1_in->Clone(Form("%s_clone1", h1_in->GetName()));
    auto *h2 = (TH2 *)h2_in->Clone(Form("%s_clone2", h2_in->GetName()));
    auto *h3 = (TH2 *)h3_in->Clone(Form("%s_clone3", h3_in->GetName()));
    auto *h4 = (TH2 *)h4_in->Clone(Form("%s_clone4", h4_in->GetName()));

    h1->SetDirectory(nullptr);
    h2->SetDirectory(nullptr);
    h3->SetDirectory(nullptr);
    h4->SetDirectory(nullptr);

    // Optional normalization (global integral)
    if (normalize)
    {
        if (h1->Integral() > 0)
            h1->Scale(1.0 / h1->Integral());
        if (h2->Integral() > 0)
            h2->Scale(1.0 / h2->Integral());
        if (h3->Integral() > 0)
            h3->Scale(1.0 / h3->Integral());
        if (h4->Integral() > 0)
            h4->Scale(1.0 / h4->Integral());
    }

    // gPad->Modified();
    // gPad->Update();

    auto *c = new TCanvas(Form("c2D_%s", outname), title, 1800, 900);
    c->Divide(2, 2);

    c->cd(1);
    gPad->SetRightMargin(0.14);
    // set z max range user 10^3
    if (logz)
        gPad->SetLogz();
    h1->SetMaximum(1e3);
    h1->SetStats(0);
    h1->Draw("colz");

    c->cd(2);
    gPad->SetRightMargin(0.14);
    if (logz)
        gPad->SetLogz();
    h2->SetMaximum(1e3);
    h2->SetStats(0);
    h2->Draw("colz");

    c->cd(3);
    gPad->SetRightMargin(0.14);
    if (logz)
        gPad->SetLogz();
    h3->SetMaximum(1e3);
    h3->SetStats(0);
    h3->Draw("colz");

    c->cd(4);
    gPad->SetRightMargin(0.14);
    if (logz)
        gPad->SetLogz();
    h4->SetMaximum(1e3);
    h4->SetStats(0);
    h4->Draw("colz");

    c->SaveAs(outname);
}

void print_counts(const char *tag, const Hists &h)
{
    auto ll = [](double x)
    { return (long long)std::llround(x); };

    std::cout << "\n==================== " << tag << " ====================\n";

    std::cout << "REC::Particle neutrons (entries):\n";
    std::cout << "  hP     : " << ll(h.hP->GetEntries()) << "\n";
    std::cout << "  hTheta : " << ll(h.hTheta->GetEntries()) << "\n";
    std::cout << "  hPhi   : " << ll(h.hPhi->GetEntries()) << "\n";

    std::cout << "MC::Particle matched neutrons (entries):\n";
    std::cout << "  hPMC     : " << ll(h.hPMC->GetEntries()) << "\n";
    std::cout << "  hThetaMC : " << ll(h.hThetaMC->GetEntries()) << "\n";
    std::cout << "  hPhiMC   : " << ll(h.hPhiMC->GetEntries()) << "\n";

    std::cout << "Neutrons with a matched CND cluster (entries):\n";
    std::cout << "  hEnergy_CND : " << ll(h.hEnergy_CND->GetEntries()) << "\n";
    std::cout << "  hTheta_CND  : " << ll(h.hTheta_CND->GetEntries()) << "\n";
    std::cout << "  hPhi_CND    : " << ll(h.hPhi_CND->GetEntries()) << "\n";

    std::cout << "Residuals (entries):\n";
    std::cout << "  hDTheta_CND     : " << ll(h.hDTheta_CND->GetEntries()) << "\n";
    std::cout << "  hDPhi_CND       : " << ll(h.hDPhi_CND->GetEntries()) << "\n";
    std::cout << "  hDTheta_REC_MC  : " << ll(h.hDTheta_REC_MC->GetEntries()) << "\n";
    std::cout << "  hDPhi_REC_MC    : " << ll(h.hDPhi_REC_MC->GetEntries()) << "\n";
    std::cout << "  hAngle_REC_MC   : " << ll(h.hAngle_REC_MC->GetEntries()) << "\n";
    std::cout << "  hDP_REC_MC      : " << ll(h.hDP_REC_MC->GetEntries()) << "\n";
}

TH1F *make_diff_hist(const TH1F *hOSG, const TH1F *hCJ, const char *name)
{
    if (!hOSG || !hCJ)
        return nullptr;
    auto *h = (TH1F *)hOSG->Clone(name);
    h->SetDirectory(nullptr);
    h->Add((TH1 *)hCJ, -1.0);
    return h;
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
    TChain *chCJ2 = make_chain(BASE, "testCJ2");

    // Book hists
    Hists hOSG = book_hists("OSG");
    Hists hCJ0 = book_hists("CJ0");
    Hists hCJ1 = book_hists("CJ1");
    Hists hCJ2 = book_hists("CJ2");

    std::cout << "OSG chain files = "
              << (chOSG->GetListOfFiles() ? chOSG->GetListOfFiles()->GetEntries() : -1)
              << std::endl;

    // fill dtheta with a loose matching to get a rough cut, then apply a Gaussian cut to clean it up and see the effect on the other variables
    Hists hOSG_cal = book_hists("OSG_cal");
    process_chain(chOSG, hOSG_cal, "OSG_cal", maxEvents, 3, /*maxAngleDeg=*/180.0, /*dThetaCut=*/nullptr);

    GaussCut thCut = make_gauss_cut(hOSG_cal.hDTheta_REC_MC, /*nsigma=*/3.0, /*fitHalfRange=*/4.0);
    std::cout << "Derived dTheta cut from OSG: mean=" << thCut.mean
              << " sigma=" << thCut.sigma
              << " => |dTheta-mean| < " << (thCut.nsigma * thCut.sigma) << " deg\n";
    if (!thCut.valid)
    {
        std::cout << "WARNING: thCut invalid (fit failed or low stats). No dTheta cut will be applied.\n";
    }

    // GaussCut thCut1;
    // thCut1.mean = -6.28171;
    // thCut1.sigma = 20.1061; // default to 2 deg if fit failed or low stats
    // thCut1.nsigma = 3.0;
    // thCut1.valid = true;

    process_chain(chOSG, hOSG, "OSG", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, &thCut);
    process_chain(chCJ0, hCJ0, "CJ0", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, &thCut);
    process_chain(chCJ1, hCJ1, "CJ1", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, &thCut);
    process_chain(chCJ2, hCJ2, "CJ2", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, &thCut);

    // process_chain(chOSG, hOSG, "OSG", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, /*dThetaCut=*/nullptr);
    // process_chain(chCJ0, hCJ0, "CJ0", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, /*dThetaCut=*/nullptr);
    // process_chain(chCJ1, hCJ1, "CJ1", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, /*dThetaCut=*/nullptr);
    // process_chain(chCJ2, hCJ2, "CJ2", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, /*dThetaCut=*/nullptr);

    // --- Bin-by-bin differences (OSG - CJ)
    auto *dE_CJ0 = make_diff_hist(hOSG.hEnergy_CND, hCJ0.hEnergy_CND, "dE_OSG_minus_CJ0");
    auto *dE_CJ1 = make_diff_hist(hOSG.hEnergy_CND, hCJ1.hEnergy_CND, "dE_OSG_minus_CJ1");
    auto *dE_CJ2 = make_diff_hist(hOSG.hEnergy_CND, hCJ2.hEnergy_CND, "dE_OSG_minus_CJ2");

    auto *dTh_CJ0 = make_diff_hist(hOSG.hTheta_CND, hCJ0.hTheta_CND, "dTh_OSG_minus_CJ0");
    auto *dTh_CJ1 = make_diff_hist(hOSG.hTheta_CND, hCJ1.hTheta_CND, "dTh_OSG_minus_CJ1");
    auto *dTh_CJ2 = make_diff_hist(hOSG.hTheta_CND, hCJ2.hTheta_CND, "dTh_OSG_minus_CJ2");

    auto *dPh_CJ0 = make_diff_hist(hOSG.hPhi_CND, hCJ0.hPhi_CND, "dPh_OSG_minus_CJ0");
    auto *dPh_CJ1 = make_diff_hist(hOSG.hPhi_CND, hCJ1.hPhi_CND, "dPh_OSG_minus_CJ1");
    auto *dPh_CJ2 = make_diff_hist(hOSG.hPhi_CND, hCJ2.hPhi_CND, "dPh_OSG_minus_CJ2");

    // print hOSG.hTheta_CND minumum angle bin with non-zero content
    printf("OSG hTheta_CND min angle with content: %f deg\n", hOSG.hTheta_CND->GetXaxis()->GetBinLowEdge(hOSG.hTheta_CND->FindFirstBinAbove(0)));
    printf("CJ0 hTheta_CND min angle with content: %f deg\n", hCJ0.hTheta_CND->GetXaxis()->GetBinLowEdge(hCJ0.hTheta_CND->FindFirstBinAbove(0)));
    printf("CJ1 hTheta_CND min angle with content: %f deg\n", hCJ1.hTheta_CND->GetXaxis()->GetBinLowEdge(hCJ1.hTheta_CND->FindFirstBinAbove(0)));
    printf("CJ2 hTheta_CND min angle with content: %f deg\n", hCJ2.hTheta_CND->GetXaxis()->GetBinLowEdge(hCJ2.hTheta_CND->FindFirstBinAbove(0)));

    {
        auto *c = new TCanvas("c_diff", "OSG - CJ bin-by-bin differences and originals", 1800, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hEnergy_CND, hCJ0.hEnergy_CND, hCJ1.hEnergy_CND, hCJ2.hEnergy_CND, dE_CJ0, dE_CJ1, dE_CJ2}, {"OSG", "CJ0", "CJ1", "CJ2", "OSG - CJ0", "OSG - CJ1", "OSG - CJ2"}, /*normalize=*/false, /*logy=*/true, /*legend box*/ 0.12, 0.75, 0.55, 0.90, /*ncols=*/2, "CND E difference (bin-by-bin);E [GeV];counts");

        c->cd(2);
        draw_overlay_1D_N({hOSG.hTheta_CND, hCJ0.hTheta_CND, hCJ1.hTheta_CND, hCJ2.hTheta_CND, dTh_CJ0, dTh_CJ1, dTh_CJ2}, {"OSG", "CJ0", "CJ1", "CJ2", "OSG - CJ0", "OSG - CJ1", "OSG - CJ2"}, /*normalize=*/false, /*logy=*/true, /*legend box*/ 0.12, 0.75, 0.55, 0.90, /*ncols=*/2, "CND #theta difference (bin-by-bin);#theta [deg];counts");

        c->cd(3);
        draw_overlay_1D_N({hOSG.hPhi_CND, hCJ0.hPhi_CND, hCJ1.hPhi_CND, hCJ2.hPhi_CND, dPh_CJ0, dPh_CJ1, dPh_CJ2}, {"OSG", "CJ0", "CJ1", "CJ2", "OSG - CJ0", "OSG - CJ1", "OSG - CJ2"}, /*normalize=*/false, /*logy=*/true, /*legend box*/ 0.12, 0.75, 0.55, 0.90, /*ncols=*/2, "CND #phi difference (bin-by-bin);#phi [deg];counts");

        c->SaveAs("Histograms/cmp_CND_scint_diff_binbybin_all.png");
    }
    // ------------------------------------------------------------
    // 1) REC::Particle neutrons: p, theta, phi
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_particle_1D", "REC::Particle neutrons and MC::Particle", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hP, hCJ0.hP, hCJ1.hP, hCJ2.hP, hOSG.hPMC, hCJ0.hPMC, hCJ1.hPMC, hCJ2.hPMC}, {"OSG REC::Particle", "CJ0 REC::Particle", "CJ1 REC::Particle", "CJ2 REC::Particle", "OSG MC::Particle", "CJ0 MC::Particle", "CJ1 MC::Particle", "CJ2 MC::Particle"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.12, 0.75, 0.83, 0.90, /*ncols=*/2, "Neutron p from REC::Particle and MC::Particle;p [GeV];Counts");
        c->cd(2);
        draw_overlay_1D_N({hOSG.hTheta, hCJ0.hTheta, hCJ1.hTheta, hCJ2.hTheta, hOSG.hThetaMC, hCJ0.hThetaMC, hCJ1.hThetaMC, hCJ2.hThetaMC}, {"OSG REC::Particle", "CJ0 REC::Particle", "CJ1 REC::Particle", "CJ2 REC::Particle", "OSG MC::Particle", "CJ0 MC::Particle", "CJ1 MC::Particle", "CJ2 MC::Particle"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.12, 0.75, 0.83, 0.90, /*ncols=*/2, "Neutron #theta from REC::Particle and MC::Particle;#theta [deg];Counts");
        c->cd(3);
        draw_overlay_1D_N({hOSG.hPhi, hCJ0.hPhi, hCJ1.hPhi, hCJ2.hPhi, hOSG.hPhiMC, hCJ0.hPhiMC, hCJ1.hPhiMC, hCJ2.hPhiMC}, {"OSG REC::Particle", "CJ0 REC::Particle", "CJ1 REC::Particle", "CJ2 REC::Particle", "OSG MC::Particle", "CJ0 MC::Particle", "CJ1 MC::Particle", "CJ2 MC::Particle"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.12, 0.75, 0.83, 0.90, /*ncols=*/2, "Neutron #phi from REC::Particle and MC::Particle;#phi [deg];Counts");
        c->SaveAs("Histograms/cmp_RECParticle_MCParticle_neutrons_1D.png");
    }

    // ------------------------------------------------------------
    // 2) CND-matched scintillator info: energy, theta_cluster, phi_cluster
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_cndcluster_1D", "CND matched scintillator (REC::Scintillator, detector==CND)", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hEnergy_CND, hCJ0.hEnergy_CND, hCJ1.hEnergy_CND, hCJ2.hEnergy_CND}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1, "Neutron E from CND-cluster;E [GeV];Counts");
        c->cd(2);
        draw_overlay_1D_N({hOSG.hTheta_CND, hCJ0.hTheta_CND, hCJ1.hTheta_CND, hCJ2.hTheta_CND}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1, "Neutron #theta from CND-cluster;#theta [deg];Counts");
        c->cd(3);
        draw_overlay_1D_N({hOSG.hPhi_CND, hCJ0.hPhi_CND, hCJ1.hPhi_CND, hCJ2.hPhi_CND}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/true, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1, "Neutron #phi from CND-cluster;#phi [deg];Counts");
        c->SaveAs("Histograms/cmp_CND_scint_1D.png");
    }

    // ------------------------------------------------------------
    // 3) Residuals: Δtheta, Δphi  (cluster - particle)
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_residuals_1D", "CND residuals", 1500, 450);
        c->Divide(2, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hDTheta_CND, hCJ0.hDTheta_CND, hCJ1.hDTheta_CND, hCJ2.hDTheta_CND}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/false, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1);
        c->cd(2);
        draw_overlay_1D_N({hOSG.hDPhi_CND, hCJ0.hDPhi_CND, hCJ1.hDPhi_CND, hCJ2.hDPhi_CND}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/false, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1);

        c->SaveAs("Histograms/cmp_CND_residuals_1D.png");
    }

    // ------------------------------------------------------------
    // 4) Residuals: Δp, Δtheta, Δphi (particle(REC) - particle(MC))
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_residuals_MC_1D", "CND residuals MC", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hDP_REC_MC, hCJ0.hDP_REC_MC, hCJ1.hDP_REC_MC, hCJ2.hDP_REC_MC}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/false, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1);
        c->cd(2);
        draw_overlay_1D_N({hOSG.hDTheta_REC_MC, hCJ0.hDTheta_REC_MC, hCJ1.hDTheta_REC_MC, hCJ2.hDTheta_REC_MC}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/false, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1);
        c->cd(3);
        draw_overlay_1D_N({hOSG.hDPhi_REC_MC, hCJ0.hDPhi_REC_MC, hCJ1.hDPhi_REC_MC, hCJ2.hDPhi_REC_MC}, {"OSG", "CJ0", "CJ1", "CJ2"}, /*normalize =*/false, /*logy =*/false, /*legend box*/ 0.15, 0.72, 0.35, 0.90, /*ncols=*/1);
        c->SaveAs("Histograms/cmp_residuals_MC_1D.png");
    }

    // ------------------------------------------------------------
    // 5) CND layer for matched neutrons
    // ------------------------------------------------------------
    {
        auto *c = new TCanvas("c_layer_1D", "CND layer", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hOSG.hCNDLayer, hCJ0.hCNDLayer, hCJ1.hCNDLayer, hCJ2.hCNDLayer},
                          {"OSG", "CJ0", "CJ1", "CJ2"},
                          /*normalize=*/false,
                          /*logy=*/true,
                          0.15, 0.72, 0.35, 0.90, 1,
                          "CND layer for matched neutrons;layer;Counts");

        c->cd(2);
        draw_overlay_1D_N({hOSG.hCNDLayer_occupancy, hCJ0.hCNDLayer_occupancy, hCJ1.hCNDLayer_occupancy, hCJ2.hCNDLayer_occupancy},
                          {"OSG", "CJ0", "CJ1", "CJ2"},
                          /*normalize=*/false,
                          /*logy=*/true,
                          0.15, 0.72, 0.35, 0.90, 1,
                          "CND layer occupancy for matched neutrons;layer;Counts");
        c->cd(3);
        draw_overlay_1D_N({hOSG.hCND_NLayers, hCJ0.hCND_NLayers, hCJ1.hCND_NLayers, hCJ2.hCND_NLayers},
                          {"OSG", "CJ0", "CJ1", "CJ2"},
                          /*normalize=*/false,
                          /*logy=*/true,
                          0.15, 0.72, 0.35, 0.90, 1,
                          "Number of CND layers hit per neutron;layers;Counts");
        c->SaveAs("Histograms/cmp_CND_layer_1D.png");
    }

    // ------------------------------------------------------------
    // 6) CND distributions for OSG events missing in CJ0/CJ1
    // ------------------------------------------------------------
    auto keysOSG = collect_keys(chOSG, maxEvents);

    Hists hCJ0_inOSG = book_hists("CJ0_inOSG");
    Hists hCJ1_inOSG = book_hists("CJ1_inOSG");

    process_chain_allow_keys(chCJ0, hCJ0_inOSG, "CJ0_inOSG", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, &thCut, keysOSG);
    process_chain_allow_keys(chCJ1, hCJ1_inOSG, "CJ1_inOSG", maxEvents, /*cnd_id*/ 3, /*maxAngleDeg=*/10, &thCut, keysOSG);

    {
        auto *c = new TCanvas("c_missing_events_cnd",
                              "CND distributions for OSG events missing in CJ0/CJ1", 1500, 450);
        c->Divide(3, 1);

        c->cd(1);
        draw_overlay_1D_N({hCJ0_inOSG.hEnergy_CND, hCJ1_inOSG.hEnergy_CND},
                          {"OSG events missing in CJ0", "OSG events missing in CJ1"},
                          /*normalize=*/false, /*logy=*/true,
                          0.15, 0.72, 0.55, 0.90, 1,
                          "CND E for OSG events missing in CJ0/CJ1;E [GeV];Counts");

        c->cd(2);
        draw_overlay_1D_N({hCJ0_inOSG.hTheta_CND, hCJ1_inOSG.hTheta_CND},
                          {"OSG events missing in CJ0", "OSG events missing in CJ1"},
                          /*normalize=*/false, /*logy=*/true,
                          0.15, 0.72, 0.55, 0.90, 1,
                          "CND #theta for OSG events missing in CJ0/CJ1;#theta [deg];Counts");

        c->cd(3);
        draw_overlay_1D_N({hCJ0_inOSG.hPhi_CND, hCJ1_inOSG.hPhi_CND},
                          {"OSG events missing in CJ0", "OSG events missing in CJ1"},
                          /*normalize=*/false, /*logy=*/true,
                          0.15, 0.72, 0.55, 0.90, 1,
                          "CND #phi for OSG events missing in CJ0/CJ1;#phi [deg];Counts");

        c->SaveAs("Histograms/cmp_CND_missing_events_1D.png");
    }

    // ------------------------------------------------------------
    // 7) 2D histograms (triptychs)
    // ------------------------------------------------------------
    draw_triptych_2D(hOSG.hPTheta, hCJ0.hPTheta, hCJ1.hPTheta, hCJ2.hPTheta,
                     "p vs theta (REC::Particle)", "Histograms/cmp_2D_p_vs_theta.png",
                     /*normalize =*/false, /*logz =*/true);

    draw_triptych_2D(hOSG.hPPhi, hCJ0.hPPhi, hCJ1.hPPhi, hCJ2.hPPhi,
                     "p vs phi (REC::Particle)", "Histograms/cmp_2D_p_vs_phi.png",
                     /*normalize =*/false, /*logz =*/true);

    draw_triptych_2D(hOSG.hThetaPhi, hCJ0.hThetaPhi, hCJ1.hThetaPhi, hCJ2.hThetaPhi,
                     "theta vs phi (REC::Particle)", "Histograms/cmp_2D_theta_vs_phi.png",
                     /*normalize =*/false, /*logz =*/true);

    // ------------------------------------------------------------
    // 8) Save all histograms into a ROOT file
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

        fout->mkdir("CJ2");
        fout->cd("CJ2");
        write_hists(hCJ2);

        fout->mkdir("CJ0_inOSG");
        fout->cd("CJ0_inOSG");
        write_hists(hCJ0_inOSG);

        fout->mkdir("CJ1_inOSG");
        fout->cd("CJ1_inOSG");
        write_hists(hCJ1_inOSG);

        fout->mkdir("DIFF");
        fout->cd("DIFF");
        if (dE_CJ0)
            dE_CJ0->Write();
        if (dE_CJ1)
            dE_CJ1->Write();
        if (dE_CJ2)
            dE_CJ2->Write();
        if (dTh_CJ0)
            dTh_CJ0->Write();
        if (dTh_CJ1)
            dTh_CJ1->Write();
        if (dTh_CJ2)
            dTh_CJ2->Write();
        if (dPh_CJ0)
            dPh_CJ0->Write();
        if (dPh_CJ1)
            dPh_CJ1->Write();
        if (dPh_CJ2)
            dPh_CJ2->Write();

        fout->Write();

        // std::cout << "hPhi_OSG entries = " << hOSG.hPhi->GetEntries()
        //           << "  bin(0) = " << hOSG.hPhi->GetBinContent(hOSG.hPhi->FindBin(0.0))
        //           << std::endl;

        // int bx = hOSG.hPPhi->GetXaxis()->FindBin(1e-6); // near p=0
        // int by = hOSG.hPPhi->GetYaxis()->FindBin(0.0);  // phi=0
        // std::cout << "2D bin(p~0,phi=0) = " << hOSG.hPPhi->GetBinContent(bx, by) << "\n";
        // std::cout << "2D integral = " << hOSG.hPPhi->Integral() << "\n";

        print_counts("OSG", hOSG);
        print_counts("CJ0", hCJ0);
        print_counts("CJ1", hCJ1);
        print_counts("CJ2", hCJ2);
        fout->Close();
    }
    std::cout << "Done. Wrote PNGs + cmp_cnd_versions_hists.root" << std::endl;
}