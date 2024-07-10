#include <iostream>
#include <string>
#include <fstream>

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TPave.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

#include "drawF.h"

int main() {

    Pythia8::Pythia pythia; //Creates Pythia object
    pythia.readFile("../config1.cmnd"); //Reads the config file
    pythia.init(); //Inits the Pythia object

    //some parameters variables used for configuration
    int nXBins = 400 / 2, nYBins = 314 / 2; //resolutions of 2D histogram
    double nXMax = 4; //maximal rapidity value


    setUpRootStyle(); //Sets up the root style of the resulted histograms
    auto canvas = new TCanvas(); //creates a new canvas
    canvas->SetMargin(0.06, 0.02, 0.08, 0.06); //sets the margins of the canvas
    auto pTflow = createTH2D(nXBins, nYBins, nXMax); //creates a 2D histogram
    //creates array of 2D histograms


    TString pdf = "../results/"; //some string used to save the results - used as a file, where the results are saved

    std::map<TString, fastjet::JetDefinition> jetDefs; //map of jet definitions

    int i = 0; //iterartor
    //some parameters used for configuration
    double R = 0.3;
    double pTmin_jet = 5;
    double pTmin_hadron = 1, yMax = 4;
    TString description = "Number of events: " + std::to_string(pythia.mode("Main:numberOfEvents"));

    //define jet finding algorithms here:
    jetDefs["Anti-#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::antikt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);

    jetDefs["Anti-#it{k_{t}} jets, #it{R} = " + std::to_string(2 * R)] = fastjet::JetDefinition(
            fastjet::antikt_algorithm, 2 * R, fastjet::E_scheme, fastjet::Best);

    jetDefs["#it{k_{t}} jets, #it{R} = " + std::to_string(2 * R)] = fastjet::JetDefinition(
            fastjet::kt_algorithm, 2 * R, fastjet::E_scheme, fastjet::Best);
    //CaCambridge-Aachen example
    //        jetDefs["Cambridge-Aachen jets,  #it{R} = " + std::to_string(R)] = fastjet::JetDefinition(
    //                fastjet::cambridge_algorithm, R, fastjet::E_scheme, fastjet::Best);
    //till here


    TH2D *pTflows[jetDefs.size()];
    for (int i = 0; i < jetDefs.size(); i++) {
        pTflows[i] = createTH2D(nXBins, nYBins, nXMax);
    }

    auto &event = pythia.event;
    std::vector<Pythia8::Particle> particles_histogram; //vector of particles
    std::vector<Pythia8::Particle> particles_histogram_to_draw; //vector of particles
    std::vector<fastjet::PseudoJet> stable_particles; //vector of stable particles

    for (int iEvent = 0; iEvent < pythia.mode("Main:numberOfEvents"); ++iEvent) { //choosing final particles only
        if (!pythia.next()) continue;
        stable_particles.clear();
        particles_histogram.clear();

        for (int i = 0; i < event.size(); ++i) {
            auto &p = event[i];
            if (not p.isFinal()) continue;
            stable_particles.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
            //pushes the particle to the vector of stable particles
            particles_histogram.push_back(p); //pushes the particle to the vector of particles
            particles_histogram_to_draw.push_back(p);
        }

        //Ghost are needed otherwise jet images is bad or not possible to find
        fastjet::PseudoJet ghost; //creates a ghost particle with 0 momentum
        double pTghost = 1e-100;
        for (int iy = 1; iy <= 400 / 2; ++iy) {
            for (int iphi = 1; iphi <= 314 / 2; ++iphi) {
                double y = pTflows[0]->GetXaxis()->GetBinCenter(iy);
                double phi = pTflows[0]->GetYaxis()->GetBinCenter(iphi);
                ghost.reset_momentum_PtYPhiM(pTghost, y, phi, 0);
                stable_particles.push_back(ghost);
            }
        }


        canvas->SetLogz(); //log the z axis, so jets are more clearly seen
        canvas->SetRightMargin(0.14);

        i = 0;
        for (auto jetDef: jetDefs) {
            fastjet::ClusterSequence clustSeq(stable_particles, jetDef.second);
            auto jets = sorted_by_pt(clustSeq.inclusive_jets(pTmin_jet));
            // Fill the pT flow.
            // For each jet:

            for (auto jet: jets) {
                // For each particle:
                for (auto c: jet.constituents()) {
                    if (c.pt() > 1e-50) continue; //gets rid of bubbles in jets
                    pTflows[i]->Fill(c.rap(), c.phi_std(), jet.pt());
                }
            }
            ++i;
        }

    } //move it to the end in order to split events



    //here '}' must be added in order to split events
    //part of code to turn off hello notifications

    i = 0;
    for (auto jetDef: jetDefs) { //!Draw
        pTflows[i]->GetZaxis()->SetRangeUser(pTmin_jet / 4, pTflows[i]->GetBinContent(pTflows[i]->GetMaximumBin()) * 4);
        pTflows[i]->GetZaxis()->SetMoreLogLabels();
        pTflows[i]->Draw("colz");

        drawParticles_histogram(particles_histogram_to_draw, pTmin_hadron);

        drawText(0.06, 0.96, description);
        drawText(0.87, 0.96, jetDef.first +
                             Form(", #it{p}_{T} > %.0f GeV", pTmin_jet), 31);
        drawdrawLegend();
        canvas->Print(pdf + "[" + description + "] " + jetDef.first + ".pdf");;
        printf("Produced %s\n\n", pdf.Data());
        ++i;
    }

    delete pTflow;
    for (i = 0; i < jetDefs.size(); i++)
        delete pTflows[i];
    delete canvas;

    return 0;
}