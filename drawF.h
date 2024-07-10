//
// Created by nikol on 2/21/2024.
//

#ifndef PYTHIAPROJECT_DRAWF_H
#define PYTHIAPROJECT_DRAWF_H

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






void drawParticles_histogram(std::vector<Pythia8::Particle> & particles_histogram, double minpT);
TH2D * createTH2D(int nXBins, int nYBins, double nXMax);
void drawdrawLegend();
void setUpRootStyle();
void drawText(double x, double y, TString txt, int align= 11, double tsize= 0.032);
void drawParticleMarker(const Pythia8::Particle &p, int style, int col, double size= 1.0);
void drawParticleText(const Pythia8::Particle &p, int colourHS);
void drawLegendBox(double x1, double y1, double x2, double y2);
void drawMarker(double x, double y, int style, int col, double size= 1.0);

#endif //PYTHIAPROJECT_DRAWF_H
