#include "drawF.h"



void drawParticles_histogram(std::vector<Pythia8::Particle> & particles_histogram, double minpT){
    int colPos = kRed, colNeg = kBlue, colNeut = kGreen + 3;
    double yMax = 4;
    for (auto &p: particles_histogram) {
        if (!( std::abs(p.y()) < yMax && p.pT() > minpT )) continue; //gets things only in canvas
        if (p.charge() > 0) {
            drawParticleMarker(p, 5, colPos, 0.8);
        } else if (p.charge() < 0) {
            drawParticleMarker(p, 5, colNeg, 0.8);
        } else {
            drawParticleMarker(p, 21, colNeut, 0.4);
            drawParticleMarker(p, 5, colNeut, 0.8);
        }
    }
}

TH2D * createTH2D(int nXBins, int nYBins, double nXMax){ //don't forget to free the memory
    auto result = new TH2D("", ";Rapidity #it{y};Azimuth #it{#phi};Jet #it{p}_{T} [GeV]",
                           nXBins, -nXMax, nXMax, nYBins, -TMath::Pi(), TMath::Pi());

    result->GetYaxis()->SetTitleOffset(0.5);//offset of y-axis label "Phi"
    result->GetZaxis()->SetTitleOffset(1.3); //offset of z-axis label "pT"
    return result;

}

void drawdrawLegend() {
    int colPos = kRed, colNeg = kBlue, colNeut = kGreen + 3;
    double x = 0.67, y = 0.814, dx = 0.19, dy = 0.126; //parameters for the size of the legend box
    double dxt = 0.001, yl = 0.92, dyl =0.03;  //parameters for text

    drawLegendBox(x, y, x + dx, y + dy);

    drawText(x + dxt, yl, "Stable particles", 12);
    drawText(x + dxt, yl - dyl, "    #bf{#minus}    #scale[0.9]{posetive}", 12);
    drawText(x + dxt, yl - (2*dyl), "    #bf{#minus}    #scale[0.9]{neutral}", 12);
    drawText(x + dxt, yl - (3*dyl), "    #bf{#minus}    #scale[0.9]{negative}", 12);
    drawMarker(0.685, yl - dyl, 5, colPos, 0.8);
    drawMarker(0.685, yl - (2*dyl), 21, colNeut, 0.4);
    drawMarker(0.685, yl - (2*dyl), 5, colNeut, 0.8);
    drawMarker(0.685, yl - (3*dyl), 5, colNeg, 0.8);
}

void setUpRootStyle() {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02, "x");
    gStyle->SetTickLength(0.015, "y");
    gStyle->SetPalette(55);
}

void drawText(double x, double y, TString txt, int align,
              double tsize) {
    static auto tex = new TLatex();
    tex->SetTextAlign(align);
    tex->SetTextSize(tsize);
    tex->SetTextFont(42);
    tex->SetNDC();
    tex->DrawLatex(x, y, txt);
}

//==========================================================================
// Text to draw a marker at the (y, phi) coordinates of a particle.
// Absolute coordinates.

void drawParticleMarker(const Pythia8::Particle &p, int style, int col,
                        double size) {
    static auto m = new TMarker();
    m->SetMarkerStyle(style);
    m->SetMarkerSize(size);
    m->SetMarkerColor(col);
    m->DrawMarker(p.y(), p.phi());
}

//==========================================================================
// Method to draw a marker+text of a particle.

void drawParticleText(const Pythia8::Particle &p, int colourHS) {
    // Draws a marker at (y, phi) of particle. Circle for parton, star
    // for boson.
    bool isParton = (std::abs(p.id()) <= 5 || p.id() == 21);
    int col = colourHS;
    drawParticleMarker(p, isParton ? 20 : 29, col, isParton ? 0.8 : 1.2);

    // Format the name-string of the particle according to ROOT's TLatex.
    // Print the text right under the marker.
    TString name = p.name();
    if (name.Contains("bar")) name = "#bar{" + name.ReplaceAll("bar", "") + "}";
    name.ReplaceAll("+", "^{+}").ReplaceAll("-", "^{-}").ReplaceAll("h0", "H");
    static auto tex = new TLatex();
    tex->SetTextSize(0.03);
    tex->SetTextFont(42);
    tex->SetTextAlign(11);
    tex->SetTextColor(col);
    tex->DrawLatex(p.y() + 0.1, p.phi() - 0.1, "#it{" + name + "}");
}

//==========================================================================
// Draws a box for text to appear.

void drawLegendBox(double x1, double y1, double x2, double y2) {
    static auto *box = new TPave(x1, y1, x2, y2, 1, "ndc");
    box->SetFillColor(kWhite);
    box->Draw();
}

//==========================================================================
// Draw a marker for legend.

void drawMarker(double x, double y, int style, int col, double size) {
    auto m = new TMarker(x, y, style);
    m->SetMarkerSize(size);
    m->SetMarkerColor(col);
    m->SetNDC(true);
    m->Draw();
}
