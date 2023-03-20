// Draw a graph with text attached to each point.
// The text is drawn in a TExec function, therefore if the text is
// moved interactively, it will be automatically updated.
// Author: Olivier Couet

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <deque>
#include <cstring>
#include <algorithm>
#include <map>
#include <iterator>
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TColor.h"
#include "TApplication.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TRootCanvas.h"
#include "TExec.h"
#include "snprintf.h"
 

void graphtext() {
   TCanvas *c = new TCanvas("c","A Simple Graph Example with Text",700,500);
   c->SetGrid();
 
   const Int_t n = 10;
   auto gr = new TGraph(n);
   gr->SetTitle("A Simple Graph Example with Text");
   gr->SetMarkerStyle(20);
   auto ex = new TExec("ex","drawtext();");
   gr->GetListOfFunctions()->Add(ex);
 
   Double_t x, y;
   for (Int_t i=0;i<n;i++) {
      x = i*0.1;
      y = 10*sin(x+0.2);
      gr->SetPoint(i,x,y);
 
   }
   gr->Draw("ALP");
}
 
void drawtext()
{
   Int_t i,n;
   Double_t x,y;
   TLatex l;
 
   l.SetTextSize(0.025);
   l.SetTextFont(42);
   l.SetTextAlign(21);
   l.SetTextColor(kBlue);
 
   auto g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   n = g->GetN();
 
   for (i=0; i<n; i++) {
      g->GetPoint(i,x,y);
      l.PaintText(x,y+0.2,Form("(%4.2f,%4.2f)",x,y));
   }
}

