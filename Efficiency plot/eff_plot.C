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

void eff_plot()
{

   TCanvas *c1 = new TCanvas("c3", "c3",72,64,1188,673);
   TMultiGraph  *mg  = new TMultiGraph(); 
   float e1[3], e2[8], de1[3], de2[8], dp1[3] = {0}, dp2[8] = {0};                                     //10,10    5,10       5,5          5,8,           4,10         8, 10
   //float p1[3] = {0}; 
   float p2[8] = {0};
     
    // p1[0] = 10.0; p1[1] = 5.0;  p1[2] = 4.0; //  keeping zoomed parameter constant at 10
    // p2[0] = 10.0; p2[1] = 8.0;  p2[2] = 7.0; p2[3] = 6.0; p2[4] = 5.0;  // keeping unzoomed parameter constant at 5
     
    // e1[0] = 0.944; e1[1] = 0.946; e1[2] = 0.945;
    // e2[0] = 0.948; e2[1] = 0.946; e2[2] = 0.946; e2[3] = 0.933; e2[4] = 0.919;
     
    // de1[0] = 0.00407; de1[1] = 0.004002; de1[2] = 0.004037;
    // de2[0] = 0.00393; de2[1] = 0.00400; de2[2] = 0.00400; de2[3] = 0.00442; de2[4] = 0.00483;
    
     //p2[0] = 6.0; p2[1] = 12.0;  
     p2[0] = 19.0; p2[1] = 25.0; p2[2] = 32.0, p2[3] = 38.0, p2[4] = 45.0, p2[5] = 51.0, p2[6] = 58.0, p2[7] = 64.0;  
   //e2[0] = 0.946; e2[1] = 0.942; 
   e2[0] = 0.430; e2[1] = 0.448; e2[2] = 0.469, e2[3] = 0.496, e2[4] = 0.534, e2[5] = 0.571, e2[6] = 0.621, e2[7] = 0.664;
  // de2[0] = 0.00400; de2[1] = 0.00414; 
   de2[0] = 0.00239; de2[1] = 0.00531; de2[2] = 0.00464; de2[3] = 0.00809; de2[4] = 0.00694; de2[5] = 0.00604; de2[6] = 0.01597; de2[7] = 0.00661;
    
  //   e2[0] = 0.410; e2[1] = 0.417; e2[2] = 0.433; e2[3] = 0.452; e2[4] = 0.496, e2[5] = 0.510, e2[6] = 0.545, e2[7] = 0.572, e2[8] = 0.642, e2[9] = 0.671;
 //    de2[0] = 0.0; de2[1] = 0.0; de2[2] = 0.0; de2[3] = 0.0; de2[4] = 0.0, de2[5] = 0.0, de2[6] = 0.0, de2[7] = 0.0, de2[8] = 0.0, de2[9] = 0.0;
    
     
   //  TGraphErrors *gr1; 
   //  gr1 = new TGraphErrors(3, p1, e1, dp1, de1);
   //  gr1->SetLineColor(4);   
  //   mg->Add(gr1);
     
     TGraphErrors *gr2;
     gr2 = new TGraphErrors(8, p2, e2, dp2, de2);
     gr2->SetLineColor(kBlue+2);   
     mg->Add(gr2);
 
     c1->cd();
     mg->SetTitle("Algorithm runtime with increasing noise hits");
     mg->GetXaxis()->SetTitle("Number of noise hits");
     mg->GetYaxis()->SetTitle("Algorithm runtime (s/event)");
     mg->GetYaxis()->SetLimits(0.40,0.70);
     //mg->GetXaxis()->SetLimits(-1000, 1000);
     mg->Draw("A*");
     
    
}
