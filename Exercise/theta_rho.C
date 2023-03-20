#include "Math/Integrator.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TPave.h"
void theta_rho()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);

     TCanvas *c1 = new TCanvas("c1", "XZ",72,64,1188,673);
     TH2D* h1 = new TH2D("h1", "Accumulator; #theta; #rho", 180, -90.0, 90.0, 100, -2000.0, 2000.0);
     TCanvas *c2 = new TCanvas("c2", "YZ",72,64,1188,673);
     TH2D* h2 = new TH2D("h1", "Accumulator; #theta; #rho", 180, -90.0, 90.0, 50, -1000.0, 1000.0);
          
     
     const int c = 6, r = 180, theta_max = 180, rho_max = 2000;
     double_t x[c], z[c], y[c], rho[r], theta[r], rh, th;
     
     x[0] = -2000.0, x[1] = 0.0, x[2] = 2000.0;
     z[0] =  150.0,  z[1] = 0.0, z[2] = -150.0;
     x[3] = 2000.0, x[4] = 0.0, x[5] = -2000.0;
     z[3] =  150.0,  z[4] = 0.0, z[5] = -150.0;
    
     y[0] = 1000.0, y[1] = 0.0, y[2] = -1000.0;
     y[3] = -1000.0, y[4] = 0.0, y[5] = 1000.0;
     
  /*-------------------------for XZ case -----------------------*/    
          
       theta[0] = -90.0;
       double h = (double) theta_max/r;   //  theta_max = 180
    
   
     for(int j = 0; j < c; j ++)
     {           
      for(int l = 0; l < r; l ++)
     {     
         rh = x[j] * cos(theta[l]* M_PI/180.0) + z[j] * sin(theta[l]* M_PI/180.0);    
         th = theta[l];
         
         h1->Fill(th,rh);
         theta[l+1] = theta[l] + h;
     }
    }
     
     c1->cd();
     h1->SetStats(0);
     gStyle->SetPalette(1);
     h1->Draw("COLZ");
     gPad->Update();
     TPaletteAxis *palette1 = (TPaletteAxis*)h1->GetListOfFunctions()->FindObject("palette");
     palette1->SetY2NDC(0.7);
     
   
 /*-------------------------for YZ case -----------------------*/
    
    for(int j = 0; j < c; j ++)
     {           
      for(int l = 0; l < r; l ++)
     {     
         rh = y[j] * cos(theta[l]* M_PI/180.0) + z[j] * sin(theta[l]* M_PI/180.0);    
         th = theta[l];
         
         h2->Fill(th,rh);
         theta[l+1] = theta[l] + h;
     }
    }
     
     c2->cd();
     h2->SetStats(0);
     //gStyle->SetPalette(1);
     h2->Draw("COLZ");
     gPad->Update();
     TPaletteAxis *palette2 = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
     palette2->SetY2NDC(0.7);
     
     return 0;
     
 }
 
 
 
