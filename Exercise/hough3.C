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

int hough3()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);

     TRandom* Uni= new TRandom();
     TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
     TCanvas *c2 = new TCanvas("c2", "c1",72,64,1188,673);
     TCanvas *c3 = new TCanvas("c3", "c3",72,64,1188,673);
     TH2D* h1 = new TH2D("h1", "Accumulator; #theta; #rho", 36, 0.0, 180.0, 36, -5.0, 9.0);
     
     
     TMultiGraph  *mg  = new TMultiGraph();
     
     
     
     const int c = 5, r = 36, theta_max = 180, rho_max = 10;
     double_t accumulator[theta_max][rho_max];
     double_t x[c], y[c], rho[r], theta[r], m = 1.0, b = 2.0, rh, th;
     int k = 1;   
     for(int i = 0; i< c; i++)
     {
          x[i] = 1 * k;
          y[i] = m*x[i] + b;
           
           k++;
     
     }
         
     c1->cd();
     TGraph* gr1 = new TGraph(c,x,y);
     gr1->SetTitle("Coordinate Space");
     gr1->GetXaxis()->SetTitle("x axis");
     gr1->GetYaxis()->SetTitle("y axis");
     gr1->Draw("AP*");

       theta[0] = 0.0;
       double h = (double) theta_max/r;   //  theta_max = 180
         
         
   
     for(int j = 0; j < c; j ++)
     {   
      rho[0] = x[j] * cos(theta[0]* M_PI/180.0) + y[j] * sin(theta[0]* M_PI/180.0);
        
      for(int l = 1; l < r; l ++)
     {     
           
         theta[l] = theta[l-1] + h;
         rho[l] = x[j] * cos(theta[l]* M_PI/180.0) + y[j] * sin(theta[l]* M_PI/180.0);
         
         th = theta[l];
         rh = rho[l];      
         h1->Fill(th,rh);
     }
     
     TGraph* gr2 = new TGraph(r,theta,rho);   
     mg->Add(gr2);
    
    }
     
     
     c2->cd();
     mg->SetTitle("Parameter Space");
     mg->GetXaxis()->SetTitle("#theta");
     mg->GetYaxis()->SetTitle("#rho");
     mg->Draw("AC");
     
     c3->cd();
     h1->SetStats(0);
     h1->DrawCopy("COLZ");
     
     
     return 0;
     
 }
 
 
 
