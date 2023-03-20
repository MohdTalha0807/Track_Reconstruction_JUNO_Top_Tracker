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

int hough1()
{    

     TRandom* Uni= new TRandom();
     TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
     TCanvas *c2 = new TCanvas("c2", "c1",72,64,1188,673);
     TCanvas *c3 = new TCanvas("c3", "c3",72,64,1188,673);
     TH1D* h1 = new TH1D("h1", "Hough Space; theta;rho;", 100, 0.0, 2.0);
     TH1D *h2[20];
     
     TMultiGraph  *mg  = new TMultiGraph();
     
     //double_t accumulator[slope][12];
     
     const int c = 5, r = 10;
     double_t x[c], y[c], intercept[r], slope[r], m = 1.0, b = 2.0;
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

       slope[0] = 0.0;
       double h = (double) 2.0/10.0;
       cout<<h<<endl;  
         
   /*     
     for(int j = 0; j < c; j ++)
     {   
      for(int k = 0; k < 36; k ++)
     {    
       
         rho[k] = x[j] * cos(theta[k]) + y[j] * sin(theta[k]);
                 
         theta[++k] = theta[k] + h;
     
     }
     
    }
     
   */  
     
     for(int j = 0; j < c; j ++)
     {   
      
      intercept[0] = -1.0 * slope[0]*x[j] + y[j];
      for(int l = 1; l < r; l ++)
     {     
           
         slope[l] = slope[l-1] + h;
         intercept[l] = -1.0*slope[l]*x[j] + y[j];
                
     }
     
     TGraph* gr2 = new TGraph(r,slope,intercept);   
     mg->Add(gr2);
    }
     
     
     c2->cd();
     mg->SetTitle("Parameter Space");
     mg->GetXaxis()->SetTitle("Slope");
     mg->GetYaxis()->SetTitle("Intercept");
     mg->Draw("AC");
     
     
     return 0;
     
 }
