//#include "binsize_rho_x.h"

#include "TROOT.h"
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
#include <vector>
#include <algorithm>


vector<float> binsize_rho_x()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);
     TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);

     
     const int c = 2, r = 720, theta_max = 90, rho_max = 2000;
     float x[c], z[c], theta[r], rh1, rh2, th, drho, min;    
     std::vector<float> delta_rhoXZ;
     std::vector<float> theta_x; 
     int cntr = 0;
          x[0] = 1988.13; x[1] = 1990.77;
          z[0] = 153.785; z[1] = 153.785;
           
       theta[0] = -90.0;
       double h = (double) 180/r;   //  theta_max = 180
         
       
     for(int l = 0; l <= r; l++)
      {    
         rh1 = x[0] * cos(theta[l]* M_PI/180.0) + z[0] * sin(theta[l]* M_PI/180.0);
         rh2 = x[1] * cos(theta[l]* M_PI/180.0) + z[1] * sin(theta[l]* M_PI/180.0);
         
         th = theta[l];
         drho = abs(rh1 - rh2);
      
         if(th >= -1*86.0 && th <= 86.0)
       { 
         delta_rhoXZ.insert(delta_rhoXZ.begin() + cntr, drho);
         cntr++;      
        }
         
       theta[l+1] = theta[l] + h; 
      }
    
 /*      for(int l = 0; l < delta_rho.size(); l++)
       {
         cout<<delta_rho[l]<<endl;
         
       }
       min = *min_element(delta_rho.begin(), delta_rho.end());
       
       cout<<"----------------------------------------"<<endl;
       cout<<"The minimum value of the vector XZ is : "<<min<<endl;
       cout<<"----------------------------------------"<<endl;
   */    
       int ic = 0;
       
        for(int i = 0; i<=r; i++)
       {  
         if(theta[i] >= -1*86.0 && theta[i] <= 86.0)
         {
          theta_x.insert(theta_x.begin() + ic, theta[i]);
          ic++;
         } 
       }
       
      // delta_rhoXZ.insert(delta_rhoXZ.end(), theta_x.begin(), theta_x.end());
   
   
       int dim = theta_x.size();
      // gPad->SetGrid();
       TGraph *gr = new TGraph(dim, &theta_x[0], &delta_rhoXZ[0]);
      // gr->SetTitle("Delta Rho vs theta for XZ projection");
       gr->GetXaxis()->SetTitle("#theta (degrees)");
       gr->GetYaxis()->SetTitle("#delta #rho (cm)");
       gr->GetYaxis()->SetLimits(0.0,2.64);
       gr->GetXaxis()->SetLimits(-90.0, 90.0);
       gr->SetLineColor(kBlue);
       gr->Draw("AC");
         
     return delta_rhoXZ;
     
 }
 
 
