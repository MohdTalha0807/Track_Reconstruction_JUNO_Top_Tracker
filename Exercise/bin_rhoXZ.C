//#include "binsize_theta_x.h"

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
#include <cmath>



void bin_rhoXZ()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);
     TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
     
     Double_t strip = 2.64;
     const int c = 720;
     Double_t  min, theta_radian, theta_degree, thx_radian, thx_degree, thx_phi, theta_phi; 
     std::vector<float> delta_rhoXZ;
     std::vector<float> theta_x; 
     int cntr =0;
     
        
        
      Double_t z = 150.0;
      Double_t theta[c];
      theta[0] = -90.0;
      
      for(int i =0; i < 720; i++)
      {  
            
         Double_t strip = 2.64;
         Double_t delta_rho = strip * cos(theta[i] * M_PI/180.0); 
   
         
         delta_rhoXZ.insert(delta_rhoXZ.begin() + i, delta_rho);
         theta_x.insert(theta_x.begin() + i, theta[i]);
         theta[i+1] = theta[i] + 0.25;
         
        
      }
        
     
       cout<<"delta_thetaXZ[l]"<<"\t"<<"theta_x"<<endl;
       cout<<"----------------------------------------------------"<<endl;
       
       for(int l = 0; l < delta_rhoXZ.size(); l++)
       {
         cout<<delta_rhoXZ[l]<<"\t"<<theta_x[l]<<endl;
         
       }
       min = *min_element(delta_rhoXZ.begin(), delta_rhoXZ.end());
       
       cout<<"----------------------------------------"<<endl;
       cout<<"The minimum value of the delta theta for XZ projection is : "<<min<<endl;
       cout<<"----------------------------------------"<<endl;
       
       int dim = theta_x.size();
       
       c1->cd();
       gPad->SetGrid();
       TGraph *gr = new TGraph(dim, &theta_x[0], &delta_rhoXZ[0]);
       gr->SetTitle("Delta theta vs theta for XZ projection");
       gr->GetXaxis()->SetTitle("#theta (degrees)");
       gr->GetYaxis()->SetTitle("#delta#theta (degrees)");
       gr->GetYaxis()->SetLimits(-2.64,2.64);
       gr->GetXaxis()->SetLimits(-90.0, 90.0);
       gr->Draw("AC");
       
      
        
     return 0;
     
 }
 

 
