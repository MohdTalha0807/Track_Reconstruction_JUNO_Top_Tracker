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



void bin_thetaXZ()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);
     TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
     
     Double_t strip = 2.64;
     Double_t numerator, denominator, th =0.0;
     const int c = 720;
     Double_t  min, theta_radian, theta_degree, thx_radian, thx_degree, thx_phi, theta_phi; 
     std::vector<float> delta_thetaXZ;
     std::vector<float> theta_x; 
     int cntr =0;
     
        
        
      Double_t z = 150.0;
      Double_t theta[c];
      theta[0] = -90.0;
      
      for(int i =0; i < 720; i++)
      {  
         
         Double_t rho1 = sqrt( z*z + (pow((z/tan((90.0 - theta[i]) * M_PI/180.0)),2)));
         Double_t rho2 = sqrt( z*z + (pow((z/tan((90.0 - theta[i]) * M_PI/180.0)) + strip,2))); 
         Double_t delta_theta = acos((rho1*rho1 + rho2*rho2 - strip*strip)/(2.0 * rho1 * rho2));
     
         theta_degree = delta_theta * 180.0/M_PI;
         
         
         delta_thetaXZ.insert(delta_thetaXZ.begin() + i, theta_degree);
         theta_x.insert(theta_x.begin() + i, theta[i]);
         theta[i+1] = theta[i] + 0.25;
         
        
      }
        
     
       cout<<"delta_thetaXZ[l]"<<"\t"<<"theta_x"<<endl;
       cout<<"----------------------------------------------------"<<endl;
       
       for(int l = 0; l < delta_thetaXZ.size(); l++)
       {
         cout<<delta_thetaXZ[l]<<"\t"<<theta_x[l]<<endl;
         
       }
       min = *min_element(delta_thetaXZ.begin(), delta_thetaXZ.end());
       
       cout<<"----------------------------------------"<<endl;
       cout<<"The minimum value of the delta theta for XZ projection is : "<<min<<endl;
       cout<<"----------------------------------------"<<endl;
       
       int dim = theta_x.size();
       
       c1->cd();
       gPad->SetGrid();
       TGraph *gr = new TGraph(dim, &theta_x[0], &delta_thetaXZ[0]);
       gr->SetTitle("Delta theta vs theta for XZ projection");
       gr->GetXaxis()->SetTitle("#theta (degrees)");
       gr->GetYaxis()->SetTitle("#delta#theta (degrees)");
       gr->GetYaxis()->SetLimits(-2.64,2.64);
       gr->GetXaxis()->SetLimits(-90.0, 90.0);
       gr->Draw("AC");
       
        delta_thetaXZ.insert(delta_thetaXZ.end(), theta_x.begin(), theta_x.end());
        
     return 0;
     
 }
 

 
