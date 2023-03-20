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



void binsize_theta_x()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);
     TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
     
     Double_t strip = 2.64;
     Double_t numerator, denominator, th =0.0;
     const int c = 4500;
     Double_t x[c], z, min, theta_radian, theta_degree, thx_radian, thx_degree, thx_phi, theta_phi; 
     std::vector<float> delta_thetaXZ;
     std::vector<float> theta_x; 
     int cntr =0;
     
          x[0] = -2250.0;
          z = 150;
        
      for(int i =0; i < 1704; i++)
      {
         numerator   = pow(x[i],2) + z*z + x[i]*strip;                              
         denominator = sqrt(pow(x[i],2) + z*z) * sqrt(pow(x[i] + strip,2) + z*z);
      
         th = numerator/denominator;
    
         thx_radian = atan(z/x[i]);
         thx_degree = thx_radian * 180.0/M_PI;
         if(thx_degree > 0)
         {
           thx_phi = 90.0 - thx_degree;
          }
         else
         {
           thx_phi = -90.0 - thx_degree;
         }  
         
         theta_radian = acos(th);
         theta_degree = theta_radian * 180.0/M_PI;
         
         
         delta_thetaXZ.insert(delta_thetaXZ.begin() + i, theta_degree);
         theta_x.insert(theta_x.begin() + i, thx_phi);
         x[i+1] = x[i] + 2.64;
        
      }
        
     
     //  cout<<"delta_thetaXZ[l]"<<"\t"<<"theta_x"<<endl;
     //  cout<<"----------------------------------------------------"<<endl;
       
     //  for(int l = 0; l < delta_thetaXZ.size(); l++)
      // {
      //   cout<<delta_thetaXZ[l]<<"\t"<<theta_x[l]<<endl;
         
     //  }
      // min = *min_element(delta_thetaXZ.begin(), delta_thetaXZ.end());
       
     //  cout<<"----------------------------------------"<<endl;
     //  cout<<"The minimum value of the delta theta for XZ projection is : "<<min<<endl;
      // cout<<"----------------------------------------"<<endl;
       
       int dim = theta_x.size();
       
       c1->cd();
      // gPad->SetGrid();
       TGraph *gr = new TGraph(dim, &theta_x[0], &delta_thetaXZ[0]);
       gr->SetTitle("Delta theta vs theta for XZ projection");
       gr->GetXaxis()->SetTitle("#theta (degrees)");
       gr->GetYaxis()->SetTitle("#delta#theta (degrees)");
       gr->SetLineColor(kBlue);
       gr->GetYaxis()->SetLimits(-0.5,1.5);
       gr->GetXaxis()->SetLimits(-90.0, 90.0);
       gr->Draw("AC");
       
     //   delta_thetaXZ.insert(delta_thetaXZ.end(), theta_x.begin(), theta_x.end());
        
     return 0;
     
 }
 

 
