#include "binsize_theta_y.h"

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
//#include <vector>
#include <algorithm>
#include <cmath>



vector<float> binsize_theta_y()
{    
     TLatex Tl;
     Tl.SetTextSize(0.03);
    // TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
     
     double strip = 2.64;
     double numerator, denominator;
     const int c = 2000;
     double y[c], z, th, min, theta_radian, theta_degree, thy_radian, thy_degree; 
     std::vector<float> delta_thetaYZ;
     std::vector<float> theta_y; 
 //  std::vector<float> thbinyz_th;
     int cntr =0;
     
          y[0] = -1000.0;
          z = 150;
        
      for(int i =0; i < 758; i++)
      {
         numerator   = pow(y[i],2) + z*z + y[i]*strip;
         denominator = sqrt(pow(y[i],2) + z*z) * sqrt(pow(y[i] + strip,2) + z*z);
      
         th = numerator/denominator;
         
         thy_radian = atan(z/y[i]);
         thy_degree = thy_radian * 180.0/M_PI;
         
         theta_radian = acos(th);
         theta_degree = theta_radian * 180.0/M_PI;
         
         delta_thetaYZ.insert(delta_thetaYZ.begin() + i, theta_degree);
         theta_y.insert(theta_y.begin() + i, thy_degree);
           
         y[i+1] = y[i] + 2.64;
      }
        
  
     delta_thetaYZ.insert(delta_thetaYZ.end(), theta_y.begin(), theta_y.end());
  
  
  //     cout<<"delta_theta"<<"\t"<<"theta_y"<<endl;
  //     cout<<"-------------------------------------------"<<endl;
        
  /*      for(int l = 0; l < delta_thetaYZ.size(); l++)
       {
         cout<<delta_thetaYZ[l]<<endl;
         
       } */
   //    min = *min_element(delta_theta.begin(), delta_theta.end());
       
    //   cout<<"----------------------------------------"<<endl;
    //   cout<<"The minimum value of the delta theta for YZ projection is : "<<min<<endl;
    //   cout<<"----------------------------------------"<<endl;
       
     //  int dim = theta_y.size();
      /* 
       c1->cd();
       gPad->SetGrid();
       TGraph *gr = new TGraph(dim, &theta_y[0], &delta_theta[0]);
       gr->SetTitle("Delta theta vs theta for YZ projection");
       gr->GetXaxis()->SetTitle("#theta (degrees)");
       gr->GetYaxis()->SetTitle("#delta#theta (degrees)");
       //gr->GetYaxis()->SetLimits(0.0,2.64);
       //gr->GetXaxis()->SetLimits(-90.0, 90.0);
       gr->Draw("AC");
       */
     return delta_thetaYZ;
     
 }
 
 
 /*
 
         vector<float> binsize_thetaxz;
        vector<float> binsize_thetayz;
        vector<float> binsize_rhoxz;
        vector<float> binsize_rhoyz;
        
        vector<float> delta_thetaxz;
        vector<float> delta_thetayz;
        vector<float> delta_rhoxz;
        vector<float> delta_rhoyz;

        
        binsize_thetaxz = binsize_theta_x();
        binsize_thetayz = binsize_theta_y();
        binsize_rhoxz   = binsize_rho_x();
        binsize_rhoyz   = binsize_rho_y();
          
        int cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
*/

 
