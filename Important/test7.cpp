#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstring>

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TRootCanvas.h"

using namespace std;

int main(int argc, char** argv) 
{  
   TApplication app("app", &argc, argv);
   
   
   
   float x = -15.0, y = -15.0, zx = 0.0, zy = 0.0, m = 1.0, c = 2.0;
   int dim;
   TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
   TCanvas *c2 = new TCanvas("c2", "c2",72,64,1188,673);
   TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
   TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();
   
  
   std::vector<float> xpos;
   std::vector<float> zxpos;
   std::vector<float> ypos;
   std::vector<float> zypos;
   int g = 0, h = 0;
   for(int i = -15; i < 15; i++)
   {
       x++;
       y++;
       zx = m*x + c; 
       zy = m*x + c;
       xpos.insert(xpos.begin() + g, x);
       zxpos.insert(zxpos.begin() + g, zx);
       ypos.insert(ypos.begin() + h, y);
       zypos.insert(zypos.begin() + h, zy);
       
   }
    
    for(int i = 0; i < 15; i++)
    {
    
     cout<<xpos[i]<<"\t"<<ypos[i]<<"\t"<<zxpos[i]<<"\t"<<zypos[i]<<endl;   
    
    }  
                          
       dim = xpos.size();
                             
       c1->cd();                
       TGraph* gr1 = new TGraph(dim, &xpos[0], &zxpos[0]);
       gr1->SetTitle("XZ Projection");
       gr1->GetXaxis()->SetTitle("X pos of hits");
       gr1->GetYaxis()->SetTitle("Z pos of hits");
       
       gr1->Draw("AP*");
       c1->Modified(); c1->Update();
       rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       
       
       
       c2->cd();
       TGraph* gr2 = new TGraph(dim, &ypos[0], &zypos[0]);
       gr2->SetTitle("YZ Projection");
       gr2->GetXaxis()->SetTitle("Y pos of hits");
       gr2->GetYaxis()->SetTitle("Z pos of hits");
      
       gr2->Draw("AP*");           
       c2->Modified(); c2->Update();
       rc2->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       app.Run();
       
  return 0;
}







