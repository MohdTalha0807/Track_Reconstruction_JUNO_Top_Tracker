#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstring>

/**** Including the required ROOT Classes  ******/

#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TRootCanvas.h"

//#include "geninfo.h"
#include "TTDigit.h"
#include "CHANNEL_ID_TOOLS.h"


using namespace std;

float Hough_rho(float coord, float z, float theta)
{
  float theta_rad =(M_PI *theta)/(180);
  float c_pos = coord;
  float z_pos = z;
  float rho = c_pos * cos(theta_rad) + z_pos * sin(theta_rad);
  return rho;
}

void remove_duplicate(vector<float> &c_pos)
{
    auto end = c_pos.end();
    for (auto it = c_pos.begin(); it != end; ++it) {
        end = remove(it + 1, end, *it);
    }
 
    c_pos.erase(end, c_pos.end());
}



int main(int argc, char** argv)
{     
   TApplication app("app", &argc, argv);
   
   vector<float> xpos;
   vector<float> zxpos;
   vector<float> ypos;
   vector<float> zypos;
   
  
   TCanvas *c1 = new TCanvas();
   c1->SetCanvasSize(800, 800);
   c1->SetWindowSize(800, 800);
   c1->Divide(1,2);
     
   //TCanvas *c2 = new TCanvas("c2", "c2",72,64,1188,673);
   TCanvas *c3 = new TCanvas();
   c3->SetCanvasSize(800, 800);
   c3->SetWindowSize(800, 800);
   c3->Divide(1,2);
     
   TCanvas *c4 = new TCanvas();
   c4->SetCanvasSize(800, 800);
   c4->SetWindowSize(800, 800);
   c4->Divide(1,2);
  
   TRootCanvas *rc1 = (TRootCanvas *)c1->GetCanvasImp();
   //TRootCanvas *rc2 = (TRootCanvas *)c2->GetCanvasImp();
   TRootCanvas *rc3 = (TRootCanvas *)c3->GetCanvasImp();
   TRootCanvas *rc4 = (TRootCanvas *)c4->GetCanvasImp();
  
   TH2F* hough1 = new TH2F("hough1", "Accumulator; #theta; #rho", 90, -90.0, 90.0, 90, -5.0, 25.0);
   TH2F* hough2 = new TH2F("hough2", "Accumulator; #theta; #rho", 90, -90.0, 90.0, 90, -10.0, 10.0);
   TMultiGraph  *mg1  = new TMultiGraph();
   TMultiGraph  *mg2  = new TMultiGraph();
   
   
   TFile *fileIn1 = new TFile("muon_265_user.root","READ");
   
   TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");
   
   int ch_id[1000], No_Ch, no_hits = 0;
   double charge_pe[1000], threshold_pe = 0.33, theta_max = 90.0;
   float cord_x[1000], cord_y[1000], cord_z[1000]; 
   
   treeIn1->SetBranchAddress("TB_DMchannel", &ch_id);
   treeIn1->SetBranchAddress("NTouchedChannel",&No_Ch);
   treeIn1->SetBranchAddress("TB_pe",&charge_pe); 
   treeIn1->SetBranchAddress("TB_xcC",&cord_x);
   treeIn1->SetBranchAddress("TB_ycC",&cord_y);
   treeIn1->SetBranchAddress("TB_zcC",&cord_z);

   int nentries = int(treeIn1->GetEntries());
    
   
   
   for(int k=0; k<15; k++)
   {
      
      treeIn1->GetEntry(k);
      
      if(No_Ch==0)continue;
      Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0, no_hits_pe = 0;
      Double_t no_pe = 0.0;
      Int_t Wall_idx[3][3][7] = {0}, Wall_idy[3][3][7] = {0};     // 3-d array to count the xy coincidence and store the location of the rejected hits
      
     
      cout<<"Event number : "<<k<<"\t"<<No_Ch<<endl;
      cout<<"-----------------------------------------------------------------------------------"<<endl;
      
      Int_t xcount[3] = {0}, ycount[3] = {0}; 
      
    
          
          for(int hit=0; hit < No_Ch; hit++)
          {
            no_pe = charge_pe[hit];
            if(no_pe >= threshold_pe)
            {
                no_hits_pe ++ ;  // storing the number of hits passed the pe threshold  
            }
          }
          
          
          
          for(int i = 0; i < 3 ; i++)
          {
            for (int j = 0; j < 7 ; j++)
            {
             
             for(int hit=0; hit < No_Ch; hit++)
          {
          
            leaf_channel = ch_id[hit];
            no_pe = charge_pe[hit];
          
            if(no_pe >= threshold_pe)
            {
             leaf_layer = getLayerID(leaf_channel);  
             leaf_row= getRowID(leaf_channel);
             leaf_column=getColumnID(leaf_channel); 
          
              for(int lyr= 0; lyr < 3; lyr ++)
              {
                   if (leaf_layer == lyr)
                  {
                      if(leaf_row == j && leaf_column == i)
                   { leaf_pmt = getPMTID(leaf_channel);   
                 
                      if(leaf_pmt == 0 or leaf_pmt == 1 or leaf_pmt == 2 or leaf_pmt == 3 or leaf_pmt == 8 or leaf_pmt == 9 or leaf_pmt == 10 or leaf_pmt == 11)
                        { xcount[lyr] = xcount[lyr] + 1; Wall_idx[lyr][leaf_column][leaf_row] = Wall_idx[lyr][leaf_column][leaf_row] + 1; }
            
                      if(leaf_pmt == 4 or leaf_pmt == 5 or leaf_pmt == 6 or leaf_pmt == 7 or leaf_pmt == 12 or leaf_pmt == 13 or leaf_pmt == 14 or leaf_pmt == 15)
                        { ycount[lyr] = ycount[lyr] + 1; Wall_idy[lyr][leaf_column][leaf_row] = Wall_idy[lyr][leaf_column][leaf_row] + 1; }
             
                   }
                  }
               } 
             }
           }
          } 
         }
             
 
 /*------------------------- Checking for the XY coincidences on all the three layers and rejecting hits with no such coincidence ----------------------------*/
 
 
 
           
            if(xcount[0] > 0 && xcount[1] > 0 && xcount[2] > 0 && ycount[0] > 0 && ycount[1] > 0 && ycount[2] > 0)
          {   
             
             Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0;
             Double_t no_pe = 0.0; 
             Float_t x_c = 0.0, y_c = 0.0, z_c = 0.0;
             int g = 0, h = 0;
             
              for(int hit=0; hit < No_Ch; hit++)
            {  
               
               leaf_channel = ch_id[hit];
               no_pe = charge_pe[hit];
          
                if(no_pe >= threshold_pe)
                 {
                   leaf_layer = getLayerID(leaf_channel);   // location of the hit.
                   leaf_pmt = getPMTID(leaf_channel);
                   leaf_strip = getStripID(leaf_channel);
                   leaf_row= getRowID(leaf_channel);
                   leaf_column=getColumnID(leaf_channel);
                  
                   x_c = cord_x[hit];
                   y_c = cord_y[hit];
                   z_c = cord_z[hit];
                    
                   if(Wall_idx[leaf_layer][leaf_column][leaf_row] == 0 || Wall_idy[leaf_layer][leaf_column][leaf_row] == 0) 
                     {
                          no_hits = no_hits + 1; 
                          // continue;
                     }
                   else 
                  { 
                  cout<<"layer"<<"\t"<<"row"<<"\t"<<"column"<<"\t"<< "pmt" <<"\t"<< "strip" << "\t"<<"X"<< "\t"<<"Y"<< "\t"<<"Z"<< endl;
                  cout<< leaf_layer<<"\t" << leaf_row<<"\t" << leaf_column<<"\t" << leaf_pmt << "\t" << leaf_strip << "\t" << x_c << "\t" << y_c << "\t" << z_c << endl;
                  cout<<"----------------------------------------------------------------------------------------"<<endl;
                  
                  if(x_c < 1e8)  
                   {  
                      xpos.insert(xpos.begin() + g, x_c);
                      zxpos.insert(zxpos.begin() + g, z_c);
                      g++;
                   }
                   
                  else
                   { ypos.insert(ypos.begin() + h, y_c);
                     zypos.insert(zypos.begin() + h, z_c);
                     h++;
                   } 
                  
                  }
                }
          
              }
           }
   
   
 /*------------------------------------ Applying the Hough Transform ---------------------------------------------*/           
    
        int dim[2] = {0};
        const int r = 90;  // number of steps for the angle theta
        double rho[r], theta[r], rh, th;
        double h = (double) 180/r;
      
        float cord, z;
        theta[0] = -90.0; 
        dim[0] = xpos.size();
        dim[1] = ypos.size();  
        if(dim[0] > 1 && dim[1] > 1)
       {
         
         
       for(int c = 0; c < 2; c++)
       {  
            
         for(int j = 0; j < dim[c]; j ++)
     {      
            if(c == 0)
           {  cord = xpos[j]/1000.0;
              z = zxpos[j]/1000.0; }
           else
           {  cord = ypos[j]/1000.0; 
              z   = zypos[j]/1000.0;  }
                            

      for(int l = 0; l < r; l ++)
     {     
           
         
         rho[l] = Hough_rho(cord,z,theta[l]);
         
         th = theta[l];
         rh = rho[l];      
         
         if(c == 0)
           {  hough1->Fill(th,rh); }
           else
           {  hough2->Fill(th,rh); }
         
         theta[l+1] = theta[l] + h;
      }
       
       if(c == 0)
           {  TGraph* gr3 = new TGraph(r,theta,rho);   
              mg1->Add(gr3); }
           else
           {  TGraph* gr4 = new TGraph(r,theta,rho);   
              mg2->Add(gr4); }
       
    
     }
     
        }
       
 /*------------------------------------------------ Plotting the results ------------------------------------------*/    
       c1->cd(1);
       TGraph* gr1 = new TGraph(dim[0], &xpos[0], &zxpos[0]);
       gr1->SetTitle("XZ Projection");
       gr1->GetXaxis()->SetTitle("X pos of hits");
       gr1->GetYaxis()->SetTitle("Z pos of hits");
       gr1->Draw("AP*");
       c1->Modified(); c1->Update();
       rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       
       c1->cd(2);
       TGraph* gr2 = new TGraph(dim[1], &ypos[0], &zypos[0]);
       gr2->SetTitle("YZ Projection");
       gr2->GetXaxis()->SetTitle("Y pos of hits");
       gr2->GetYaxis()->SetTitle("Z pos of hits");
       gr2->Draw("AP*");           
       c1->Modified(); c1->Update();
       rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       
       c3->cd(1);
       mg1->SetTitle("Parameter Space (XZ)");
       mg1->GetXaxis()->SetTitle("#theta");
       mg1->GetYaxis()->SetTitle("#rho");
       mg1->Draw("AC");
       c3->Modified(); c3->Update();
       rc3->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       
       c3->cd(2);
       mg2->SetTitle("Parameter Space (YZ)");
       mg2->GetXaxis()->SetTitle("#theta");
       mg2->GetYaxis()->SetTitle("#rho");
       mg2->Draw("AC");
       c3->Modified(); c3->Update();
       rc3->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       
       
       c4->cd(1);
       hough1->SetStats(0);
       hough1->DrawCopy("COLZ");
       c4->Modified(); c4->Update();
       rc4->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
       
       c4->cd(2);
       hough2->SetStats(0);
       hough2->DrawCopy("COLZ");
       c4->Modified(); c4->Update();
       rc4->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
        
       app.Run();
      
      
      }
      
   
   
   }

      cout<<"Total number of hits that were rejected : "<<no_hits<<endl;
}









