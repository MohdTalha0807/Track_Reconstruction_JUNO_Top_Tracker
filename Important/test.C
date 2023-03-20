#include <vector>


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

#include "geninfo.h"
#include "TTDigit.h"
#include "CHANNEL_ID_TOOLS.h"

void test()
{
   std::vector<float> xpos;
   std::vector<float> zxpos;
   std::vector<float> ypos;
   std::vector<float> zypos;
   
   
   TCanvas *c1 = new TCanvas("c1", "c1",72,64,1188,673);
   TCanvas *c2 = new TCanvas("c2", "c2",72,64,1188,673);
   
   TFile *fileIn1 = new TFile("muon_265_user.root","READ");
   
   TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");
   
   int ch_id[1000], No_Ch, no_hits = 0;
   double charge_pe[1000], threshold_pe = 0.33;
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
      //int g = 0;
      //int h = 0;
      //float xpos[6] = {0.0}, ypos[7] = {0.0}, zxpos[6] = {0.0}, zypos[7] = {0.0};
     
      cout<<"Event number : "<<k<<endl;
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
   
           
            if(xcount[0] > 0 && xcount[1] > 0 && xcount[2] > 0 && ycount[0] > 0 && ycount[1] > 0 && ycount[2] > 0)
          {   
             
             Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0;
             Double_t no_pe = 0.0; 
             Float_t x_c = 0.0, y_c = 0.0, z_c = 0.0;
             
             
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
                   { xpos.insert(xpos.begin() + i, x_c);
                     zxpos.insert(zxpos.begin() + hit, z_c);
                   }
                   
                  else
                   { ypos.insert(ypos.begin() + hit, y_c);
                     zypos.insert(zypos.begin() + hit, z_c);
                     
                   } 
                  
                  }
                }
          
              }
           }
           
       
       
       
       
       
       c1->cd();
       TGraph* gr1 = new TGraph(dim1, &xpos[0], &zxpos[0]);
       gr1->SetTitle("XZ Projection");
       gr1->GetXaxis()->SetTitle("X pos of hits");
       gr1->GetYaxis()->SetTitle("Z pos of hits");
       gr1->GetYaxis()->SetLimits(-10000,10000);
       gr1->GetXaxis()->SetLimits(-15000, 15000);
       gr1->Draw("AP*");
       
       c2->cd();
       TGraph* gr2 = new TGraph(dim2, &ypos, &zypos);
       gr2->SetTitle("YZ Projection");
       gr2->GetXaxis()->SetTitle("Y pos of hits");
       gr2->GetYaxis()->SetTitle("Z pos of hits");
       gr2->GetYaxis()->SetLimits(-10000,10000);
       gr2->GetXaxis()->SetLimits(-15000, 15000);
       gr2->Draw("AP*");           
      
       for(int i = 0; i < no_hits_pe; i++)
       {
           cout<<xpos[i]<<"\t"<<zxpos[i]<<"\t"<<ypos[i]<<"\t"<<zypos[i]<<endl;
       }
   
   
   
   }

      cout<<"Total number of hits that were rejected : "<<no_hits<<endl;
}









