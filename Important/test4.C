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

void test4()
{
   TFile *fileIn1 = new TFile("muon_265_user.root","READ");
   
   TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");
   
   int ch_id[1000], No_Ch;
   double charge_pe[1000];
   Int_t Wall_id1[21][21], Wall_id2[21][21], Wall_id3[21][21];
   
   for(int i = 0; i < 21 ; i++){
   for (int j = 0; j < 21 ; j++){
        Wall_id1[i][j] = 0;
        Wall_id2[i][j] = 0;
        Wall_id3[i][j] = 0;
    }
    }
   
   treeIn1->SetBranchAddress("TB_DMchannel", &ch_id);
   treeIn1->SetBranchAddress("NTouchedChannel",&No_Ch);
   treeIn1->SetBranchAddress("TB_pe",&charge_pe);


   int nentries = int(treeIn1->GetEntries());
   
   
   for(int k=0; k<1000; k++)
   {
      
      treeIn1->GetEntry(k);
      
      if(No_Ch==0)continue;
      Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0;
      Double_t no_pe = 0.0;
    
      cout<<"Event number : "<<k<<endl;
      cout<<"-----------------------------------------------------------------------------------"<<endl;
      
      int x1count = 0, x2count = 0, x3count = 0, y1count = 0, y2count = 0, y3count = 0; 
 
    
          for(int i = 0; i < 3 ; i++)
          {
            for (int j = 0; j < 7 ; j++)
            {
             
             for(int hit=0; hit < No_Ch; hit++)
          {
          
            leaf_channel = ch_id[hit];
            no_pe = charge_pe[hit];
          
            if(no_pe >= 0.33)
            {
             leaf_layer = getLayerID(leaf_channel);   // location of the hit.
          
          
                if(leaf_layer == 0)
            {
               leaf_row= getRowID(leaf_channel);
               leaf_column=getColumnID(leaf_channel);
             
               if(leaf_row == j && leaf_column == i)
               { leaf_pmt = getPMTID(leaf_channel);   
                 
                if(leaf_pmt == 0 or leaf_pmt == 1 or leaf_pmt == 2 or leaf_pmt == 3 or leaf_pmt == 8 or leaf_pmt == 9 or leaf_pmt == 10 or leaf_pmt == 11)
                { x1count = x1count + 1;}
            
                if(leaf_pmt == 4 or leaf_pmt == 5 or leaf_pmt == 6 or leaf_pmt == 7 or leaf_pmt == 12 or leaf_pmt == 13 or leaf_pmt == 14 or leaf_pmt == 15)
                { y1count = y1count + 1;}
             
               }     
            }
          
          else if(leaf_layer == 1)
          {
             leaf_row= getRowID(leaf_channel);
             leaf_column=getColumnID(leaf_channel);
         
             if(leaf_row == j && leaf_column == i)
             { leaf_pmt = getPMTID(leaf_channel);   
                 
                if(leaf_pmt == 0 or leaf_pmt == 1 or leaf_pmt == 2 or leaf_pmt == 3 or leaf_pmt == 8 or leaf_pmt == 9 or leaf_pmt == 10 or leaf_pmt == 11)
                { x2count = x2count + 1;}
            
                if(leaf_pmt == 4 or leaf_pmt == 5 or leaf_pmt == 6 or leaf_pmt == 7 or leaf_pmt == 12 or leaf_pmt == 13 or leaf_pmt == 14 or leaf_pmt == 15)
                { y2count = y2count + 1; }
             
             }      
          
          }
          
          else
          {
             leaf_row= getRowID(leaf_channel);
             leaf_column=getColumnID(leaf_channel);
              
           if(leaf_row == j && leaf_column == i)
             { leaf_pmt = getPMTID(leaf_channel);   
                 
                if(leaf_pmt == 0 or leaf_pmt == 1 or leaf_pmt == 2 or leaf_pmt == 3 or leaf_pmt == 8 or leaf_pmt == 9 or leaf_pmt == 10 or leaf_pmt == 11)
                 { x3count = x3count + 1; }
            
                if(leaf_pmt == 4 or leaf_pmt == 5 or leaf_pmt == 6 or leaf_pmt == 7 or leaf_pmt == 12 or leaf_pmt == 13 or leaf_pmt == 14 or leaf_pmt == 15)
                 { y3count = y3count + 1; }
             
             }   
            
          }
          }
         } 
             
             
            
          }
          }
   
          if(x1count > 0 && x2count > 0 && x3count > 0 && y1count > 0 && y2count > 0 && y3count > 0)
          {   
             Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0;
             Double_t no_pe = 0.0;
             
              for(int hit=0; hit < No_Ch; hit++)
            {   
               leaf_channel = ch_id[hit];
               no_pe = charge_pe[hit];
          
                if(no_pe >= 0.33)
                 {
                   leaf_layer = getLayerID(leaf_channel);   // location of the hit.
                   leaf_pmt = getPMTID(leaf_channel);
                   leaf_strip = getStripID(leaf_channel);
                   leaf_row= getRowID(leaf_channel);
                   leaf_column=getColumnID(leaf_channel);

          
          
          cout<<"layer"<<"\t"<<"row"<<"\t"<<"column"<<"\t"<<"pmt"<<endl;
          cout<< leaf_layer<<"\t" << leaf_row<<"\t" << leaf_column<<"\t" << leaf_pmt <<endl;
          cout<<"----------------------------------------------------------------------------------------"<<endl;
        
        }
          
         }
            
        }
     
}
}







