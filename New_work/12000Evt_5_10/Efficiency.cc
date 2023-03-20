#include "Libraries2.h"
#include <bits/stdc++.h>

using namespace std;

int main()
{
  gROOT->SetStyle("Plain");
  
  TFile *fileIn = new TFile("output_canvas.root","READ");
  TTree *treeIn = (TTree*)fileIn->Get("TT_Tracks");
  
  Int_t num_trXZ_gen, num_trYZ_gen, num_trXZ_hough, num_trYZ_hough;
  vector<double> * th_rh_XZ_gen = 0; 
  vector<double> * th_rh_YZ_gen = 0; 
  vector<double> * th_rh_XZ_hough = 0; 
  vector<double> * th_rh_YZ_hough = 0;
  
  TBranch *bxz_gen = 0;
  TBranch *byz_gen = 0;
  TBranch *bxz_hough = 0;
  TBranch *byz_hough = 0; 
  
   //treeIn->SetBranchAddress("Evt_ID",&Ev_Id);
   treeIn->SetBranchAddress("Num_TracksXZ_geninfo", &num_trXZ_gen); // contains number of tracks for each event for XZ projection from geninfo
   treeIn->SetBranchAddress("Num_TracksYZ_geninfo",&num_trYZ_gen);  // contains number of tracks for each event for YZ projection from geninfo
   treeIn->SetBranchAddress("ThetaRhoXZ_gen",&th_rh_XZ_gen, &bxz_gen);    // contains theta and rho for each event for XZ projection from geninfo
   treeIn->SetBranchAddress("ThetaRhoYZ_gen",&th_rh_YZ_gen, &byz_gen);    // contains theta and rho for each event for YZ projection from geninfo
   
   treeIn->SetBranchAddress("Num_TracksXZ_hough",&num_trXZ_hough); // contains number of tracks for each event for XZ projection from hough
   treeIn->SetBranchAddress("Num_TracksYZ_hough",&num_trYZ_hough); // contains number of tracks for each event for YZ projection from hough  
   treeIn->SetBranchAddress("ThetaRhoXZ_hough",&th_rh_XZ_hough, &bxz_hough); // contains theta and rho for each event for XZ projection from Hough
   treeIn->SetBranchAddress("ThetaRhoYZ_hough",&th_rh_YZ_hough, &byz_hough); // contains theta and rho for each event for YZ projection from Hough
  /* treeIn->SetBranchAddress("CoordX_hough",&cx_hough);
   treeIn->SetBranchAddress("CoordX",&cx);
   treeIn->SetBranchAddress("CoordZX",&zx);
   treeIn->SetBranchAddress("CoordY_hough",&cy_hough);
   treeIn->SetBranchAddress("CoordY",&cy);
   treeIn->SetBranchAddress("CoordZY",&zy);
   */
   TFile fout("output_eff_plots2.root","recreate");
   
   
   TH1D* h1 = new TH1D("h1", "Resolution XZ; | #Delta (#theta_{hough} - #theta_{gen}) |; Number of single Tracks;", 500, 0.0, 10.0);
   TH1D* h2 = new TH1D("h2", "Resolution YZ; | #Delta (#theta_{hough} - #theta_{gen}) |; Number of single Tracks;", 500, 0.0, 10.0);
      
   int nentries = int(treeIn->GetEntries());
   cout<<"number of entries in TTree = "<<nentries<<endl;
   
   vector<float> d_thXZ, d_thYZ;
   float th_genXZ, th_houghXZ, delta_thXZ, th_genYZ, th_houghYZ, delta_thYZ;
   long  tot_XZ_hough = 0, tot_YZ_hough = 0, tot_XZ_gen = 0, tot_YZ_gen = 0;
   UInt_t j = 0;
  for(int l = 0; l < nentries; l++)
 { 
 
   treeIn -> GetEntry(l); 
   
    Int_t tentry = treeIn->LoadTree(l);
    bxz_gen->GetEntry(tentry);
    byz_gen->GetEntry(tentry);
    bxz_hough->GetEntry(tentry);
    byz_hough->GetEntry(tentry);
    
    
   // cout<<"Ev num = "<<l<<endl;
   // cout<<"----------------------------------"<<endl;
   if(num_trXZ_gen == 1)
   {   th_genXZ = th_rh_XZ_gen -> at(j);
         
      if(num_trXZ_hough == 1)
      {   
         th_houghXZ = th_rh_XZ_hough -> at(j);
         
         tot_XZ_hough += num_trXZ_hough;
       
       //cout<<th_rh_XZ_gen -> at(j)<<"\t"<<th_rh_XZ_hough -> at(j)<<endl; 
           
         delta_thXZ = abs(th_genXZ - th_houghXZ); 
         d_thXZ.push_back(delta_thXZ);   
         h1 ->Fill(delta_thXZ);
      }
    
      tot_XZ_gen += num_trXZ_gen;
          
   }
  
  //cout<<"-------------------------------------------------------------"<<endl;
   if(num_trYZ_gen == 1)
   {  //cout<<"num tr gen = "<<num_trYZ_gen<<endl; 
      if(num_trYZ_hough == 1)
      {  // cout<<"num tr hough = "<<num_trYZ_hough<<endl;
         th_genYZ = th_rh_YZ_gen  -> at(j);
         th_houghYZ = th_rh_YZ_hough  -> at(j);
         tot_YZ_hough += num_trYZ_hough; 
      
       //cout<<th_rh_YZ_gen -> at(j)<<"\t"<<th_rh_YZ_hough -> at(j)<<endl; 
         
         delta_thYZ = abs(th_genYZ - th_houghYZ); 
         d_thYZ.push_back(delta_thYZ);   
         h2 ->Fill(delta_thYZ);
      }  
      tot_YZ_gen += num_trYZ_gen;
          
   }
  
  }
         
      
      cout<<tot_XZ_hough<<"\t"<<tot_XZ_gen<<"\t"<<tot_YZ_hough<<"\t"<<tot_YZ_gen<<endl;
      float effXZ =  ((double)tot_XZ_hough)/tot_XZ_gen;
      float effYZ =  ((double)tot_YZ_hough)/tot_YZ_gen;
      sort(d_thXZ.begin(), d_thXZ.end());
      sort(d_thYZ.begin(), d_thYZ.end());
     /* for(int i = 0; i < d_thXZ.size(); i++)
      {
      cout<<d_thXZ[i]<<endl;
      } */
      cout<<"size of the vector = "<<d_thXZ.size()<<endl;
      int xz1 = d_thXZ.size()/2;
      int yz1 = d_thYZ.size()/2;
      int xz2 = (d_thXZ.size()*9)/10;
      int yz2 = (d_thYZ.size()*9)/10;
      
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"efficiency in XZ Proj: "<< effXZ <<endl;   
      cout<<"efficiency in YZ Proj: "<< effYZ <<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"Resolution in XZ Proj for 50 percent entries: "<< d_thXZ[xz1] <<endl;   
      cout<<"Resolution in YZ Proj for 50 percent entries: "<< d_thYZ[yz1] <<endl;
      cout<<"Resolution in XZ Proj for 90 percent entries: "<< d_thXZ[xz2] <<endl;   
      cout<<"Resolution in YZ Proj for 90 percent entries: "<< d_thYZ[yz2] <<endl;
      cout<<"-------------------------------------------------------------"<<endl;
       
         
  
            fout.cd();
            TCanvas *c = new TCanvas();
            c->SetCanvasSize(1800, 1200);
            c->SetWindowSize(1800, 1000);
            c->Divide(1,2);
           
            c->cd(1);
            gPad->SetGrid();
            h1->SetLineColor(kBlue);
            h1->Draw("hist");
            c->Modified(); c->Update();
            
            c->cd(2);
            gPad->SetGrid();
            h2->SetLineColor(kBlue);
            h2->Draw("hist");
            c->Modified(); c->Update();
            
            c->Write();
            c->Clear();
            
            fout.Close();
  
      return 0;
 }//main 
  
  
  
  
  
  
  
  
  
  
  
  
  
