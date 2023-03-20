#include "Libraries.h"
#include <bits/stdc++.h>
using namespace std;

//Hough Parametrization

float Hough_rho(float coord, float z, float theta)
{
  float theta_rad =(M_PI *theta)/(180);  // converting from degree to radian
  float c_pos = coord;   // coordinate position
  float z_pos = z;      
  float rho = c_pos * cos(theta_rad) + z_pos * sin(theta_rad);
  return rho;
}

//Hough Accumulator (or Accumulator Matrix)

TH2F* Hough_Acc(int nevt , int nhit_coord , Double_t *coord, Double_t *z , int number_bin_theta, int number_bin_rho, std::string coord_name, int rho)
{
  Double_t bin_size_theta = (180.0/(number_bin_theta-1))  ;
  Double_t theta_min = -90 - bin_size_theta/2.0;
  Double_t theta_max = 90 + bin_size_theta/2.0;
  
  Double_t bin_size_rho   = (rho/(number_bin_rho-1)) ; // if we take rho to be in cm
  Double_t rho_min   = -rho/2.0 - bin_size_rho/2.0;
  Double_t rho_max   =  rho/2.0 + bin_size_rho/2.0;

  std::string name_acc = "Acc_evt_" + coord_name + std::to_string(nevt) +"_hits_"+ std::to_string(nhit_coord);
  TH2F* Acc = new TH2F(name_acc.c_str(),"Acc", number_bin_theta,theta_min,theta_max ,number_bin_rho,rho_min,rho_max);
  Double_t rh = 0.0;
  Double_t theta= 0.0;

  for(int nhit = 0 ; nhit < nhit_coord ; nhit++)
  {
    for(int theta_deg = 0; theta_deg < number_bin_theta; theta_deg++) // theta in degrees (-90, 90)
    {
      theta = theta_min + (bin_size_theta * (theta_deg+0.5)); 
      rh = Hough_rho(coord[nhit], z[nhit], theta);
 
        Acc->Fill(theta , rh);
    
     }
   }
    return(Acc);
 }


// Removing identical points from the hough analysis
vector<float> remove_same_points(int* dim, vector<float> &xpos, vector<float> &zxpos, vector<float> &ypos, vector<float> &zypos)
{     
      vector<float> cxpos;
      vector<float> cypos; 
      
      int npos_x = dim[0];
      int npos_y = dim[1];
      int same_coord_counterXZ[npos_x] = {0};
      int same_coord_counterYZ[npos_y] = {0};

      
          for(int v = 0; v < npos_x; v++)
        {
          for(int b = 0; b < npos_x; b++)
          {
              if(xpos[b]== xpos[v] && zxpos[v]==zxpos[b])
              {
                 same_coord_counterXZ[v] = same_coord_counterXZ[v] + 1;
              
              }
          }
        }
        
           for(int v = 0; v < npos_y; v++ )
        {
          for(int b = 0; b < npos_y; b++)
          {
              if(ypos[v] == ypos[b] && zypos[v]==zypos[b])
              {
                 same_coord_counterYZ[v] = same_coord_counterYZ[v] + 1;
              
              }
          }
        }
      
       int gh = 0, hi = 0, xsize =0, ysize = 0;
       
          int ic = 0, nc = 0, reject = 0;
            
         for(int j = 0; j < dim[0]; j ++)
     {      
            if(same_coord_counterXZ[j] >= 2)
           {   
               ic = ic + 1;
               if(ic > 1)
               { reject++; 
                 ic = 0;
                 continue;
               }
               else
              {  
               cxpos.insert(cxpos.begin() + gh, xpos[j]);
               cxpos.insert(cxpos.begin() + gh + 1, zxpos[j]);
               gh = gh + 2;
              }
            }
              else
              { 
                cxpos.insert(cxpos.begin() + gh, xpos[j]);
                cxpos.insert(cxpos.begin() + gh + 1, zxpos[j]);
                gh = gh + 2;
              }
       }   
       
        for(int j = 0; j < dim[1]; j++)
       {      
              if(same_coord_counterYZ[j] >= 2)
           {   nc = nc + 1;
               if(nc > 1)
               { reject++; 
                 nc = 0;
                 continue;
               }
               else
              {
               cypos.insert(cypos.begin() + hi, ypos[j]);
               cypos.insert(cypos.begin() + hi + 1, zypos[j]);
               hi = hi + 2;  
              }
           }
             else
             {  
                 cypos.insert(cypos.begin() + hi, ypos[j]);
                 cypos.insert(cypos.begin() + hi + 1, zypos[j]);
                 hi = hi + 2;
                
             }
             

         }
         
     xsize = cxpos.size();
     ysize = cypos.size();

     cxpos.insert(cxpos.end(), cypos.begin(), cypos.end());
     cxpos.push_back(xsize);
     cxpos.push_back(ysize);

     return cxpos;
}



//Main function

int main()
{  
   gROOT->SetStyle("Plain"); // what does this mean and why we are using it?

   vector<float> xpos;
   vector<float> zxpos;
   vector<float> ypos;
   vector<float> zypos;

   TFile *fileIn1 = new TFile("muon_265_user.root","READ");
   TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");

   int ch_id[1000], No_Channels, no_hits = 0;
   float charge_pe[1000];
   float threshold_pe =1.0/3.0;
   
   float cord_x[1000], cord_y[1000], cord_z[1000];
   
   treeIn1->SetBranchAddress("TB_DMchannel", &ch_id);
   treeIn1->SetBranchAddress("NTouchedChannel",&No_Channels);
   treeIn1->SetBranchAddress("TB_pe",&charge_pe);
   treeIn1->SetBranchAddress("TB_xcC",&cord_x);
   treeIn1->SetBranchAddress("TB_ycC",&cord_y);
   treeIn1->SetBranchAddress("TB_zcC",&cord_z);
   //int nentries = int(treeIn1->GetEntries());

   TFile f("output_canvas.root","recreate");
   f.cd();
   int lim2=0;
   cout<<"Please insert the number of events to be analyzed (remember > 0) : ";
   cin>>lim2;
   if(lim2>0)
   {    cout<<"merci :) "<<endl;
        cout<<"Analyzing .... "<<endl;
   }

   for(int k=0; k < lim2; k++)
   { 
      treeIn1->GetEntry(k);

      if(No_Channels==0)continue;
      Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0, no_hits_pe = 0;
      //Double_t no_pe = 0.0;
      Int_t Wall_idx[3][3][7] = {0}, Wall_idy[3][3][7] = {0};     // 3-d array to count the xy coincidence and store the location of the rejected hits

      cout<<"--------------------------------------------------------------------------------- "<<endl;
      cout<<"****  Event number : "<<k<<"\t"<<No_Channels<<endl;
      cout<<"-----------------------------------------------------------------------------------"<<endl;
      Int_t xcount[3] = {0}; //hits in each layer for XZ
      Int_t ycount[3] = {0}; //hits in each layer for YZ

      for(int hit=0; hit < No_Channels; hit++)
      {
        if(charge_pe[hit] >= threshold_pe)
        {
          no_hits_pe++ ;  // storing the number of hits passed the pe threshold
        }
      }
      cout<<"Number of hits with larger charge than 1/3 p.e : "<<no_hits_pe<<endl;
      cout<<"-------------------------------------------------------------"<<endl;

      int number_column=3;
      int number_row=7;
      //int number_layer=3;

      for(int i = 0; i < number_column ; i++)
      { 
        for (int j = 0; j < number_row ; j++)
        { 
          for(int hit=0; hit < No_Channels; hit++)
          { 
            leaf_channel = ch_id[hit];
            

            if(charge_pe[hit] >= threshold_pe)
            { 
               leaf_layer  = getLayerID(leaf_channel);
               leaf_row    = getRowID(leaf_channel);
               leaf_column = getColumnID(leaf_channel);

               for(int lyr= 0; lyr < 3; lyr ++)
               { 
                   if (leaf_layer == lyr)
                   {
                     if(leaf_row == j && leaf_column == i)
                     {
                       leaf_pmt = getPMTID(leaf_channel);
                       if(leaf_pmt == 0 or leaf_pmt == 1 or leaf_pmt == 2 or leaf_pmt == 3 or leaf_pmt == 8 or leaf_pmt == 9 or leaf_pmt == 10 or leaf_pmt == 11)
                       {
                         xcount[lyr] = xcount[lyr] + 1;
                         Wall_idx[lyr][leaf_column][leaf_row] = Wall_idx[lyr][leaf_column][leaf_row] + 1;
                       }
                       else
                       {
                         ycount[lyr] = ycount[lyr] + 1; 
                         Wall_idy[lyr][leaf_column][leaf_row] = Wall_idy[lyr][leaf_column][leaf_row] + 1;
                       }
                      }
                     }
                   }
                }
              }
            }
           }
                      
/*------------------------- Checking for the XY coincidences on all the three layers and rejecting hits with no such coincidence ----------------------------*/
            
           
        

            float x_c = 0.0, y_c = 0.0, z_c = 0.0;
            int g = 0, h = 0;
           if(xcount[0] > 0 && xcount[1] > 0 && xcount[2] > 0 && ycount[0] > 0 && ycount[1] > 0 && ycount[2] > 0)
           { 


             for(int hit=0; hit < No_Channels; hit++)
             { 
               leaf_channel = ch_id[hit];

               if(charge_pe[hit] >= threshold_pe)
               { 
                   leaf_layer = getLayerID(leaf_channel);   // location of the hit.
                   leaf_row= getRowID(leaf_channel);
                   leaf_column=getColumnID(leaf_channel);

                   x_c = cord_x[hit]/10;
                   y_c = cord_y[hit]/10;
                   z_c = cord_z[hit]/10;

                   if(Wall_idx[leaf_layer][leaf_column][leaf_row] == 0 || Wall_idy[leaf_layer][leaf_column][leaf_row] == 0)
                   {
                     no_hits++;
                     continue;
                   }
                   else
                   { 
                     leaf_pmt = getPMTID(leaf_channel);
                     leaf_strip = getStripID(leaf_channel);

                     cout<< "layer" << "\t" << "row" << "\t" << "column" << "\t" << "pmt" << "\t" << "strip" << "\t" << "X" << "\t" << "Y" << "\t" << "Z" << "\t" << "charge" << endl;
                     cout<< leaf_layer<< "\t" << leaf_row <<"\t" << leaf_column<<"\t" << leaf_pmt << "\t" << leaf_strip << "\t" << x_c << "\t" << y_c << "\t" << z_c << "\t" << charge_pe[hit] << endl;
                     cout<<"----------------------------------------------------------------------------------------"<<endl;

                     if(x_c < 1e8)
                     {
                        xpos.insert(xpos.begin() + g, x_c);
                        z_c = z_c + 100.0;
                        zxpos.insert(zxpos.begin() + g, z_c);
                        g++;
                     }
                      else
                     {
                       ypos.insert(ypos.begin() + h, y_c);
                       z_c = z_c + 100.0;
                       zypos.insert(zypos.begin() + h, z_c);
                       h++;
                     }
                   }
                }
              }
            }




/*------------------------------------------------ Applying the Hough Transform --------------------------------------------------*/

        const int nbintheta = 901;    // binsize theta = 0.2 degree
        //const int nbinthetay = 1801;
        const int nbinrhox = 9001;    // binsize X = 0.50 cm
        const int nbinrhoy = 4001;     // binsize Y = 0.50 cm
        const int rho_x = 4500;
        const int rho_y = 2000;
        
         Double_t X[100];
         Double_t ZX[100];
         Double_t Y[100];
         Double_t ZY[100];
        
        int dim[2] = {0}; 
        dim[0] = xpos.size();
        dim[1] = ypos.size();
        vector<float> cpos;
        int x_size = 0, y_size = 0;
        int jk = 0, ji = 0; // ii = 0;
   
           cpos = remove_same_points(dim, xpos, zxpos, ypos, zypos);    // removing identical points from the hough analysis

        y_size = cpos.back();
        cpos.pop_back();
        x_size = cpos.back();
        cpos.pop_back();
     /*   
       cout<<"-----------------------------------------------------------------"<<endl;
       cout << "X coord" <<"\t" << "\t" << "Z coord" << endl;
       cout<<"--------------------------"<<endl;  
         for(auto i = 0; i < x_size; i=i+2)
         {  
            cout << cpos[i] <<"\t" << "\t" << cpos[i+1] << endl;            
         }
       cout<<"--------------------------"<<endl;
           ii = x_size; 
       cout << "Y coord" <<"\t" << "\t" << "Z coord" << endl; 
       cout<<"--------------------------"<<endl;  
        for(auto i = 0; i < y_size; i=i+2)
         {  
           cout << cpos[ii] <<"\t" << "\t" << cpos[ii+1] << endl;            
           ii = ii + 2;
         }
       
      */  
        for(int i = 0; i < x_size; i = i + 2)
        {
            X[jk] = cpos[i];
           ZX[jk] = cpos[i+1];
             jk++;
        }
        
          ji = x_size;   
        for(int i = 0; i < y_size; i++)
        {
            Y[i] = cpos[ji];
           ZY[i] = cpos[ji + 1];
             ji = ji+2;
        }
        
        int point_XZ = x_size/2;
        int point_YZ = y_size/2;
        
        TH2F* houghXZ;
        houghXZ = Hough_Acc(k, point_XZ , X, ZX, nbintheta, nbinrhox , "histo_ZX_", rho_x);

        TH2F* houghYZ;
        houghYZ = Hough_Acc(k, point_YZ , Y, ZY, nbintheta, nbinrhoy , "histo_ZY_", rho_y);
        f.cd();

         std::string coord_name_canvas =  "Hough_simu"+ std::to_string(k);
         TCanvas *c1 = new TCanvas("c1",coord_name_canvas.c_str());
         c1->SetCanvasSize(1800, 1200);
         c1->SetWindowSize(1800, 1000);
         c1->Divide(2,2);
         c1->SetName(coord_name_canvas.c_str());


/*----------------------------------------------------------- Plotting the results -----------------------------------------------------------*/

       cout<<"Canvas Started !!!!"<<endl;

       c1->cd(1);
       gPad->SetGrid();
       gROOT->ForceStyle(); // What is this used for?
       TGraph *gr1;
       gr1 = new TGraph(dim[0], &xpos[0], &zxpos[0]);
       std::string name_xz = "XZ : Evt "+ std::to_string(k);
       gr1->SetName(name_xz.c_str());
       gr1->SetTitle(name_xz.c_str());
       gr1->GetXaxis()->SetTitle("X pos of hits");
       gr1->GetYaxis()->SetTitle("Z pos of hits");
       gr1->SetLineColor(kRed);
       gr1->SetMarkerColor(kRed);
       gr1->SetMarkerSize(1.0);
       gr1->SetMarkerStyle(20);
       gr1->GetYaxis()->SetLimits(-500,500);
       gr1->GetXaxis()->SetLimits(-2000, 2000);
       gr1->Draw("AP*");
       c1->Modified(); c1->Update();
      

       c1->cd(3);
       gPad->SetGrid();
       gROOT->ForceStyle();
       TGraph *gr2;
       gr2 = new TGraph(dim[1], &ypos[0], &zypos[0]);
       std::string name_yz = "YZ : Evt "+ std::to_string(k);
       gr2->SetName(name_yz.c_str());
       gr2->SetTitle(name_yz.c_str());
       gr2->GetXaxis()->SetTitle("Y pos of hits");
       gr2->GetYaxis()->SetTitle("Z pos of hits");
       gr2->SetLineColor(kBlue);
       gr2->SetMarkerColor(kBlue);
       gr2->SetMarkerSize(1.0);
       gr2->SetMarkerStyle(20);
       gr2->GetYaxis()->SetLimits(-500,500);
       gr2->GetXaxis()->SetLimits(-1000, 1000);
       gr2->Draw("AP*");
       c1->Modified(); c1->Update();
       

       c1->cd(2);
       gPad->SetGrid();
       
       houghXZ->GetXaxis()->SetLimits(-90, 90);
       std::string name_acc_xz = "ACC_XZ : Evt "+ std::to_string(k);
       houghXZ->SetName(name_acc_xz.c_str());
       houghXZ->SetTitle(name_acc_xz.c_str());
       houghXZ->SetStats(0);
       gROOT->ForceStyle();
       gStyle->SetPalette(55);
       houghXZ->Draw("COLZ");
       c1->Modified(); c1->Update();
      

       c1->cd(4);
       gPad->SetGrid();
      
       std::string name_acc_yz = "ACC_YZ : Evt "+ std::to_string(k);
       houghYZ->SetName(name_acc_yz.c_str());
       houghYZ->SetTitle(name_acc_yz.c_str());
       houghYZ->GetXaxis()->SetLimits(-90, 90);
       houghYZ->SetStats(0);
       gROOT->ForceStyle();
       gStyle->SetPalette(55);
       houghYZ->Draw("COLZ");
       c1->Modified(); c1->Update();
      
       c1->Write();
       c1->Clear();
       int nhits_rejected = 0;
       nhits_rejected = no_hits;
       cout<<"Total number of hits that were rejected : "<<nhits_rejected<<endl;

     xpos.clear();
     zxpos.clear();
     ypos.clear();
     zypos.clear();

}
   f.Close();
   return 0;
}
