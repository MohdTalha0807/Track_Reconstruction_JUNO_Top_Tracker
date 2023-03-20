#include "Libraries.h"
#include "TGraphErrors.h"

using namespace std;
#define ncounts 256

//Hough Parametrization
float Hough_rho(float coord, float z, float theta)
{
  float theta_rad =(M_PI *theta)/(180);
  float c_pos = coord;
  float z_pos = z;
  float rho = c_pos * cos(theta_rad) + z_pos * sin(theta_rad);
  return rho;
}

//Hough Accumulator (or Accumulator Matrix)
TH2F* Hough_Acc(int nevt , int nhit_coord , Double_t *coord, Double_t *z , int number_bin_theta, int number_bin_rho,std::string coord_name)
{
  Double_t bin_size_theta = (180.0/(number_bin_theta-1))  ;
  Double_t theta_min = -90 - bin_size_theta/2.0;
  Double_t theta_max = 90 + bin_size_theta/2.0;
  bin_size_theta = ((theta_max-theta_min)/number_bin_theta)  ;
  Double_t bin_size_rho   = (8000/(number_bin_rho-1)) ; //100.0/(number_bin_rho - 1 );
  Double_t rho_min   = -8000 - bin_size_rho/2.0;
  Double_t rho_max   = 8000 + bin_size_rho/2.0;
  bin_size_rho = ((rho_max-rho_min)/number_bin_rho)  ;

  std::string name_acc = "Acc_evt_" + coord_name + std::to_string(nevt) +"_hits_"+ std::to_string(nhit_coord);
  TH2F* Acc = new TH2F(name_acc.c_str(),"Acc", number_bin_theta,theta_min,theta_max ,number_bin_rho,rho_min,rho_max);
  Double_t pho = 0.0;
  Double_t theta= 0.0;

  for(int nhit = 0 ; nhit<nhit_coord ; nhit++)
  {
    for(int theta_deg = 0; theta_deg<number_bin_theta; theta_deg++) // theta in degrees (-90, 90)
    {
      theta = theta_min + (bin_size_theta*(theta_deg+0.5));
      pho = Hough_rho(coord[nhit], z[nhit], theta);
      if(theta>=-80 and theta<=80)
      {
        Acc->Fill(theta , pho);
      }
     }
   }
    return(Acc);
 }

//Main function
int main()
{
   gROOT->SetStyle("Plain");
   vector<float> xpos;
   vector<float> zxpos;
   vector<float> ypos;
   vector<float> zypos;

   TFile *fileIn1 = new TFile("muon_265_user.root","READ");
   TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");
   int ch_id[1000], No_Channels, no_hits = 0;
   float charge_pe[1000];
   float threshold_pe =1.0/3.0;
   double theta_max = 90.0;
   float cord_x[1000], cord_y[1000], cord_z[1000];
   treeIn1->SetBranchAddress("TB_DMchannel", &ch_id);
   treeIn1->SetBranchAddress("NTouchedChannel",&No_Channels);
   treeIn1->SetBranchAddress("TB_pe",&charge_pe);
   treeIn1->SetBranchAddress("TB_xcC",&cord_x);
   treeIn1->SetBranchAddress("TB_ycC",&cord_y);
   treeIn1->SetBranchAddress("TB_zcC",&cord_z);
   int nentries = int(treeIn1->GetEntries());

   TFile f("demo.root","recreate");
   f.cd();
   int nlim2=0;
   cout<<"Please insert, the numbers of events will be analyzed (remember > 0) : ";
   cin>>nlim2;
   if(nlim2>0){cout<<"merci :) "<<endl;
  cout<<"Analyzing .... "<<endl;
  }

   for(int k=0; k<nlim2; k++)
   {
      treeIn1->GetEntry(k);

      if(No_Channels==0)continue;
      Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0, no_hits_pe = 0;
      Double_t no_pe = 0.0;
      Int_t Wall_idx[3][3][7] = {0}, Wall_idy[3][3][7] = {0};     // 3-d array to count the xy coincidence and store the location of the rejected hits

      cout<<"--------------------------------------------------------------------------------- "<<endl;
      cout<<"****  Event number : "<<k<<"\t"<<No_Channels<<endl;
      cout<<"-----------------------------------------------------------------------------------"<<endl;
      Int_t xcount[3] = {0}; //hits in each layer for XZ
      Int_t ycount[3] = {0}; //hits in each layer for YZ

      for(int hit=0; hit< No_Channels; hit++)
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
      int number_layer=3;
      for(int i = 0; i < number_column ; i++)
      {
        for (int j = 0; j < number_row ; j++)
        {
          for(int hit=0; hit < No_Channels; hit++)
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
                     {
                       leaf_pmt = getPMTID(leaf_channel);
                       if(leaf_pmt == 0 or leaf_pmt == 1 or leaf_pmt == 2 or leaf_pmt == 3 or leaf_pmt == 8 or leaf_pmt == 9 or leaf_pmt == 10 or leaf_pmt == 11)
                       {
                         xcount[lyr] = xcount[lyr] + 1;
                         Wall_idx[lyr][leaf_column][leaf_row] = Wall_idx[lyr][leaf_column][leaf_row] + 1;
                       }
                       if(leaf_pmt == 4 or leaf_pmt == 5 or leaf_pmt == 6 or leaf_pmt == 7 or leaf_pmt == 12 or leaf_pmt == 13 or leaf_pmt == 14 or leaf_pmt == 15)
                       {
                         ycount[lyr] = ycount[lyr] + 1; Wall_idy[lyr][leaf_column][leaf_row] = Wall_idy[lyr][leaf_column][leaf_row] + 1;
                       }
                      }
                     }
                   }
                }
              }
            }
           }
/*------------------------- Checking for the XY coincidences on all the three layers and rejecting hits with no such coincidence ----------------------------*/
            Double_t X[100];
            Double_t ZX[100];
            Double_t errorX[100];
            Double_t errorZX[100];
            Double_t Y[100];
            Double_t ZY[100];
            Double_t errorY[100];
            Double_t errorZY[100];
            //Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0;

            Float_t x_c = 0.0, y_c = 0.0, z_c = 0.0;
            int g = 0, h = 0;
           if(xcount[0] > 0 && xcount[1] > 0 && xcount[2] > 0 && ycount[0] > 0 && ycount[1] > 0 && ycount[2] > 0)
           {


             for(int hit=0; hit < No_Channels; hit++)
             {
               leaf_channel = ch_id[hit];

               if(charge_pe[hit] >= threshold_pe)
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
                     no_hits++;
                     continue;
                   }
                   else
                   {
                     cout<<"*layer"<<"\t"<<"*row"<<"\t"<<"*column"<<"\t"<< "*pmt" <<"\t"<< "*strip" << "\t"<<"*X"<< "\t"<<"*Y"<< "\t"<<"*Z"<< endl;
                     cout<< leaf_layer<<"\t" << leaf_row<<"\t" << leaf_column<<"\t" << leaf_pmt << "\t" << leaf_strip << "\t" << x_c << "\t" << y_c << "\t" << z_c << endl;
                     cout<<"----------------------------------------------------------------------------------------"<<endl;

                     if(x_c < 1e8)
                     {
                        xpos.insert(xpos.begin() + g, x_c);
                        zxpos.insert(zxpos.begin() + g, z_c);

                        X[g]=cord_x[hit];
                        ZX[g]=cord_z[hit];
                        errorX[g] = 2.64/(sqrt(12.0));
                        errorZX[g] = 1.06/(sqrt(12.0));
                        g++;
                     }

                     if(y_c < 1e8)
                     {
                       ypos.insert(ypos.begin() + h, y_c);
                       zypos.insert(zypos.begin() + h, z_c);
                       Y[h]=cord_y[hit];
                       ZY[h]=cord_z[hit];
                       errorX[h] = 2.64/(sqrt(12.0));
                       errorZX[h] = 1.06/(sqrt(12.0));
                       h++;
                     }
                   }
                }
              }
            }




/*------------------------------------------------ Applying the Hough Transform --------------------------------------------------*/
        const int nbintheta=61;
        const int nbinrhox=1201;
        const int nbinrhoy=1001;
        const int point_XZ = g;
        const int point_YZ = h;

        TH2F* houghXZ;
        houghXZ = Hough_Acc(k, point_XZ , X, ZX, nbintheta, nbinrhox , "histo_ZX_");


        TH2F* houghYZ;
        houghYZ = Hough_Acc(k, point_YZ , Y, ZY, nbintheta, nbinrhoy , "histo_ZY_");
        f.cd();

         std::string coord_name_canvas =  "Hough_simu"+ std::to_string(k);
         TCanvas *c1 = new TCanvas("c1",coord_name_canvas.c_str());
         c1->SetCanvasSize(1800, 1200);
         c1->SetWindowSize(1800, 1000);
         c1->Divide(2,2);
         c1->SetName(coord_name_canvas.c_str());





/*------------------------------------------------ Plotting the results ------------------------------------------*/
       cout<<"Canvas Started !!!!"<<endl;



       c1->cd(1);
       gPad->SetGrid();
       gROOT->ForceStyle();
       TGraphErrors *gr1;
       gr1 = new TGraphErrors (point_XZ, X, ZX,errorX, errorZX);
       std::string name_xz = "XZ : Evt "+ std::to_string(k);
       gr1->SetName(name_xz.c_str());
       gr1->SetTitle(name_xz.c_str());
       gr1->GetXaxis()->SetTitle("X pos of hits");
       gr1->GetYaxis()->SetTitle("Z pos of hits");
       gr1->SetLineColor(kRed);
       gr1->SetMarkerColor(kRed);
       gr1->SetMarkerSize(1.0);
       gr1->SetMarkerStyle(20);
       gr1->GetYaxis()->SetLimits(-3000,3000);
       gr1->GetXaxis()->SetLimits(-10000, 10000);
       gr1->Draw("AP*");
       c1->Modified(); c1->Update();
       //rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

       c1->cd(3);
       gPad->SetGrid();
       gROOT->ForceStyle();
       TGraphErrors *gr2;
       gr2 = new TGraphErrors (point_YZ, Y, ZY,errorY, errorZY);
       std::string name_yz = "YZ : Evt "+ std::to_string(k);
       gr2->SetName(name_yz.c_str());
       gr2->SetTitle(name_yz.c_str());
       gr2->GetXaxis()->SetTitle("Y pos of hits");
       gr2->GetYaxis()->SetTitle("Z pos of hits");
       gr2->SetLineColor(kBlue);
       gr2->SetMarkerColor(kBlue);
       gr2->SetMarkerSize(1.0);
       gr2->SetMarkerStyle(20);
       gr2->GetYaxis()->SetLimits(-3000,3000);
       gr2->GetXaxis()->SetLimits(-10000, 10000);
       gr2->Draw("AP*");
       c1->Modified(); c1->Update();
       //rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

       c1->cd(2);
       gPad->SetGrid();
       //houghXZ->GetYaxis()->SetLimits(-1000,1000);
       houghXZ->GetXaxis()->SetLimits(-100, 100);
       std::string name_acc_xz = "ACC_XZ : Evt "+ std::to_string(k);
       houghXZ->SetName(name_acc_xz.c_str());
       houghXZ->SetTitle(name_acc_xz.c_str());
       houghXZ->SetStats(0);
       gROOT->ForceStyle();
       gStyle->SetPalette(55);
       houghXZ->Draw("COLZ");
       c1->Modified(); c1->Update();
       //rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

       c1->cd(4);
       gPad->SetGrid();
       //houghYZ->GetYaxis()->SetLimits(-1000,1000);
       std::string name_acc_yz = "ACC_YZ : Evt "+ std::to_string(k);
       houghYZ->SetName(name_acc_yz.c_str());
       houghYZ->SetTitle(name_acc_yz.c_str());
       houghYZ->GetXaxis()->SetLimits(-100, 100);
       houghYZ->SetStats(0);
       gROOT->ForceStyle();
       gStyle->SetPalette(55);
       houghYZ->Draw("COLZ");
       c1->Modified(); c1->Update();
       //rc1->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");


       //app.Run();

       c1->Write();
       c1->Clear();
       int nhits_rejected=0;
       nhits_rejected=no_hits;
       cout<<"Total number of hits that were rejected : "<<nhits_rejected<<endl;

}

   f.Close();
return 0;
}
