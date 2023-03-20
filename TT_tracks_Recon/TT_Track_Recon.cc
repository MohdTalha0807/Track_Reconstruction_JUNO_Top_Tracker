/* Further improvements which can be implemented:

- main function could be reduced to just minimum function calls and other primary tasks.
- plotting could be moved to a separate visualization framework  
- Better naming schemes could be implemented

*/
#include "Libraries.h"
#include <bits/stdc++.h>
#include "other_lib.h"
using namespace std;     // absolutely not a good idea to use this here 
                         // better to use scope resolution operator and define namespace individually to avoid conflict

//-----------------------------------------------------------------Main function-------------------------------------------------------------------------//

int main()
{  //1
   gROOT->SetStyle("Plain");

   vector<float> xpos;
   vector<float> zxpos;
   vector<float> ypos;
   vector<float> zypos;
   
   vector<float> Th_RhXZ_gen;
   vector<float> Th_RhYZ_gen;
   vector<float> Th_RhXZ_hough;
   vector<float> Th_RhYZ_hough;
   vector<float> Co_X_hough, Co_X, Co_ZX, Co_Y_hough, Co_Y, Co_ZY;
   int ntracksXZ_gen, ntracksYZ_gen, ntracksXZ_hough, ntracksYZ_hough;
   
    map<int, vector<float>> mymap; 
    vector<float> v;
    vector<float> tmp; 
    int tdm_ch;
    float cx, cy, cz;
    
   fstream myFile;
   myFile.open("ChannelPosition.txt", ios::in);
   
   if(myFile.is_open())
   {
     while( myFile >> tdm_ch >> cx >> cy >> cz)
     {
           v.push_back(cx);         
           v.push_back(cy);
           v.push_back(cz);
           mymap[tdm_ch] = v;
           v.clear(); 
     }
     myFile.close();
   }

   TFile *fileIn1 = new TFile("muon_265_user.root","READ");
   TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");
   TTree *treeIn2 = (TTree*)fileIn1->Get("geninfo");
    


   int ch_id[1000], No_Channels, no_hits = 0, ev_id[1000], nsim_part, gen_sizeofarr, n_tracks1, n_tracks2,  n_tracks = 0, n_tr_gen = 0; 
   float efficiency = 0.0;
   float charge_pe[1000];
   float threshold_pe =1.0/3.0;
   Double_t bc = 0.0;
   float co_x[1000], co_y[1000], co_z[1000], px[1000], py[1000], pz[1000];  // cord_x[1000], cord_y[1000], cord_z[1000]
   float WP_HEIGHT = 4350.0;
   float TT_HEIGHT =  840.0;
   float Z_GLOBAL_SHIFT = (WP_HEIGHT/2.0) + (TT_HEIGHT/2.0) - 100.0;
   int Ev_Id;
   
   treeIn1->SetBranchAddress("TB_DMchannel", &ch_id);
   treeIn1->SetBranchAddress("NTouchedChannel",&No_Channels);
   treeIn1->SetBranchAddress("TB_pe",&charge_pe);

   treeIn2->SetBranchAddress("evtID",&ev_id);
   treeIn2->SetBranchAddress("nInitParticles",&nsim_part);
   treeIn2->SetBranchAddress("InitX",&co_x);
   treeIn2->SetBranchAddress("InitY",&co_y);
   treeIn2->SetBranchAddress("InitZ",&co_z);
   treeIn2->SetBranchAddress("InitPX",&px);
   treeIn2->SetBranchAddress("InitPY",&py);
   treeIn2->SetBranchAddress("InitPZ",&pz);

 
   //int nentries = int(treeIn1->GetEntries());

   TFile f("output_canvas.root","recreate");
   
   TTree *treeOut = new TTree("TT_Tracks","TT_Tracks");
   
   treeOut -> Branch("Evt_ID",&Ev_Id, "Ev_Id/i");
   treeOut -> Branch("Num_TracksXZ_geninfo", &ntracksXZ_gen, "ntracksXZ_gen/I");
   treeOut -> Branch("Num_TracksYZ_geninfo", &ntracksYZ_gen, "ntracksYZ_gen/I");
   treeOut -> Branch("ThetaRhoXZ_gen", &Th_RhXZ_gen);
   treeOut -> Branch("ThetaRhoYZ_gen", &Th_RhYZ_gen); 
  
  
   
   treeOut -> Branch("Num_TracksXZ_hough", &ntracksXZ_hough, "ntracksXZ_hough/I");
   treeOut -> Branch("Num_TracksYZ_hough", &ntracksYZ_hough, "ntracksYZ_hough/I");
   treeOut -> Branch("ThetaRhoXZ_hough", &Th_RhXZ_hough); 
   treeOut -> Branch("ThetaRhoYZ_hough", &Th_RhYZ_hough); 
   treeOut -> Branch("CoordX_hough", &Co_X_hough);
   treeOut -> Branch("CoordX", &Co_X);
   treeOut -> Branch("CoordZX", &Co_ZX);
   treeOut -> Branch("CoordY_hough", &Co_Y_hough);
   treeOut -> Branch("CoordY", &Co_Y);
   treeOut -> Branch("CoordZY", &Co_ZY);
   
   f.cd();
   int lim1=0;                                      //-----------------------------------------Enter Limits---------------------------//
   int lim2=0;
   cout<<"Please insert the number of events to be analyzed (remember > 0) : ";
   //cin>>lim1;
   cin>>lim2;
   if(lim2>0)
   {    cout<<"merci :) "<<endl;
        cout<<"Analyzing .... "<<endl;
   }

   for(int k=lim1; k < lim2; k++)
   { //2
      
      Co_X_hough.clear(); Co_X.clear(); Co_ZX.clear(); Co_Y_hough.clear(); Co_Y.clear(); Co_ZY.clear(); Th_RhXZ_gen.clear(); Th_RhYZ_gen.clear(); Th_RhXZ_hough.clear(); Th_RhYZ_hough.clear();

      treeIn1->GetEntry(k);
      treeIn2->GetEntry(k);
      
      
       float x_gen, y_gen, slope_x, slope_y, intercept_x, intercept_y;
       vector<float> genc_x;
       vector<float> genc_y;
       vector<float> genc_z;
       vector<float> genp_x;
       vector<float> genp_y;
       vector<float> genp_z; 
       
     //cout<<"-------------------------------------------------geninfo data--------------------------------------------------------"<<endl; 

       for(int in = 0; in < nsim_part; in++)
      {  
            slope_x = pz[in]/px[in];
            slope_y = pz[in]/py[in];
            intercept_x = (co_z[in]/10.0 - Z_GLOBAL_SHIFT) - (slope_x * co_x[in]/10.0);
            intercept_y = (co_z[in]/10.0 - Z_GLOBAL_SHIFT) - (slope_y * co_y[in]/10.0); 
            x_gen = (150.0/slope_x) -1.0*intercept_x/slope_x;
            y_gen = (150.0/slope_y) -1.0*intercept_y/slope_y;
        
        if((-2350.0 < x_gen && x_gen < 2350.0) && (-1000.0 < y_gen && y_gen < 1000.0) )
       {
        genc_x.push_back(co_x[in]/10.0);
        genc_y.push_back(co_y[in]/10.0);
        float z_genc = (co_z[in]/10.0) - Z_GLOBAL_SHIFT;
        genc_z.push_back(z_genc);
        genp_x.push_back(px[in]/10.0);
        genp_y.push_back(py[in]/10.0);
        genp_z.push_back(pz[in]/10.0); 
        }
        else{ continue;}
       }
       
       gen_sizeofarr = genc_x.size();

      if(No_Channels==0)continue;
      Int_t leaf_channel=0, leaf_layer=0, leaf_pmt=0, leaf_strip=0 ,leaf_row=0, leaf_column=0, no_hits_pe = 0;
      //Double_t no_pe = 0.0;
      Int_t Wall_idx[3][3][7] = {0}, Wall_idy[3][3][7] = {0};     // 3-d array to count the xy coincidence and store the location of the rejected hits

      cout<<"--------------------------------------------------------------------------------- "<<endl;
      cout<<"****  Event number : "<<k<<"\t"<<No_Channels<<endl;
      cout<<"-----------------------------------------------------------------------------------"<<endl;
      Int_t xcount[3] = {0}; //hits in each layer for XZ
      Int_t ycount[3] = {0}; //hits in each layer for YZ

      int number_column=3;
      int number_row=7;

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
                      
/*--------------------------------------------- Checking for the XY coincidences on all the three layers and rejecting hits with no such coincidence ----------------------------------*/
            
           
        

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
                   
                   tmp = mymap[leaf_channel];

                   x_c = tmp[0]/10.0;
                   y_c = tmp[1]/10.0;
                   z_c = tmp[2]/10.0;
                   
                   tmp.clear();

                   if(Wall_idx[leaf_layer][leaf_column][leaf_row] == 0 || Wall_idy[leaf_layer][leaf_column][leaf_row] == 0)
                   {
                     no_hits++;
                     continue;
                   }
                   else
                   { 
                     leaf_pmt = getPMTID(leaf_channel);
                     leaf_strip = getStripID(leaf_channel);

                    //cout<< "layer" << "\t" << "row" << "\t" << "column" << "\t" << "pmt" << "\t" << "strip" << "\t" << "X" << "\t" << "Y" << "\t" << "Z" << "\t" << "charge" << endl;
                 //cout<< leaf_layer<< "\t" << leaf_row <<"\t" << leaf_column<<"\t" << leaf_pmt << "\t" << leaf_strip << "\t" << x_c << "\t" << y_c << "\t" << z_c << "\t" << charge_pe[hit] << endl;
                     //cout<<"-----------------------------------------------------------------------------------------"<<endl;

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

          
          
        if(xpos.empty() || ypos.empty())
        { continue;     }  

        int nbintheta = 361;    // binsize theta = 0.5 degree
                                      
        int nbinrhox = 1801;    // binsize X = 2.5 cm initial
        int nbinrhoy = 801;     // binsize Y = 2.5 cm initial
        int rho_x = 4500;
        int rho_y = 2000;
        
         Double_t X[100];
         Double_t ZX[100];
         Double_t Y[100];
         Double_t ZY[100];
        
        int dim[2] = {0}; 
        dim[0] = xpos.size();
        dim[1] = ypos.size();
        vector<float> cposxz;
        vector<float> dummy_xz;
        vector<float> dummy_yz;
        vector<float> cposyz;
        int x_size = 0, y_size = 0;
        int jk = 0, ji = 0; // ii = 0;
           
           dummy_xz = xpos;
           dummy_yz = ypos;  
           cposxz = remove_same_points(xpos, zxpos, dummy_xz);    // removing identical points for the hough analysis
           cposyz = remove_same_points(ypos, zypos, dummy_yz);
           
           
           
        y_size = cposyz.size()/3;
        x_size = cposxz.size()/3;
     /*   
       cout<<"-----------------------------------------------------------------"<<endl;
       cout << "X coord" <<"\t" << "\t" << "Z coord" << endl;
       cout<<"--------------------------"<<endl;  
         for(auto i = 0; i < y_size; i++)
         {  
            cout << cposyz[i] <<"\t"<< "\t" << cposyz[y_size + i]<<"\t"<< "\t"<< cposyz[2 * y_size + i]<< endl;            
         }
       cout<<"--------------------------"<<endl;
       cout << "Y coord" <<"\t" << "\t" << "Z coord" << endl; 
       cout<<"--------------------------"<<endl;  
        for(auto i = 0; i < dim[1]; i++)
         {  
           cout << ypos[i] <<"\t" << "\t" << zypos[i] << endl;            
      
         }
         
      */
        for(int i = 0; i < x_size; i++)
        {
            X[i] = cposxz[i];
           ZX[i] = cposxz[x_size + i];
         }
         
        for(int i = 0; i < y_size; i++)
        {
            Y[i] = cposyz[i];
           ZY[i] = cposyz[y_size + i];
        }
        
        int point_XZ = x_size;
        int point_YZ = y_size;
        Double_t th_thxz = 86.0;
        Double_t th_thyz = 82.0;
        Double_t theta_range = 180.0;
        
          
        vector<float> theta_pos_x;
        vector<float> rho_pos_x;
        vector<float> theta_pos_y;
        vector<float> rho_pos_y;
        vector<float> thetaXZ;
        vector<float> rhoXZ;
        vector<float> thetaYZ;
        vector<float> rhoYZ;
           
        
        int npeak_maxXZ = 0, npeak_maxYZ = 0;
        Double_t binc = 0.5;
        Double_t bin_size_theta = ((double)theta_range/(nbintheta-1));
        Double_t bin_size_rhoXZ   = ((double)rho_x/(nbinrhox-1));
        Double_t bin_size_rhoYZ   = ((double)rho_y/(nbinrhoy-1));
        Double_t theta_min = -90.0 - bin_size_theta/2.0;
        Double_t theta_max = 90.0 + bin_size_theta/2.0; 
         
        //cout<<bin_size_rhoXZ<<"\t"<<bin_size_rhoYZ<<endl;
    
           for(int it =0; it < 5; it++)                                       // varying the bin range ----------------------------------------------
         { //3
         
          Double_t rho_minXZ   = -rho_x/2.0 - bin_size_rhoXZ/2.0 + bin_size_rhoXZ/5.0 * it;
          Double_t rho_maxXZ   =  rho_x/2.0 + bin_size_rhoXZ/2.0 + bin_size_rhoXZ/5.0 * it;
          Double_t rho_minYZ   = -rho_y/2.0 - bin_size_rhoYZ/2.0 + bin_size_rhoYZ/5.0 * it;
          Double_t rho_maxYZ   =  rho_y/2.0 + bin_size_rhoYZ/2.0 + bin_size_rhoYZ/5.0 * it;
        
        TH2D* houghXZ;
        TH2D* houghYZ;
        
        std::string histo_zx = "Histo_ZX_" + std::to_string(it);
        std::string histo_zy = "Histo_ZY_" + std::to_string(it);

        
        houghXZ = Hough_Acc(k, point_XZ , X, ZX, nbintheta, nbinrhox , histo_zx.c_str(), rho_x, th_thxz, theta_range, theta_min, theta_max, rho_minXZ, rho_maxXZ, bin_size_theta, bin_size_rhoXZ, binc);
        
        houghYZ = Hough_Acc(k, point_YZ , Y, ZY, nbintheta, nbinrhoy , histo_zy.c_str(), rho_y, th_thyz, theta_range, theta_min, theta_max, rho_minYZ, rho_maxYZ, bin_size_theta, bin_size_rhoYZ, binc);
        f.cd();   
          
        
 /*-------------------------------------------------------------------------------algorithm for isolating the peaks-----------------------------------------------------*/
 
        
        
        int cntxz = 0, cntyz = 0, peak_counterXZ = 0, peak_counterYZ = 0;
        float thix, rhojx, thiy, rhojy, th1xz = 90.0, th2xz = 0.0, th1yz = 90.0, th2yz = 0.0, e = 5.0, e_zoom = 4.0, erhy = 80.0, rh1yz = 1200.0, rh2yz = 0.0, erhx = 80.0, rh1xz = 3000.0, rh2xz = 0.0;  
          
          for(int i= 1; i<= (nbintheta-1); i++)
        {    
             for(int j = 1; j<=(nbinrhox - 1); j++)
           {      
             float bin_contXZ = houghXZ->GetBinContent(i,j);        
                      
             if(bin_contXZ >= 3)                 // checking for the peaks in XZ Hough space ---------------------
            {  
               thix = -90.0 + (i-1)* theta_range/(double)(nbintheta-1);
               rhojx = -2250.0 + (j-1)* rho_x/(double)(nbinrhox-1);
               th2xz = thix;
               rh2xz = rhojx;        
              if((th1xz - e) < th2xz && th2xz < (th1xz + e) && (rh1xz - erhx) < rh2xz && rh2xz < (rh1xz + erhx))
              {
                if(bin_contXZ >= 4)
                 {  
                      peak_counterXZ++;                                           
                      theta_pos_x.insert(theta_pos_x.begin()+ cntxz, thix);
                      rho_pos_x.insert(rho_pos_x.begin()+ cntxz, rhojx);
                      cntxz++;
                 }
                continue;
              }
             else
              {                        
                peak_counterXZ++;                                           
                theta_pos_x.insert(theta_pos_x.begin()+ cntxz, thix);
                rho_pos_x.insert(rho_pos_x.begin()+ cntxz, rhojx);
                cntxz++;
                th1xz = th2xz;
                rh1xz = rh2xz;
                continue;
              }
            }
          }
        }  
        
   /*      cout<<"------------------------------------------------------"<<endl;
         cout<<"Theta"<<"\t"<<"Rho"<<"\t"<<"For XZ Projection"<<endl;
         cout<<"------------------------------------------------------"<<endl;
          
         for(auto i = 0; i < theta_pos_x.size(); i++)
         {  
           cout <<theta_pos_x[i]<<"\t"<<rho_pos_x[i]<<endl;            
         }
     */    
        
            for(int i= 1; i<= (nbintheta-1); i++)
        {    
             for(int j = 1; j<=(nbinrhoy - 1); j++)
           {      
             float bin_contYZ = houghYZ->GetBinContent(i,j);
             
             if(bin_contYZ >= 3)                                                      // checking for the peaks in YZ Hough space ---------------------
            {  thiy = -90.0 + (i-1)* theta_range/(nbintheta-1);
               rhojy = -1000.0 + (j-1)* rho_y/(nbinrhoy-1);
               th2yz = thiy;
               rh2yz = rhojy;

             if((th1yz - e) < th2yz && th2yz < (th1yz + e) && (rh1yz - erhy) < rh2yz && rh2yz < (rh1yz + erhy))
              { 
                if(bin_contYZ >= 4)
                 {    //cout<<bin_contYZ<<"bincont"<<endl;
                      peak_counterYZ++;                                           
                      theta_pos_y.insert(theta_pos_y.begin()+ cntyz, thiy);
                      rho_pos_y.insert(rho_pos_y.begin()+ cntyz, rhojy);
                      cntyz++;
                 }
                continue;
              }
             else
              {                        
               peak_counterYZ++;
               theta_pos_y.insert(theta_pos_y.begin()+ cntyz, thiy);
               rho_pos_y.insert(rho_pos_y.begin()+ cntyz, rhojy);
               cntyz++;
               th1yz = th2yz;
               rh1yz = rh2yz;
               continue;
              }
            }
          }
        }  
        
             
   //   cout<<"------------------------------------------------------"<<endl;
  /*       cout<<"Theta"<<"\t"<<"Rho"<<"\t"<<"For YZ Projection"<<endl;
         cout<<"------------------------------------------------------"<<endl;
      
           for(auto i = 0; i < theta_pos_y.size(); i++)
         {  
           cout <<theta_pos_y[i]<<"\t"<<rho_pos_y[i]<<endl;            
         } 
       */ 
         
         
             
 /*------------------------------------------------------------------------further analysis of the peaks------------------------------------------------------*/
     
      
     
          if(peak_counterXZ > 0 || peak_counterYZ > 0)
      { //4   
         Double_t th_newx = 0.0, th_newy = 0.0, dlt_rhoXZ, dlt_thetaXZ, dlt_rhoYZ, dlt_thetaYZ;
         int new_nbinrhoXZ, new_nbinthetaXZ, new_nbinrhoYZ, new_nbinthetaYZ;   
         vector<float> ntheta_pos_y;
         vector<float> nrho_pos_y;
         vector<float> ntheta_pos_x;
         vector<float> nrho_pos_x;
         
     //      for(int im = 0; im < theta_pos_x.size(); im++) 
     //      {
     //         cout<<theta_pos_x[im]<<endl; 
     //         //cout<<"Theta used to cal binsize = "<<th_newx<<endl;
     //      }

    //    cout<<"------------------------------------------------------"<<endl;
 /*      cout<<"New Theta"<<"\t"<<"New Rho"<<"\t"<<"For Zoomed XZ Projection"<<endl;
       cout<<"------------------------------------------------------"<<endl;
  */     
       
       for(int jt = 0; jt < 10; jt++)                                  // varying the new bin range---------------------------------------------
      { //5      
         int new_peak_counterXZ = 0, new_peak_counterYZ = 0, cntxz = 0, cntyz = 0;
         
         for(auto i = 0; i < theta_pos_x.size(); i++)
         {    
              th_newx = theta_pos_x[i];         
              dlt_rhoXZ   = delta_rho_calculation(th_newx);
              dlt_thetaXZ = delta_theta_calculation(th_newx);
              //dlt_rhoXZ = dlt_rhoXZ;               

              Double_t new_rho_minXZ = rho_pos_x[i] - 200.0 + jt * (dlt_rhoXZ/ 10.0);
              Double_t new_rho_maxXZ = rho_pos_x[i] + 200.0 + jt * (dlt_rhoXZ/ 10.0);
              Double_t new_theta_minXZ = theta_pos_x[i] - 4.0;
              Double_t new_theta_maxXZ = theta_pos_x[i] + 4.0;
              Double_t theta_rangeXZ = 8.0;
              Double_t rho_rangeXZ = 400.0;
             
              new_nbinrhoXZ   =  (rho_rangeXZ / dlt_rhoXZ);
              new_nbinthetaXZ =  (theta_rangeXZ / dlt_thetaXZ);
          
            TH2D* zoom_houghXZ;
            zoom_houghXZ = Hough_Acc(k, point_XZ , X, ZX, new_nbinthetaXZ, new_nbinrhoXZ, "histo_ZX_", rho_rangeXZ, th_thxz, theta_rangeXZ, new_theta_minXZ, new_theta_maxXZ, new_rho_minXZ, new_rho_maxXZ, dlt_thetaXZ, dlt_rhoXZ, bc); 
            
           //   for (int im = 0; im < theta_pos_x.size(); im++) 
           //   {
           //   //cout<<ntheta_pos_x[im]<<endl;
           //   cout<<"Theta used to cal binsize = "<<th_newx<<endl;
           //   }
        
             
              double thiXZ = 0.0, rhojXZ = 0.0;  
              float nth1xz = 90.0, nth2xz = 0.0, nerhx = 20.0, nrh1xz = 3000.0, nrh2xz = 0.0;

             for(int i= 1; i<=new_nbinthetaXZ; i++) 
           {    
              for(int j = 1; j<=new_nbinrhoXZ; j++)
             {      
              float bin_contxz = zoom_houghXZ-> GetBinContent(i,j);
                       
                if(bin_contxz >= 3)
              { thiXZ = new_theta_minXZ + i*dlt_thetaXZ;
                rhojXZ = new_rho_minXZ +  j*dlt_rhoXZ;
                nth2xz = thiXZ;
                nrh2xz = rhojXZ;
                if((nth1xz - e_zoom) < nth2xz && nth2xz < (nth1xz + e_zoom) && (nrh1xz - nerhx) < nrh2xz && nrh2xz < (nrh1xz + nerhx))
               { 
                  /*if(bin_contxz >= 4)
                 {  new_peak_counterXZ++; 
                    ntheta_pos_x.insert(ntheta_pos_x.begin()+ cntxz, thiXZ);
                    nrho_pos_x.insert(nrho_pos_x.begin()+ cntxz, rhojXZ);
                    cntxz++;
                  }*/
                 continue;
               }
                else
               {                        
                new_peak_counterXZ++; 
                ntheta_pos_x.insert(ntheta_pos_x.begin()+ cntxz, thiXZ);
                nrho_pos_x.insert(nrho_pos_x.begin()+ cntxz, rhojXZ);
                cntxz++;
                nth1xz = nth2xz;
                nrh1xz = nrh2xz;
                continue;
               }
              }
            }
          }  
           
           
             

            if(npeak_maxXZ < new_peak_counterXZ)
           {  
              thetaXZ.clear();
              rhoXZ.clear(); 
              npeak_maxXZ = new_peak_counterXZ;
              thetaXZ = ntheta_pos_x;
              rhoXZ   = nrho_pos_x;
         /*  cout<<"-----------------------------------------------"<<endl;
                 
              for(int i = 0; i < thetaXZ.size(); i++)
              { cout<<thetaXZ[i]<<"\t"<<rhoXZ[i]<<endl;
              } 
            cout<<"-----------------------------------------------"<<endl;
             
             f.cd();
            
            std::string coord_name_canvas2 =  "Hough_simu"+ std::to_string(k);
            TCanvas *c2 = new TCanvas("c2",coord_name_canvas2.c_str());
            c2->SetCanvasSize(1800, 1200);
            c2->SetWindowSize(1800, 1000);
            c2->SetName(coord_name_canvas2.c_str());
            
            c2->cd();
            gPad->SetGrid();
      
          std::string name_acc_xz2 = "ACC_XZ : Evt "+ std::to_string(k);
          zoom_houghXZ->SetName(name_acc_xz2.c_str());
          zoom_houghXZ->SetTitle(name_acc_xz2.c_str());
          zoom_houghXZ->GetXaxis()->SetLimits(new_theta_minXZ, new_theta_maxXZ);
          zoom_houghXZ->SetStats(0);
          gROOT->ForceStyle();
          gStyle->SetPalette(55);
          zoom_houghXZ->Draw("COLZ");
          c2->Modified(); c2->Update();  
      
          c2->Write();
          c2->Clear(); */ 

           }
          
           
          
             
     
/*           for(auto i = 0; i < ntheta_pos_x.size(); i++)
          {  
             cout <<ntheta_pos_x[i]<<"\t"<<nrho_pos_x[i]<<endl;            
          }
         
          
            
            f.cd();
            
            std::string coord_name_canvas2 =  "Hough_simu"+ std::to_string(k);
            TCanvas *c2 = new TCanvas("c2",coord_name_canvas2.c_str());
            c2->SetCanvasSize(1800, 1200);
            c2->SetWindowSize(1800, 1000);
            c2->SetName(coord_name_canvas2.c_str());
            
            c2->cd();
            gPad->SetGrid();
      
          std::string name_acc_xz2 = "ACC_XZ : Evt "+ std::to_string(k);
          zoom_houghXZ->SetName(name_acc_xz2.c_str());
          zoom_houghXZ->SetTitle(name_acc_xz2.c_str());
          zoom_houghXZ->GetXaxis()->SetLimits(new_theta_minXZ, new_theta_maxXZ);
          zoom_houghXZ->SetStats(0);
          gROOT->ForceStyle();
          gStyle->SetPalette(55);
          zoom_houghXZ->Draw("COLZ");
          c2->Modified(); c2->Update();
      
          c2->Write();
          c2->Clear();
   */      
         delete zoom_houghXZ;            
         
         }
      
    //   cout<<"-------------------------------------------------------------------------------------------------------------------"<<endl;
         
    //   cout<<"------------------------------------------------------"<<endl;
    //   cout<<"New Theta"<<"\t"<<"New Rho"<<"\t"<<"For Zoomed YZ Projection"<<endl;
    //   cout<<"------------------------------------------------------"<<endl;
         
       
          for(auto i = 0; i < theta_pos_y.size(); i++)
         {             
             th_newy = theta_pos_y[i];
            // cout<<theta_pos_y[i]<<endl;
             dlt_rhoYZ   = delta_rho_calculation(th_newy);
             dlt_thetaYZ = delta_theta_calculation(th_newy);
             //dlt_rhoYZ = dlt_rhoYZ - 0.02*dlt_rhoYZ;

             Double_t new_rho_minYZ = rho_pos_y[i] - 100.0 + jt * (dlt_rhoYZ/ 10.0);
             Double_t new_rho_maxYZ = rho_pos_y[i] + 100.0 + jt * (dlt_rhoYZ/ 10.0);
             Double_t new_theta_minYZ = theta_pos_y[i] - 4.0;
             Double_t new_theta_maxYZ = theta_pos_y[i] + 4.0;
             Double_t theta_rangeYZ = 8.0;
             Double_t rho_rangeYZ = 200.0;
             
             
             new_nbinrhoYZ   =  (rho_rangeYZ / dlt_rhoYZ);
             new_nbinthetaYZ =  (theta_rangeYZ / dlt_thetaYZ); 
             
                       
            TH2D* zoom_houghYZ;
            zoom_houghYZ = Hough_Acc(k, point_YZ , Y, ZY, new_nbinthetaYZ, new_nbinrhoYZ, "histo_ZY_", rho_rangeYZ, th_thyz, theta_rangeYZ, new_theta_minYZ, new_theta_maxYZ, new_rho_minYZ, new_rho_maxYZ, dlt_thetaYZ, dlt_rhoYZ, bc);
            
           
        
              float thiYZ = 0.0, rhojYZ = 0.0, nth1yz = 90.0, nth2yz = 0.0, nerhy = 10.0, nrh1yz = 1200.0, nrh2yz = 0.0;  
          
             for(int i= 1; i<=new_nbinthetaYZ; i++)
           {    
              for(int j = 1; j<=new_nbinrhoYZ; j++)
             {      
              float bin_contyz = zoom_houghYZ->GetBinContent(i,j);
                        
                if(bin_contyz >= 3)
              { thiYZ = new_theta_minYZ + i*dlt_thetaYZ;
                rhojYZ = new_rho_minYZ +  j*dlt_rhoYZ;
                nth2yz = thiYZ;
                nrh2yz = rhojYZ;
                if((nth1yz - e_zoom) < nth2yz && nth2yz < (nth1yz + e_zoom) && (nrh1yz - nerhy) < nrh2yz && nrh2yz < (nrh1yz + nerhy))
               {  
                 /*  if(bin_contyz >= 4)
                 {   
                  new_peak_counterYZ++;
                  ntheta_pos_y.insert(ntheta_pos_y.begin()+ cntyz, thiYZ);
                  nrho_pos_y.insert(nrho_pos_y.begin()+ cntyz, rhojYZ);
                  cntyz++; }*/
                  continue;
               }
                else
               {           
                new_peak_counterYZ++;
                ntheta_pos_y.insert(ntheta_pos_y.begin()+ cntyz, thiYZ);
                nrho_pos_y.insert(nrho_pos_y.begin()+ cntyz, rhojYZ);
                cntyz++;
                nth1yz = nth2yz;
                nrh1yz = nrh2yz;
                continue;
               }
              }
             }
          }  
              
     /*            f.cd();
            
            std::string coord_name_canvas3 =  "Hough_simu"+ std::to_string(k);
            TCanvas *c3 = new TCanvas("c3",coord_name_canvas3.c_str());
            c3->SetCanvasSize(1800, 1200);
            c3->SetWindowSize(1800, 1000);
            c3->SetName(coord_name_canvas3.c_str()); 
             
            c3->cd();
            gPad->SetGrid();
      
          std::string name_acc_yz2 = "ACC_YZ : Evt "+ std::to_string(k);
          zoom_houghYZ->SetName(name_acc_yz2.c_str());
          zoom_houghYZ->SetTitle(name_acc_yz2.c_str());
          zoom_houghYZ->GetXaxis()->SetLimits(new_theta_minYZ, new_theta_maxYZ);
          zoom_houghYZ->SetStats(0);
          gROOT->ForceStyle();
          gStyle->SetPalette(55);
          zoom_houghYZ->Draw("COLZ");
          c3->Modified(); c3->Update();
      
          c3->Write();
          c3->Clear();   */ 
                      
             if(npeak_maxYZ < new_peak_counterYZ)
           {  
              thetaYZ.clear();
              rhoYZ.clear();
              npeak_maxYZ = new_peak_counterYZ;
              thetaYZ = ntheta_pos_y;
              rhoYZ   = nrho_pos_y; 
           }
           
          
           
     //      for(auto i = 0; i < ntheta_pos_y.size(); i++)
     //     {  
     //        cout <<ntheta_pos_y[i]<<"\t"<<nrho_pos_y[i]<<endl;            
     //     }
         
           
            
         /*   f.cd();
            
            std::string coord_name_canvas3 =  "Hough_simu"+ std::to_string(k);
            TCanvas *c3 = new TCanvas("c3",coord_name_canvas3.c_str());
            c3->SetCanvasSize(1800, 1200);
            c3->SetWindowSize(1800, 1000);
            c3->SetName(coord_name_canvas3.c_str()); 
             
            c3->cd();
            gPad->SetGrid();
      
          std::string name_acc_yz2 = "ACC_YZ : Evt "+ std::to_string(k);
          zoom_houghYZ->SetName(name_acc_yz2.c_str());
          zoom_houghYZ->SetTitle(name_acc_yz2.c_str());
          zoom_houghYZ->GetXaxis()->SetLimits(new_theta_minYZ, new_theta_maxYZ);
          zoom_houghYZ->SetStats(0);
          gROOT->ForceStyle();
          gStyle->SetPalette(55);
          zoom_houghYZ->Draw("COLZ");
          c3->Modified(); c3->Update();
      
          c3->Write();
          c3->Clear();*/
         delete zoom_houghYZ;
         }
                 ntheta_pos_x.clear();
                 nrho_pos_x.clear();
                 ntheta_pos_y.clear();
                 nrho_pos_y.clear();
           } //5 
         } //4    
             
             theta_pos_x.clear();
             rho_pos_x.clear();
             theta_pos_y.clear();
             rho_pos_y.clear();

             delete houghXZ;
             delete houghYZ;
             
        } //3    
             
            
      
/*------------------------------------------------------------------------------- Plotting the results ------------------------------------------------------------*/
            
         /*   for(int i = 0; i < x_size/2; i++)
           { cout<<"-----------------------------------------------"<<endl;
             cout<<X[i]<<"\t"<<ZX[i]<<endl;
             cout<<"-----------------------------------------------"<<endl;
           } 
         */
            
            const int array_sizeXZ = dim[0];
            const int array_sizeYZ = dim[1];
            vector<float> dx;
            vector<float> dy;
            vector<float> dz1;
            vector<float> dz2;
            
            for(int i = 0; i < array_sizeXZ; i++)
            {
                dx.push_back(1.32);
                dz1.push_back(0.55);
            }
            
             for(int i = 0; i < array_sizeYZ; i++)
            {
                dy.push_back(1.32);
                dz2.push_back(0.55);
            }
            
      /*      f.cd();
            std::string coord_name_canvas2 =  "Hough_simu"+ std::to_string(k);
            TCanvas *c2 = new TCanvas("c2",coord_name_canvas2.c_str());
            c2->SetCanvasSize(1800, 1200);
            c2->SetWindowSize(1800, 1000);
            c2->Divide(2,2);
            c2->SetName(coord_name_canvas2.c_str());
            
            c2->cd(1);
            gPad->SetGrid();
            gROOT->ForceStyle(); // What is this used for?
            TGraphErrors *grxz;
            grxz = new TGraphErrors(dim[0], &xpos[0], &zxpos[0], &dx[0], &dz1[0]);
            std::string name_xz = "XZ : Evt "+ std::to_string(k);
            grxz->SetName(name_xz.c_str());
            grxz->SetTitle(name_xz.c_str());
            grxz->GetXaxis()->SetTitle("X pos of hits");
            grxz->GetYaxis()->SetTitle("Z pos of hits");
            grxz->SetLineColor(kRed);
            grxz->SetMarkerColor(kRed);
            grxz->SetMarkerSize(1.0);
            grxz->SetMarkerStyle(20);
            grxz->GetYaxis()->SetLimits(-500,500);
            grxz->GetXaxis()->SetLimits(-2250, 2250);
            grxz->Draw("AP*");
            c2->Modified(); c2->Update();
            
            
            c2->cd(1);
            gPad->SetGrid();     
            gROOT->ForceStyle(); 
            TMultiGraph *gr_tr_xz;
            cout<<X[0]<<"look"<<endl;
            gr_tr_xz = Plot_TracksXZ(thetaXZ, rhoXZ, X, ZX, point_XZ);                    // plotting tracks for the XZ case-----------------------
            gr_tr_xz->Draw("L");
            c2->Modified(); c2->Update();
            
            c2->cd(3);
            gPad->SetGrid();
            gROOT->ForceStyle();
            TGraphErrors *gryz;
            gryz = new TGraphErrors(dim[1], &ypos[0], &zypos[0], &dy[0], &dz2[0]);
            std::string name_yz = "YZ : Evt "+ std::to_string(k);
            gryz->SetName(name_yz.c_str());
            gryz->SetTitle(name_yz.c_str());
            gryz->GetXaxis()->SetTitle("Y pos of hits");
            gryz->GetYaxis()->SetTitle("Z pos of hits");
            gryz->SetLineColor(kRed);
            gryz->SetMarkerColor(kRed);
            gryz->SetMarkerSize(1.0);
            gryz->SetMarkerStyle(20);
            gryz->GetYaxis()->SetLimits(-500,500);
            gryz->GetXaxis()->SetLimits(-1000, 1000);
            gryz->Draw("AP*");
            c2->Modified(); c2->Update();
       
            c2->cd(3);
            gPad->SetGrid();     
            gROOT->ForceStyle(); 
            TMultiGraph *gr_tr_yz;
        
            gr_tr_yz = Plot_TracksYZ(thetaYZ, rhoYZ, Y, ZY, point_YZ);                   // plotting tracks for the YZ case-----------------------  
            gr_tr_yz->Draw("L");
            c2->Modified(); c2->Update();
           
            c2->cd(2);
            gPad->SetGrid();
       
       
            std::string name_acc_xz = "ACC_XZ : Evt "+ std::to_string(k);
            houghXZ->SetName(name_acc_xz.c_str());
            houghXZ->SetTitle(name_acc_xz.c_str());
            houghXZ->GetXaxis()->SetLimits(-90.4, 90.4);
            houghXZ->SetStats(0);
            gROOT->ForceStyle();
            gStyle->SetPalette(55);
            houghXZ->Draw("COLZ");
            c2->Modified(); c2->Update();
      

            c2->cd(4);
            gPad->SetGrid();
      
            std::string name_acc_yz = "ACC_YZ : Evt "+ std::to_string(k);
            houghYZ->SetName(name_acc_yz.c_str());
            houghYZ->SetTitle(name_acc_yz.c_str());
            houghYZ->GetXaxis()->SetLimits(-90.4, 90.4);
            houghYZ->SetStats(0);
            gROOT->ForceStyle();
            gStyle->SetPalette(55);
            houghYZ->Draw("COLZ");
            c2->Modified(); c2->Update();
      
            c2->Write();
            c2->Clear(); */
            
         /*------------------------------------------------------------plotting true tracks-----------------------------------*/   
            f.cd();
            std::string coord_name_canvas3 =  "Hough_simu"+ std::to_string(k);
            TCanvas *c3 = new TCanvas("c3",coord_name_canvas3.c_str());
            c3->SetCanvasSize(1800, 1200);
            c3->SetWindowSize(1800, 1000);
            c3->Divide(2,2);
            c3->SetName(coord_name_canvas3.c_str());
            
            c3->cd(1);
            gPad->SetGrid();
            gROOT->ForceStyle(); // What is this used for?
            TGraphErrors *grxz2;
            grxz2 = new TGraphErrors(dim[0], &xpos[0], &zxpos[0], &dx[0], &dz1[0]);
            std::string name_xz2 = "XZ : Evt "+ std::to_string(k);
            grxz2->SetName(name_xz2.c_str());
            grxz2->SetTitle(name_xz2.c_str());
            grxz2->GetXaxis()->SetTitle("X pos of hits");
            grxz2->GetYaxis()->SetTitle("Z pos of hits");
            grxz2->SetLineColor(kRed);
            grxz2->SetMarkerColor(kRed);
            grxz2->SetMarkerSize(1.0);
            grxz2->SetMarkerStyle(20);
            grxz2->GetYaxis()->SetLimits(-500,500);
            grxz2->GetXaxis()->SetLimits(-2250, 2250);
            grxz2->Draw("AP*");
            c3->Modified(); c3->Update();
            
            
            c3->cd(1);
            gPad->SetGrid();     
            gROOT->ForceStyle(); 
            TMultiGraph *gr_tr_xz;
            gr_tr_xz = Plot_TracksXZ(thetaXZ, rhoXZ, X, ZX, point_XZ);                    // plotting tracks for the XZ case-----------------------
            gr_tr_xz->Draw("L");
            c3->Modified(); c3->Update();
            
            c3->cd(3);
            gPad->SetGrid();
            gROOT->ForceStyle();
            TGraphErrors *gryz2;
            gryz2 = new TGraphErrors(dim[1], &ypos[0], &zypos[0], &dy[0], &dz2[0]);
            std::string name_yz2 = "YZ : Evt "+ std::to_string(k);
            gryz2->SetName(name_yz2.c_str());
            gryz2->SetTitle(name_yz2.c_str());
            gryz2->GetXaxis()->SetTitle("Y pos of hits");
            gryz2->GetYaxis()->SetTitle("Z pos of hits");
            gryz2->SetLineColor(kRed);
            gryz2->SetMarkerColor(kRed);
            gryz2->SetMarkerSize(1.0);
            gryz2->SetMarkerStyle(20);
            gryz2->GetYaxis()->SetLimits(-500,500);
            gryz2->GetXaxis()->SetLimits(-1000, 1000);
            gryz2->Draw("AP*");
            c3->Modified(); c3->Update();
       
            c3->cd(3);
            gPad->SetGrid();     
            gROOT->ForceStyle(); 
            TMultiGraph *gr_tr_yz;
            gr_tr_yz = Plot_TracksXZ(thetaYZ, rhoYZ, Y, ZY, point_YZ);                   // plotting tracks for the YZ case-----------------------  
            gr_tr_yz->Draw("L");
            c3->Modified(); c3->Update();
            
            c3->cd(2);
            gPad->SetGrid();
            gROOT->ForceStyle(); // What is this used for?
            TGraphErrors *grxz3;
            grxz3 = new TGraphErrors(dim[0], &xpos[0], &zxpos[0], &dx[0], &dz1[0]);
            std::string name_xz3 = "XZ : Evt "+ std::to_string(k);
            grxz3->SetName(name_xz3.c_str());
            grxz3->SetTitle(name_xz3.c_str());
            grxz3->GetXaxis()->SetTitle("X pos of hits");
            grxz3->GetYaxis()->SetTitle("Z pos of hits");
            grxz3->SetLineColor(kRed);
            grxz3->SetMarkerColor(kRed);
            grxz3->SetMarkerSize(1.0);
            grxz3->SetMarkerStyle(20);
            grxz3->GetYaxis()->SetLimits(-500,500);
            grxz3->GetXaxis()->SetLimits(-2250, 2250);
            grxz3->Draw("AP*");
            c3->Modified(); c3->Update();
            
            c3->cd(2);
            gPad->SetGrid();     
            gROOT->ForceStyle(); 
            TMultiGraph *gr_trtr_xz;
            gr_trtr_xz = Plot_True_Tracks(genc_x, genc_z, genp_x, genp_z, gen_sizeofarr);                    // plotting true tracks for the XZ case-----------------------
            gr_trtr_xz->Draw("L");
            c3->Modified(); c3->Update();
            
            c3->cd(4);
            gPad->SetGrid();
            gROOT->ForceStyle();
            TGraphErrors *gryz3;
            gryz3 = new TGraphErrors(dim[1], &ypos[0], &zypos[0], &dy[0], &dz2[0]);
            std::string name_yz3 = "YZ : Evt "+ std::to_string(k);
            gryz3->SetName(name_yz3.c_str());
            gryz3->SetTitle(name_yz3.c_str());
            gryz3->GetXaxis()->SetTitle("Y pos of hits");
            gryz3->GetYaxis()->SetTitle("Z pos of hits");
            gryz3->SetLineColor(kRed);
            gryz3->SetMarkerColor(kRed);
            gryz3->SetMarkerSize(1.0);
            gryz3->SetMarkerStyle(20);
            gryz3->GetYaxis()->SetLimits(-500,500);
            gryz3->GetXaxis()->SetLimits(-1000, 1000);
            gryz3->Draw("AP*");
            c3->Modified(); c3->Update();
             
            c3->cd(4);
            gPad->SetGrid();     
            gROOT->ForceStyle(); 
            TMultiGraph *gr_trtr_yz;
            gr_trtr_yz = Plot_True_Tracks(genc_y, genc_z, genp_y, genp_z, gen_sizeofarr);                   // plotting true tracks for the YZ case-----------------------  
            gr_trtr_yz->Draw("L");
            c3->Modified(); c3->Update();
      
            c3->Write();
            c3->Clear();
           
       //------------------------------------------------------------------------Storing data in a TTree---------------------------------------------------------//       
             
       
          vector<float> vec_xz_hough, vec_yz_hough; 
          int theta_sizeXZ, theta_sizeYZ;
            
          vec_xz_hough = Track_info(thetaXZ, rhoXZ, X, ZX, point_XZ);
          theta_sizeXZ = vec_xz_hough.back();
          ntracksXZ_hough = theta_sizeXZ;
          vec_xz_hough.pop_back();
          
          for(int i = (vec_xz_hough.size() - 2*theta_sizeXZ); i < vec_xz_hough.size(); i++)
          {
               Th_RhXZ_hough.push_back(vec_xz_hough[i]);
           }
          vec_xz_hough.erase(vec_xz_hough.begin() + (vec_xz_hough.size() - 2*theta_sizeXZ), vec_xz_hough.end()); 
          
          for(int i = 0; i < (vec_xz_hough.size()/3); i++)
          {     
                 Co_X_hough.push_back(vec_xz_hough[i]);
                 Co_ZX.push_back(vec_xz_hough[(vec_xz_hough.size()/3) + i]); 
                 Co_X.push_back(vec_xz_hough[2 *(vec_xz_hough.size()/3) + i]);    
          }
          
                    
          vec_yz_hough = Track_info(thetaYZ, rhoYZ, Y, ZY, point_YZ);
          theta_sizeYZ = vec_yz_hough.back();
          ntracksYZ_hough = theta_sizeYZ;
          vec_yz_hough.pop_back();
          
          for(int i = (vec_yz_hough.size() - 2*theta_sizeYZ); i < vec_yz_hough.size(); i++)
          {
               Th_RhYZ_hough.push_back(vec_yz_hough[i]);
           }
         vec_yz_hough.erase(vec_yz_hough.begin() + (vec_yz_hough.size() - 2*theta_sizeYZ), vec_yz_hough.end()); 
          
          for(int i = 0; i < (vec_yz_hough.size()/3); i++)
          {     
                 Co_Y_hough.push_back(vec_yz_hough[i]);
                 Co_ZY.push_back(vec_yz_hough[(vec_yz_hough.size()/3) + i]); 
                 Co_Y.push_back(vec_yz_hough[2 *(vec_yz_hough.size()/3) + i]);   
          }
          
      //------------------------------------------geninfo-----------------------------//    
         vector<float> vec_xz_gen, vec_yz_gen;
          int th_sizeXZ, th_sizeYZ;
            
          vec_xz_gen = True_Track_info(genc_x, genc_z, genp_x, genp_z, gen_sizeofarr);
          th_sizeXZ = vec_xz_gen.back();
          ntracksXZ_gen = th_sizeXZ;
          vec_xz_gen.pop_back();
          
          for(int i = 0; i < vec_xz_gen.size(); i++)
          {
               Th_RhXZ_gen.push_back(vec_xz_gen[i]);
           }
                    
          vec_yz_gen = True_Track_info(genc_y, genc_z, genp_y, genp_z, gen_sizeofarr);
          th_sizeYZ = vec_yz_gen.back();
          ntracksYZ_gen = th_sizeYZ;
          vec_yz_gen.pop_back();
          
          for(int i = 0; i < vec_yz_gen.size(); i++)
          {
               Th_RhYZ_gen.push_back(vec_yz_gen[i]);
           }
                  
          Ev_Id = k;
                  
          treeOut -> Fill();        
             
       int nhits_rejected = 0;
       nhits_rejected = no_hits;
       cout<<"Total number of hits that were rejected : "<<nhits_rejected<<endl;
  
     xpos.clear();
     zxpos.clear();
     ypos.clear();
     zypos.clear();
  } //2

   f.cd();
   treeOut->Write();
   f.Close();
   return 0;
}  //1




