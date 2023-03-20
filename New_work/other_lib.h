
using namespace std;
//------------------------------------------------------Hough Parametrization----------------------------------------------------//


Double_t Hough_rho(Double_t coord, Double_t z, Double_t theta)
{
  Double_t theta_rad =(M_PI *theta)/(180);  // converting from degree to radian
  Double_t c_pos = coord;   // coordinate position
  Double_t z_pos = z;      
  Double_t rho = c_pos * cos(theta_rad) + z_pos * sin(theta_rad);
  return rho;
}

//---------------------------------------------Hough Accumulator (or Accumulator Matrix)-----------------------------------------------//


TH2D* Hough_Acc(int nevt , int nhit_coord , Double_t *coord, Double_t *z , int number_bin_theta, int number_bin_rho, std::string coord_name, int rho, Double_t theta_thresh, Double_t theta_range, Double_t theta_min, Double_t theta_max, Double_t rho_min, Double_t rho_max, Double_t bin_size_theta, Double_t bin_size_rho, Double_t bincent)
{
  //Double_t bin_size_theta = (theta_range/(number_bin_theta-1)); 
  //Double_t bin_size_rho   = (rho/(number_bin_rho-1)) ; // if we take rho to be in cm
 

  std::string name_acc = "Acc_evt_" + coord_name + std::to_string(nevt) +"_hits_"+ std::to_string(nhit_coord);
  TH2D* Acc = new TH2D(name_acc.c_str(),"Acc", number_bin_theta,theta_min,theta_max ,number_bin_rho,rho_min,rho_max);
  Double_t rh = 0.0;
  Double_t theta= 0.0;

  for(int nhit = 0 ; nhit < nhit_coord ; nhit++)
  {
    for(int theta_deg = 0; theta_deg <= number_bin_theta; theta_deg++) // theta in degrees (-90, 90)
    {
      theta = theta_min + (bin_size_theta * (theta_deg + bincent)); 
      rh = Hough_rho(coord[nhit], z[nhit], theta);
     if(theta>= -1.0*theta_thresh and theta<= theta_thresh)
      {
        Acc->Fill(theta , rh);
      }
     }
   }
    return(Acc);
 }


//----------------------------------------------Removing identical points from the hough analysis-----------------------------------------//


vector<float> remove_same_points(vector<float> &xpos, vector<float> &zxpos, vector<float> &co_pos)
{     
      vector<float> cxpos;
      vector<float> x_c;
      vector<float> z_c;
      vector<float> c_c; 
      vector<float> temp;
      vector<float> temp_v;  
      x_c.clear();
      z_c.clear();
      int index;
      int same_coord_counter;
      map<int, vector<float>> map_unique_points;
      float ep = 0.01;
       cxpos.push_back(xpos[0]);
       cxpos.push_back(zxpos[0]);
       cxpos.push_back(co_pos[0]);
       map_unique_points[0] = cxpos; 
       cxpos.clear();
       for(auto i = 1; i < xpos.size(); i++)
       {       
            same_coord_counter = 0;
          for(auto itr = map_unique_points.begin(); itr != map_unique_points.end(); ++itr)
          {        //index = itr->first;   
                       temp = itr->second;
              
             if(((temp[0] - ep)  < xpos[i] && xpos[i] < (temp[0] + ep)) and ((temp[1] - ep) < zxpos[i] && zxpos[i] < (temp[1] + ep)) and ((temp[2] - ep) < co_pos[i] && co_pos[i] < (temp[2] + ep)))
              {  
                 same_coord_counter++;  temp.clear();        
              }              
              }
            if( same_coord_counter > 0 )
            {   continue; }
            else
            {        temp_v.push_back(xpos[i]);
                     temp_v.push_back(zxpos[i]);
                     temp_v.push_back(co_pos[i]);
                     map_unique_points[i] = temp_v; 
                     temp_v.clear(); 
                     temp.clear(); }
                                  
       }   
         temp.clear();
         for(auto itr = map_unique_points.begin(); itr != map_unique_points.end(); ++itr)
              {       temp = itr -> second;
                      x_c.push_back(temp[0]);
                      z_c.push_back(temp[1]);
                      c_c.push_back(temp[2]);
                      temp.clear();
               }          

     x_c.insert(x_c.end(), z_c.begin(), z_c.end());
     x_c.insert(x_c.end(), c_c.begin(), c_c.end());
     return x_c;
}



//-------------------------------------------------------------- Calculating appropriate binsizes for the corresponding peaks--------------------------------------------//

Double_t delta_rho_calculation(Double_t theta)
 { 
   
   Double_t strip = 2.64;
   Double_t delta_rho = strip * cos(theta * M_PI/180.0); 
   
   return delta_rho;
 }

Double_t delta_theta_calculation(Double_t theta)
 { 
   Double_t z = 150.0;
   Double_t strip = 2.64;
   Double_t rh1 = sqrt( z*z + (pow((z/tan((90.0 - theta) * M_PI/180.0)),2)));
   Double_t rh2 = sqrt( z*z + (pow((z/tan((90.0 - theta) * M_PI/180.0)) + strip,2))); 
   Double_t delta_theta = acos((rh1*rh1 + rh2*rh2 - strip*strip)/(2.0 * rh1 * rh2));
   Double_t delta_theta_degree = delta_theta * 180.0/M_PI;
   
   return delta_theta_degree;
 }

/*--------------------------------------------------------------------------------------------Plotting the tracks-----------------------------------------------------------------*/


TMultiGraph* Plot_TracksXZ(vector<float> &ntheta_x, vector<float> &nrho_x, Double_t *x_c, Double_t *z_c, int nhit_cx)
{        
         TMultiGraph  *mg1  = new TMultiGraph();
         float th_degree, rho, th_radian, x, temp_var1 = 0.0, temp_var2 = 0.0, temp_var3 = 0.0, e = 5.0, e2 = 15.0, ep = 4.0, key, key1, index;
         vector<float> xc, nxc;
         vector<float> xcc1, xcc2;
         vector<float> coord_x, ncoord_x;
         vector<float> zc1, nzc1;
         vector<float> zcc1, zcc2;
         vector<float> co_x1;
         vector<float> co_z1;
         vector<float> co_x2;
         vector<float> co1;
         vector<float> co2;
         vector<float> cc;
         vector<float> zc2, temp2, temp;
         vector<float> new_theta;
         vector<float> new_rho;
         multimap<float, vector<float>> mp_c;
         vector<float> v_c;
 
         float zc = 0.0, diff, delta_z, strip = 2.64, sigma = 2.64/sqrt(12.0), chi2_1, chi2_2, dummy1_chi2 = 0.0, dummy2_chi2 = 0.0, dummy_theta = 100.0, dth = 10.0, dummy_rho = 3000.0, drh=80.0, dummy_diff, delta_x = 5.0 * strip;
         int cnt, cnt2, cn = -1, cn2, nc = 0, counter, v_size;
         new_theta.clear();
         new_rho.clear();
       
         for(int i = 0; i < ntheta_x.size(); i++)
         {  
             // if((35.0 - 2.0) < ntheta_x[i] && ntheta_x[i] < (35.0 + 2.0)) 
             // continue;
             // cout<<"theta = "<<ntheta_x[i]<<endl;
               
            cnt = 0; cnt2 = 0; chi2_1 = 0.0; chi2_2 = 0.0;
            xc.clear();
            zc1.clear();
            zc2.clear();
            coord_x.clear(); 
            nxc.clear();
            nzc1.clear();
            ncoord_x.clear();     
           
             
              for(int j = 0; j < nhit_cx; j++)
            {    
                 th_degree  = ntheta_x[i];
                 rho = nrho_x[i];
                 th_radian = th_degree * M_PI/180.0;
                 
                x = (rho/cos(th_radian))  -  z_c[j]*(sin(th_radian)/cos(th_radian));
                dummy_diff = 500;
                temp_var1 = 0.0; temp_var2 = 0.0; temp_var3 = 0.0;
               for(int k = 0; k < nhit_cx; k++)                        //---------This loop will make sure that only the nearest points get saved
              {   
                  if(((z_c[k] - ep) < z_c[j]) && (z_c[j] < (z_c[k] + ep)))
                  {    diff = abs(x_c[k] - x);
                     if((dummy_diff > diff)  && (diff < delta_x))
                    {  dummy_diff = diff;
                       temp_var1 = x;
                       temp_var2 = x_c[k];
                       temp_var3 = z_c[j];
                     }
                     else
                     { continue;}
                   }   
                  else
                     { continue;} 
              }
                 
                 if(temp_var1 != 0.0 && temp_var2 != 0.0 && temp_var3 != 0.0)
                 {xc.insert(xc.begin() + cnt, temp_var1);
                 coord_x.insert(coord_x.begin() + cnt, temp_var2);        
                 zc1.insert(zc1.begin() + cnt, temp_var3);
                 zc2.insert(zc2.begin() + cnt, temp_var3);
                 cnt++; temp_var1 = 0.0; temp_var2 = 0.0; temp_var3 = 0.0;
                 }
           }  
             
          /* for(int v = 0; v < xc.size(); v++)
                {  cout<<xc[v]<<"\t"<<coord_x[v]<<"\t"<<zc1[v]<<"\t"<<zc2[v]<<"\t"<<ntheta_x[i]<<"\t"<<nrho_x[i]<<endl;
                }  
           cout<<"------------"<<ntheta_x[i]<<"----------------------"<<endl;
           */
           
           co1 = remove_same_points(xc, zc1, coord_x);
           
           v_size = co1.size()/3;
         /*  for(int v = 0; v < v_size; v++)
                {  cout<<co1[v]<<"\t"<<co1[v_size+v]<<"\t"<<co1[2*v_size+ v]<<"\t"<<endl;
                }  
           cout<<"------------"<<ntheta_x[i]<<"----------------------"<<endl; 
           */
           //xc.clear(); zc1.clear(); coord_x.clear();
           key = ntheta_x[i] + nrho_x[i];
           for(int v = 0; v < v_size; v++)
          {       
                  nxc.push_back(co1[v]); v_c.push_back(co1[v]);
                  nzc1.push_back(co1[v_size + v]); v_c.push_back(co1[v_size + v]);
                  ncoord_x.push_back(co1[2*v_size + v]); v_c.push_back(co1[2*v_size + v]);
                  mp_c.insert(make_pair( key, v_c)); 
                  v_c.clear();
           }
           
          /*  for(auto itr = mp_c.begin(); itr != mp_c.end(); ++itr)
              {        index = itr->first;                                          // -------  printing map
                       temp = itr->second;
                cout<<index<<"\t"<<temp[0]<<"\t"<<temp[2]<<"\t"<<temp[1]<<"\t"<<endl;       
                     temp.clear();
               }      
       
       cout<<"-----------map for theta = "<<ntheta_x[i]<<"----------------------"<<endl;  
               */ 
               if(v_size > 3) 
              {   
                  for(int l = 0; l < v_size; l++) 
                  {chi2_2 = chi2_2 + pow((ncoord_x[l] - nxc[l]),2)/pow(sigma,2);}  
                  
                 if(((dummy_theta - dth) < ntheta_x[i] && ntheta_x[i] < (dummy_theta + dth)) and ((dummy_rho - drh) < nrho_x[i] && nrho_x[i] < (dummy_rho + drh)))
                  {      
                        if(dummy2_chi2 < chi2_2)
                        {   nc++;
                           continue;       
                         }                                
                        else 
                       {
                         dummy_theta = ntheta_x[i];
                         dummy_rho   = nrho_x[i];
                         new_theta.at(cn) = dummy_theta;
                         new_rho.at(cn) = dummy_rho; 
                         dummy2_chi2 = chi2_2;  
                        }
                   }
                   else
                   { 
                      counter = 0;
                     for(int g = 0; g < new_theta.size(); g++)
                      { 
                        if(((new_theta[g] - dth) < ntheta_x[i] && ntheta_x[i] < (new_theta[g] + dth)) and ((new_rho[g] - drh) < nrho_x[i] && nrho_x[i] < (new_rho[g] + drh)))
                         {
                            counter++;
                         }      
                      } 
                      if(counter > 0){     continue; }
                      else{  cn++;
                      dummy_theta = ntheta_x[i];
                      dummy_rho   = nrho_x[i]; 
                      new_theta.insert(new_theta.begin() + cn, dummy_theta);
                      new_rho.insert(new_rho.begin() + cn, dummy_rho); 
                      dummy2_chi2 = chi2_2; }     
                   }   
               }
              else if(v_size > 2)
              {    
                      for(int l = 0; l < v_size; l++) 
                  {chi2_1 = chi2_1 + pow((ncoord_x[l] - nxc[l]),2)/pow(sigma,2);} 
                  
                 if((dummy_theta - dth) < ntheta_x[i] && ntheta_x[i] < (dummy_theta + dth) && (dummy_rho - drh) < nrho_x[i] && nrho_x[i] < (dummy_rho + drh))
                  {      
                        if(dummy1_chi2 < chi2_1)
                        {   nc++;
                           continue;       
                         }                                
                        else 
                       {
                         dummy_theta = ntheta_x[i];
                         dummy_rho   = nrho_x[i];
                         new_theta.at(cn) = dummy_theta;
                         new_rho.at(cn) = dummy_rho; 
                         dummy1_chi2 = chi2_1;   
                        }
                   }
                   else
                   { counter = 0; 
                     for(int g = 0; g < new_theta.size(); g++)
                      { 
                        if((new_theta[g] - dth) < ntheta_x[i] && ntheta_x[i] < (new_theta[g] + dth) && (new_rho[g] - drh) < nrho_x[i] && nrho_x[i] < (new_rho[g] + drh))
                         {
                            counter++;
                         }      
                      }      
                      if(counter > 0){     continue; }
                      else{  cn++;
                     dummy_theta = ntheta_x[i];
                     dummy_rho   = nrho_x[i]; 
                     new_theta.insert(new_theta.begin() + cn, dummy_theta);
                     new_rho.insert(new_rho.begin() + cn, dummy_rho);
                     dummy1_chi2 = chi2_1; }       
                    }  
              
               }      
            
          }//closer of theta_loop
            
        
             for(int i = 0; i < new_theta.size(); i++)
            {        key1 = new_theta[i] + new_rho[i];
                 auto range = mp_c.equal_range(key1);
               for(auto it = range.first; it!=range.second; ++it)
               {   
                     temp2 = it->second;
                     xcc1.push_back(temp2[0]);
                     zcc1.push_back(temp2[1]);
                     cc.push_back(temp2[2]); 
                     temp2.clear();       
               }
               
             int dim1 = zcc1.size();
             TGraph* gr1 = new TGraph(dim1,&xcc1[0],&zcc1[0]);
             gr1->SetLineColor(4);   
             mg1->Add(gr1);
             
             xcc1.clear();
             zcc1.clear();
             cc.clear();
             }
           
           mp_c.clear();  
        return(mg1);   
}

TMultiGraph* Plot_True_Tracks(vector<float> &gnc_x, vector<float> &gnc_z, vector<float> &gnp_x, vector<float> &gnp_z, int gn_sizeofarr)
{        
         TMultiGraph  *mg3  = new TMultiGraph();
         
       float x_gen[2] = {0.0}, z_gen[2] = {0.0}, slope, intercept, c_x, c_z, p_x, p_z;
       x_gen[0] = -2250.0; // -1000.0
       x_gen[1] = 2250.0;  // 1000.0
       for(int l = 0; l < gn_sizeofarr; l++)
       {    
            c_x = gnc_x[l];
            c_z = gnc_z[l];
            p_x = gnp_x[l];
            p_z = gnp_z[l];
            slope = p_z/p_x;
            intercept = c_z - (slope * c_x); 
            z_gen[0] = (slope * x_gen[0]) + intercept;
            z_gen[1] = (slope * x_gen[1]) + intercept;
        
            TGraph* gr3 = new TGraph(2, x_gen, z_gen);
            gr3->SetLineColor(kGreen+2);   
            mg3->Add(gr3);
       }
       
      
      return(mg3);   
}


//--------------------------------------------------------------Calculating the Efficiency --------------------------------------------------------------//

vector<float> Track_info(vector<float> &ntheta_x, vector<float> &nrho_x, Double_t *x_c, Double_t *z_c, int nhit_cx)
{        
        
         float th_degree, rho, th_radian, x, temp_var1 = 0.0, temp_var2 = 0.0, temp_var3 = 0.0, e = 5.0, e2 = 15.0, ep = 4.0, key, key1, index;
         vector<float> xc, nxc;
        // vector<float> xcc1, xcc2;
         vector<float> coord_x, ncoord_x;
         vector<float> zc1, nzc1;
        // vector<float> zcc1, zcc2;
         vector<float> co1;
         vector<float> co2;
       //  vector<float> cc;
         vector<float> zc2, temp2, temp;
         vector<float> new_theta;
         vector<float> new_rho;
         multimap<float, vector<float>> mp_c;
         vector<float> v_c;
 
         float zc = 0.0, diff, delta_z, strip = 2.64, sigma = 2.64/sqrt(12.0), chi2_1, chi2_2, dummy1_chi2 = 0.0, dummy2_chi2 = 0.0, dummy_theta = 100.0, dth = 10.0, dummy_rho = 3000.0, drh=80.0, dummy_diff, delta_x = 4.0 * strip;
         int cnt, cnt2, cn = -1, cn2, nc = 0, counter, v_size;
         new_theta.clear();
         new_rho.clear();
       
         for(int i = 0; i < ntheta_x.size(); i++)
         {  cnt = 0; cnt2 = 0; chi2_1 = 0.0; chi2_2 = 0.0;
            xc.clear(); zc1.clear(); zc2.clear(); coord_x.clear(); nxc.clear(); nzc1.clear(); ncoord_x.clear(); //xcc1.clear(); zcc1.clear(); cc.clear();    
             
              // cout<<"initial hough theta = "<<ntheta_x[i]<<"ntheta size = "<<ntheta_x.size()<<endl;
              for(int j = 0; j < nhit_cx; j++)
            {    
                 th_degree  = ntheta_x[i];
                 rho = nrho_x[i];
                 th_radian = th_degree * M_PI/180.0;
                 
                x = (rho/cos(th_radian))  -  z_c[j]*(sin(th_radian)/cos(th_radian));
                dummy_diff = 500;
                temp_var1 = 0.0; temp_var2 = 0.0; temp_var3 = 0.0;
               for(int k = 0; k < nhit_cx; k++)                        //---------This loop will make sure that only the nearest points get saved
              {   
                  if(((z_c[k] - ep) < z_c[j]) && (z_c[j] < (z_c[k] + ep)))
                  {    diff = abs(x_c[k] - x);
                     if((dummy_diff > diff)  && (diff < delta_x))
                    {  dummy_diff = diff;
                       temp_var1 = x;
                       temp_var2 = x_c[k];
                       temp_var3 = z_c[j];
                     }
                    // else
                    // { continue;}
                   }   
                 // else
                   //  { continue;} 
              }
                 
                 if(temp_var1 != 0.0 && temp_var2 != 0.0 && temp_var3 != 0.0)
                 {xc.insert(xc.begin() + cnt, temp_var1);
                 coord_x.insert(coord_x.begin() + cnt, temp_var2);        
                 zc1.insert(zc1.begin() + cnt, temp_var3);
                 zc2.insert(zc2.begin() + cnt, temp_var3);
                 cnt++; temp_var1 = 0.0; temp_var2 = 0.0; temp_var3 = 0.0;
                 }
           }  
                
           co1 = remove_same_points(xc, zc1, coord_x);
           
           v_size = co1.size()/3;
           
           key = ntheta_x[i] + nrho_x[i];
           for(int v = 0; v < v_size; v++)
          {       
                  nxc.push_back(co1[v]); v_c.push_back(co1[v]);
                  nzc1.push_back(co1[v_size + v]); v_c.push_back(co1[v_size + v]);
                  ncoord_x.push_back(co1[2*v_size + v]); v_c.push_back(co1[2*v_size + v]);
                  mp_c.insert(make_pair( key, v_c)); 
                  v_c.clear();
           }
           
              // cout<<"--------------------"<<ntheta_x[i]<<"-------------------"<<endl;
              // for(int l = 0; l < v_size; l++) 
               //   {
                //    cout<<nxc[l]<<"\t"<<ncoord_x[l]<<endl; 
                    
               //   }  
                 
            
               if(v_size > 3) 
              {   
                  for(int l = 0; l < v_size; l++) 
                  {chi2_2 = chi2_2 + pow((ncoord_x[l] - nxc[l]),2)/pow(sigma,2);}  
                  
                 if(((dummy_theta - dth) < ntheta_x[i] && ntheta_x[i] < (dummy_theta + dth)) and ((dummy_rho - drh) < nrho_x[i] && nrho_x[i] < (dummy_rho + drh)))
                  {      
                        if(dummy2_chi2 < chi2_2)
                        {   nc++;
                           continue;       
                         }                                
                        else 
                       {
                         dummy_theta = ntheta_x[i];
                         dummy_rho   = nrho_x[i];
                         new_theta.at(cn) = dummy_theta;
                         new_rho.at(cn) = dummy_rho; 
                         dummy2_chi2 = chi2_2;  
                        }
                   }
                   else
                   { 
                      counter = 0;
                     for(int g = 0; g < new_theta.size(); g++)
                      { 
                        if(((new_theta[g] - dth) < ntheta_x[i] && ntheta_x[i] < (new_theta[g] + dth)) and ((new_rho[g] - drh) < nrho_x[i] && nrho_x[i] < (new_rho[g] + drh)))
                         {
                            counter++;
                         }      
                      } 
                      if(counter > 0){     continue; }
                      else{  cn++;
                      dummy_theta = ntheta_x[i];
                      dummy_rho   = nrho_x[i]; 
                      new_theta.insert(new_theta.begin() + cn, dummy_theta);
                      new_rho.insert(new_rho.begin() + cn, dummy_rho); 
                      dummy2_chi2 = chi2_2; }     
                   }   
               }
              else if(v_size > 2)
              {    
                      for(int l = 0; l < v_size; l++) 
                  {chi2_1 = chi2_1 + pow((ncoord_x[l] - nxc[l]),2)/pow(sigma,2);} 
                  
                 if((dummy_theta - dth) < ntheta_x[i] && ntheta_x[i] < (dummy_theta + dth) && (dummy_rho - drh) < nrho_x[i] && nrho_x[i] < (dummy_rho + drh))
                  {      
                        if(dummy1_chi2 < chi2_1)
                        {   nc++;
                           continue;       
                         }                                
                        else 
                       {
                         dummy_theta = ntheta_x[i];
                         dummy_rho   = nrho_x[i];
                         new_theta.at(cn) = dummy_theta;
                         new_rho.at(cn) = dummy_rho; 
                         dummy1_chi2 = chi2_1;   
                        }
                   }
                   else
                   { counter = 0; 
                     for(int g = 0; g < new_theta.size(); g++)
                      { 
                        if((new_theta[g] - dth) < ntheta_x[i] && ntheta_x[i] < (new_theta[g] + dth) && (new_rho[g] - drh) < nrho_x[i] && nrho_x[i] < (new_rho[g] + drh))
                         {
                            counter++;
                         }      
                      }      
                      if(counter > 0){     continue; }
                      else{  cn++;
                     dummy_theta = ntheta_x[i];
                     dummy_rho   = nrho_x[i]; 
                     new_theta.insert(new_theta.begin() + cn, dummy_theta);
                     new_rho.insert(new_rho.begin() + cn, dummy_rho);
                     dummy1_chi2 = chi2_1; }       
                    }  
              
               }      
            
          }//closer of theta_loop
            
       // cout<<"------------------theta---------------"<<endl;
          /*   for(int i = 0; i < new_theta.size(); i++)
            {        key1 = new_theta[i] + new_rho[i];
                 auto range = mp_c.equal_range(key1);
               for(auto it = range.first; it!=range.second; ++it)
               {   
                     temp2 = it->second;
                     xcc1.push_back(temp2[0]);
                     zcc1.push_back(temp2[1]);
                     cc.push_back(temp2[2]); 
                     temp2.clear();       
               }
             }
           */
           /* for(int v = 0; v < zcc1.size(); v++)
                {  cout<<xcc1[v]<<"\t"<<zcc1[v]<<endl;
                }  */
          // cout<<"-------------xcc and zcc--for theta = "<<new_theta[i]<<"\t"<<new_theta.size()<<"---------------------"<<endl; 
             
         /*  xcc1.insert(xcc1.end(), zcc1.begin(), zcc1.end());
           xcc1.insert(xcc1.end(), cc.begin(), cc.end());
           xcc1.insert(xcc1.end(), new_theta.begin(), new_theta.end());
           xcc1.insert(xcc1.end(), new_rho.begin(), new_rho.end());
           */
         //  int num_tracks = new_theta.size();
           new_theta.insert(new_theta.end(), new_rho.begin(), new_rho.end());
         //  cout<<"hough theta = "<<new_theta[0]<<endl;
           mp_c.clear();
             
        return new_theta;   
}

vector<float> True_Track_info(vector<float> &gnc_x, vector<float> &gnc_z, vector<float> &gnp_x, vector<float> &gnp_z, int gn_sizeofarr)
{            
       float x_gen[2] = {0.0}, z_gen[2] = {0.0}, slope, intercept, c_x, c_z, p_x, p_z, phi_radian, phi_degree, theta_gen, rho_gen;
       vector<float> th_gen, rh_gen;
       
       for(int l = 0; l < gn_sizeofarr; l++)
       {    
            c_x = gnc_x[l];
            c_z = gnc_z[l];
            p_x = gnp_x[l];
            p_z = gnp_z[l];
            slope = p_z/p_x;
            phi_radian = atan(slope);
            phi_degree = phi_radian * (180.0/M_PI);
            
            if(slope < 0)
            {
              theta_gen = 90.0 + phi_degree;  
            }
            else 
            {
              theta_gen = -90.0 + phi_degree;
            }
            
            rho_gen = c_x * cos(theta_gen) +  c_z * sin(theta_gen);
            th_gen.push_back(theta_gen);
            rh_gen.push_back(rho_gen);
            
           // cout<<"geninfo theta = "<<theta_gen<<endl;
        }
        
         th_gen.insert(th_gen.end(), rh_gen.begin(), rh_gen.end());
         th_gen.push_back(gn_sizeofarr);
         
      return th_gen;   
}


