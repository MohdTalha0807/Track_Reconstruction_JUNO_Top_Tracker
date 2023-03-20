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
    TFile *fileIn1 = new TFile("muon_265_user.root","READ");

TTree *tree_TTDigit=(TTree*)fileIn1->Get("TTDigit"); //Read TTDigit : Simulation file
int evtID_TTDigit;
tree_TTDigit->SetBranchAddress("evtID",&evtID_TTDigit);
int NChannels;
tree_TTDigit->SetBranchAddress("NTouchedChannel",&NChannels);
int Channel_ID[1000];
tree_TTDigit->SetBranchAddress("TB_DMchannel",&Channel_ID);
int Charge_ADC[1000];
tree_TTDigit->SetBranchAddress("TB_ADC",&Charge_ADC);
float npe[1000];
tree_TTDigit->SetBranchAddress("TB_pe",&npe);
float XcPos[1000];
tree_TTDigit->SetBranchAddress("TB_xcC",&XcPos);
float YcPos[1000];
tree_TTDigit->SetBranchAddress("TB_ycC",&YcPos);
float ZcPos[1000];
tree_TTDigit->SetBranchAddress("TB_zcC",&ZcPos);

int nentries = tree_TTDigit->GetEntries();

string root_ext = ".root";
string filename_root = outfile + root_ext;
TFile * Outfile	  = new TFile(filename_root.c_str(),"RECREATE");

Outfile->cd();
TTree *tree = new TTree("Hough_Tree", "Hough_Tree"); //Defining a Tree

double area_column = 2.64*(64*(3.0/2.0));
double area_row = 2.64*(64*(7.0/2.0));

Int_t max_column=3;
Int_t max_row=7;
Int_t max_layers=3;
Float_t charge_condition= 1.0/3.0; //1.5; ///0.333;
double pi_factor=3.141516;


for(int n=0; n<nentries;n++){//Open Nentries  : numbers of events

  tree_TTDigit->GetEntry(n);
  if(NChannels==0)continue;
  cout<<" Event number  "<<  n<<"  "<<"and nhits started with a number of hits equal to : "<<NChannels<<endl;
  if(NChannels>100)continue;

  //1. Saving hits larger than 1/3 p.e

  Int_t nhits_counted=100;
  Int_t lfchannel=0;
  Int_t lflayers[nhits_counted];
  Int_t lfpmt[nhits_counted];
  Int_t lfstrip[nhits_counted];
  Int_t lfrow[nhits_counted];
  Int_t lfcolumn[nhits_counted];
  Float_t lfnpe[nhits_counted];
  Float_t lfxpos[nhits_counted];
  Float_t lfypos[nhits_counted];
  Float_t lfzpos[nhits_counted];

  Int_t nhits_condition = 0;
  for(int initial_hits=0; initial_hits<NChannels; initial_hits++)
  {
     if(npe[initial_hits] > charge_condition){
       lfchannel = Channel_ID[initial_hits];
       lflayers[nhits_condition] = LayerFromJUNOid(lfchannel);
       lfpmt[nhits_condition] = PMTsideFromJUNOid(lfchannel);
       lfstrip[nhits_condition] = StripFromJUNOid(lfchannel);
       lfrow[nhits_condition]= RowFromJUNOid(lfchannel);
       lfcolumn[nhits_condition]= ColumnFromJUNOid(lfchannel);
       lfnpe[nhits_condition]= npe[initial_hits];
       lfxpos[nhits_condition]=XcPos[initial_hits];
       lfypos[nhits_condition]=YcPos[initial_hits];
       lfzpos[nhits_condition]=ZcPos[initial_hits];
       nhits_condition++;
     }
   }

   cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
   cout<<"Number of hits larger than 1/3 p.e = "<<nhits_condition<<endl;
   cout<<"The events saved are the following : "<<endl;
   cout<<"layers"<<"  "<<"column"<<"  "<<"row"<<"  "<<"pmt"<<"  "<<"strip"<<"  "<<"npe"<<" "<<"xpos"<<" "<<"ypos"<<" "<<"zpos"<<endl;
   cout<<"--------------------------------------------------------------------------------------------------------"<<endl;

   for(int new_hits=0; new_hits<nhits_condition; new_hits++)
   {
     cout<<lflayers[new_hits]<<"  "<<lfcolumn[new_hits]<<"  "<<lfrow[new_hits]<<"  "<<lfpmt[new_hits]<<"  "<<lfstrip[new_hits]<<"  "<<lfnpe[new_hits]<<" "<<lfxpos[new_hits]<<" "<<lfypos[new_hits]<<" "<<lfzpos[new_hits]<<endl;
     cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
   }

   //2. Check XY Coincidences
   int xcount[3];
   int ycount[3];
   int n0x=0, n1x=0, n2x=0, ni=0;
   int n0y=0, n1y=0, n2y=0, nj=0;
   int nt=0;

   //cout<<"Next Steps "<<endl;
   for(int new_hits=0; new_hits<nhits_condition; new_hits++)
   {
     for(int l=0; l< max_layers; l++)
     {
      for(int i=0; i <max_column; i++)
      {
       for (int j=0; j<max_row; j++)
       {
           if(lflayers[new_hits]==l && lfrow[new_hits] == j && lfcolumn[new_hits] == i)
           {
             if(lfpmt[new_hits]== 0 or lfpmt[new_hits]== 1 or lfpmt[new_hits]== 2 or lfpmt[new_hits]== 3 or lfpmt[new_hits]== 8 or lfpmt[new_hits]== 9 or lfpmt[new_hits]== 10 or lfpmt[new_hits]== 11)
             {
               if(lflayers[new_hits]==0){n0x++;}
               if(lflayers[new_hits]==1){n1x++;}
               if(lflayers[new_hits]==2){n2x++;}
               cout<<"X yes :"<<l<<" "<<i<<" "<<j<<" "<<lfpmt[new_hits]<<endl;
               ni++;
             }

             if(lfpmt[new_hits]== 4 or lfpmt[new_hits]== 5 or lfpmt[new_hits]== 6 or lfpmt[new_hits]== 7 or lfpmt[new_hits]== 12 or lfpmt[new_hits]== 13 or lfpmt[new_hits]== 14 or lfpmt[new_hits]== 15)
             {
               if(lflayers[new_hits]==0){n0y++;}
               if(lflayers[new_hits]==1){n1y++;}
               if(lflayers[new_hits]==2){n2y++;}
               cout<<"Y yes :"<<l<<" "<<i<<" "<<j<<" "<<lfpmt[new_hits]<<endl;
               nj++;
             }
             nt++;
            }
          }
        }
       }
     }

     cout<<" number of hits in x proj : "<<ni<<endl;
     cout<<" layer 0 : "<<n0x<<endl;
     cout<<" layer 1 : "<<n1x<<endl;
     cout<<" layer 2 : "<<n2x<<endl;
     cout<<" total number of hits in x proj : "<<n0x+n1x+n2x<<endl;

     cout<<" number of hits in y proj : "<<nj<<endl;
     cout<<" layer 0 : "<<n0y<<endl;
     cout<<" layer 1 : "<<n1y<<endl;
     cout<<" layer 2 : "<<n2y<<endl;
     cout<<" total number of hits in y proj : "<<n0y+n1y+n2y<<endl;

     cout<<" number of hits in x and y proj : "<<ni+nj<<endl;
     cout<<" total number of hits in x and y proj : "<<nt<<endl;

     //Defining hits for each layer proj
     xcount[0]=n0x;
     xcount[1]=n1x;
     xcount[2]=n2x;
     ycount[0]=n0y;
     ycount[1]=n1y;
     ycount[2]=n2y;
     int nh=0;

     Double_t X[ncounts];
     Double_t Y[ncounts];
     Double_t ZX[ncounts];
     Double_t ZY[ncounts];

     Double_t errorX[ncounts];
     Double_t errorY[ncounts];
     Double_t errorZX[ncounts];
     Double_t errorZY[ncounts];

     int numberofhitsin_x=0;
     int numberofhitsin_y=0;

     //3. Minimum 1 point per layer  (one event with 3 layers points)

     if((n0x>0 and n1x>0 and n2x>0 and n0y>0 and n1y>0 and n2y>0)>2){
       /*cout<<"Event number : "<<n<<endl;
       cout<<" Layers "<<" Column "<<" Row "<<" PMT "<<" Strip "<< " Charge "<<endl;
       cout<<"**************************************************************"<<endl;*/
       for(int new_hits=0; new_hits<nhits_condition; new_hits++)
       {
         if(lfpmt[new_hits]== 0 or lfpmt[new_hits]== 1 or lfpmt[new_hits]== 2 or lfpmt[new_hits]== 3 or lfpmt[new_hits]== 8 or lfpmt[new_hits]== 9 or lfpmt[new_hits]== 10 or lfpmt[new_hits]== 11)
         {
           X[numberofhitsin_x]=  lfxpos[new_hits];
           ZX[numberofhitsin_x]= lfzpos[new_hits];
           errorX[numberofhitsin_x] = 2.64/(sqrt(12.0));
           errorZX[numberofhitsin_x] = 1.06/(sqrt(12.0));
           numberofhitsin_x++;
         }
         if(lfpmt[new_hits]== 4 or lfpmt[new_hits]== 5 or lfpmt[new_hits]== 6 or lfpmt[new_hits]== 7 or lfpmt[new_hits]== 12 or lfpmt[new_hits]== 13 or lfpmt[new_hits]== 14 or lfpmt[new_hits]== 15){
           Y[numberofhitsin_y]= lfypos[new_hits];
           ZY[numberofhitsin_y]=  lfzpos[new_hits];
           errorY[numberofhitsin_y] = 2.64/(sqrt(12.0));
           errorZY[numberofhitsin_y] = 1.06/(sqrt(12.0));
           numberofhitsin_y++;
         }
       }
     }

}
}
     /* I prepare this for you : If you would like to save the plots in a root file,
     please check :)

     std::string coord_namesimu =  "Hits_Proj_"+ std::to_string(n)+"_"+std::to_string(NChannels);
     TGraphErrors *grx_point;
     TGraphErrors *gry_point;
     TCanvas* c1 = new TCanvas("c1",coord_namesimu.c_str());
     c1->SetCanvasSize(800, 800);
     c1->SetWindowSize(800, 800);
     c1->Divide(1,2);
     c1->SetName(coord_namesimu.c_str());
     grx_point = new TGraphErrors (nhit_left, ZX, X,errorZX, errorX);
     std::string name_x = "XZ : Evt "+ std::to_string(n)+"_"+std::to_string(ni);
     grx_point->SetName(name_x.c_str());
     grx_point->SetTitle(name_x.c_str());
     grx_point->GetXaxis()->SetTitle("X(cm)");
     grx_point->GetYaxis()->SetTitle("Z(cm)");
     grx_point->GetXaxis()->SetLimits(-200,200);
     grx_point->GetYaxis()->SetRangeUser(-3000, 3000);
     grx_point->SetLineColor(kRed);
     grx_point->SetMarkerColor(kRed);
     grx_point->SetMarkerSize(1.0);
     grx_point->SetMarkerStyle(20);

     gry_point = new TGraphErrors (nhit_right, Y, ZY,errorZY, errorY);
     std::string name_y = "YZ : Evt "+ std::to_string(n)+"_"+std::to_string(nj);
     gry_point->SetName(name_y.c_str());
     gry_point->SetTitle(name_y.c_str());
     gry_point->GetXaxis()->SetTitle("Y(cm)");
     gry_point->GetYaxis()->SetTitle("Z(cm)");
     gry_point->GetYaxis()->SetLimits(-200,200);
     gry_point->GetXaxis()->SetLimits(-1500, 1500);
     gry_point->SetLineColor(kBlue);
     gry_point->SetMarkerColor(kBlue);
     gry_point->SetMarkerSize(1.0);
     gry_point->SetMarkerStyle(20);

     c1->SetGridx();
     c1->SetGridy();
     c1->cd(1);
     gPad->SetGrid();
     grx_point->Draw("AP");

     c1->SetGridx();
     c1->SetGridy();
     c1->cd(1);
     gPad->SetGrid();
     gry_point->Draw("AP");

     c1->SaveAs("output_root_file.root"); //Here you need to save the ouput root file
     c1->Clear();*/
