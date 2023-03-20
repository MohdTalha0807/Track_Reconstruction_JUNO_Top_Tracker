//here write the libraries
#include "C_Libraries.h"
#include "Root_Libraries.h"
#include "Telescope_utils.h"
#include "Hough_utils.h"
#include "Telescope_Modules.h"
#include "compa.h"
#include <string>
using namespace std;


int main(int argc, char ** argv){ //Open main
 
  if(argc!=5){
    std::cout<<"Insufficient no. of arguments detected"<<std::endl;
    std::cout<<"Follow this:"<<std::endl;
    std::cout<<"./HOUGH_TOP_TRACKER lim1 lim2 coinc.root outfile.root "<<std::endl;
    return -1;
  }

  //Inputs and Outputs
  int lim1 = std::stoi(argv[1]) ;         //Interval of entries first entries //example 0
  int lim2 = std::stoi(argv[2]) ;         //Interval of entries last entries  // 1000
  string infile1    = argv[3];            // Data.root (simulation file )
  string outfile   = argv[4];             //Output file , please write just the name without extension
 

  TFile *f = new TFile(infile1.c_str());   //Simulation open
  Bool_t stat= f->IsOpen();

  if(stat){  //Open the Simulation file
    cout << "File is now open!\nContains:\n";
    f->ls();

  TTree *tree_TTDigit=(TTree*)f->Get("TTDigit"); //Read TTDigit : Simulation file
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

    string root_ext = ".root";
    string filename_root = outfile + root_ext;
    TFile * Outfile          = new TFile(filename_root.c_str(),"RECREATE"); //Define your root file (output)
    Outfile->cd();

    //read event by event on the TTDigit

    for(int n=lim1; n<lim2;n++){//Open Nentries  : numbers of events
      tree_TTDigit->GetEntry(n);
      if(NChannels==0)continue;
      Int_t lflayers=0, lfpmt=0, lfstrip=0 , lfchannel=0 ,lfrow=0, lfcolumn=0;
    
      cout<<"number of event : "<<n<<endl;
      cout<<"---------------------------------------------------"<<endl;
 
      for(int ihit=0; ihit<NChannels; ihit++){//Open NChannels  : numbers of hits in one event
          //--------- Read the JUNOID channel for each hit ----------
          lfchannel = Channel_ID[ihit];
          lflayers = LayerFromJUNOid(lfchannel);
          lfpmt = PMTsideFromJUNOid(lfchannel);
          lfstrip = StripFromJUNOid(lfchannel);
          lfrow=RowFromJUNOid(lfchannel);
          lfcolumn=ColumnFromJUNOid(lfchannel);

          cout<<"layer:"<<lflayers <<endl;
          cout<<"pmt:"<< lfpmt  <<endl;
          cout<<"column:"<< lfcolumn <<endl;
          cout<<"strip:"<< lfstrip <<" : "<<position_xaxis(lfpmt, lfcolumn, lfstrip)<< endl;
          cout<<"---------------------------------------------------"<<endl;
        }
     }

   }
      cout<<endl;
      cout<<" Finnished loop entries "<<endl;
      Outfile->Close();

   }

  return 0;
  }//Close main

