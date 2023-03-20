//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar  7 10:45:24 2022 by ROOT version 6.24/02
// from TTree TTDigit/PE TT
// found on file: muon_265_user.root
//////////////////////////////////////////////////////////

#ifndef TTDigit_h
#define TTDigit_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class TTDigit {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evtID;
   Int_t           NTouchedChannel;
   Int_t           TB_channelC[92];   //[NTouchedChannel]
   Int_t           TB_DMchannel[92];   //[NTouchedChannel]
   Float_t         TB_pe[92];   //[NTouchedChannel]
   Float_t         TB_time[92];   //[NTouchedChannel]
   Int_t           TB_ADC[92];   //[NTouchedChannel]
   Float_t         TB_xcC[92];   //[NTouchedChannel]
   Float_t         TB_ycC[92];   //[NTouchedChannel]
   Float_t         TB_zcC[92];   //[NTouchedChannel]
   Int_t           TB_is_ctC[92];   //[NTouchedChannel]
   Int_t           TB_isMuonDepositsC[92];   //[NTouchedChannel]

   // List of branches
   TBranch        *b_evtID;   //!
   TBranch        *b_NTouchedChannel;   //!
   TBranch        *b_TB_channelC;   //!
   TBranch        *b_TB_DMchannel;   //!
   TBranch        *b_TB_pe;   //!
   TBranch        *b_TB_time;   //!
   TBranch        *b_TB_ADC;   //!
   TBranch        *b_TB_xcC;   //!
   TBranch        *b_TB_ycC;   //!
   TBranch        *b_TB_zcC;   //!
   TBranch        *b_TB_is_ctC;   //!
   TBranch        *b_TB_isMuonDepositsC;   //!

   TTDigit(TTree *tree=0);
   virtual ~TTDigit();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TTDigit_cxx
TTDigit::TTDigit(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("muon_265_user.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("muon_265_user.root");
      }
      f->GetObject("TTDigit",tree);

   }
   Init(tree);
}

TTDigit::~TTDigit()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TTDigit::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TTDigit::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TTDigit::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evtID", &evtID, &b_evtID);
   fChain->SetBranchAddress("NTouchedChannel", &NTouchedChannel, &b_NTouchedChannel);
   fChain->SetBranchAddress("TB_channelC", TB_channelC, &b_TB_channelC);
   fChain->SetBranchAddress("TB_DMchannel", TB_DMchannel, &b_TB_DMchannel);
   fChain->SetBranchAddress("TB_pe", TB_pe, &b_TB_pe);
   fChain->SetBranchAddress("TB_time", TB_time, &b_TB_time);
   fChain->SetBranchAddress("TB_ADC", TB_ADC, &b_TB_ADC);
   fChain->SetBranchAddress("TB_xcC", TB_xcC, &b_TB_xcC);
   fChain->SetBranchAddress("TB_ycC", TB_ycC, &b_TB_ycC);
   fChain->SetBranchAddress("TB_zcC", TB_zcC, &b_TB_zcC);
   fChain->SetBranchAddress("TB_is_ctC", TB_is_ctC, &b_TB_is_ctC);
   fChain->SetBranchAddress("TB_isMuonDepositsC", TB_isMuonDepositsC, &b_TB_isMuonDepositsC);
   Notify();
}

Bool_t TTDigit::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TTDigit::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TTDigit::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TTDigit_cxx
