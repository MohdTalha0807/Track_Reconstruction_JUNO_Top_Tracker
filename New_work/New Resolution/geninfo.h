//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 25 09:41:25 2022 by ROOT version 6.18/00
// from TTree geninfo/Generator Info
// found on file: muon_265_user.root
//////////////////////////////////////////////////////////

#ifndef geninfo_h
#define geninfo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class geninfo {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evtID;
   Int_t           nInitParticles;
   Int_t           InitPDGID[6];   //[nInitParticles]
   Int_t           InitTRKID[6];   //[nInitParticles]
   Float_t         InitX[6];   //[nInitParticles]
   Float_t         InitY[6];   //[nInitParticles]
   Float_t         InitZ[6];   //[nInitParticles]
   Float_t         InitPX[6];   //[nInitParticles]
   Float_t         InitPY[6];   //[nInitParticles]
   Float_t         InitPZ[6];   //[nInitParticles]
   Float_t         InitMass[6];   //[nInitParticles]
   Double_t        InitTime[6];   //[nInitParticles]
   Float_t         ExitX[6];   //[nInitParticles]
   Float_t         ExitY[6];   //[nInitParticles]
   Float_t         ExitZ[6];   //[nInitParticles]
   Double_t        ExitT[6];   //[nInitParticles]
   Float_t         ExitPX[6];   //[nInitParticles]
   Float_t         ExitPY[6];   //[nInitParticles]
   Float_t         ExitPZ[6];   //[nInitParticles]
   Float_t         TrackLength[6];   //[nInitParticles]

   // List of branches
   TBranch        *b_evtID;   //!
   TBranch        *b_nInitParticles;   //!
   TBranch        *b_InitPDGID;   //!
   TBranch        *b_InitTRKID;   //!
   TBranch        *b_InitX;   //!
   TBranch        *b_InitY;   //!
   TBranch        *b_InitZ;   //!
   TBranch        *b_InitPX;   //!
   TBranch        *b_InitPY;   //!
   TBranch        *b_InitPZ;   //!
   TBranch        *b_InitMass;   //!
   TBranch        *b_InitTime;   //!
   TBranch        *b_ExitX;   //!
   TBranch        *b_ExitY;   //!
   TBranch        *b_ExitZ;   //!
   TBranch        *b_ExitT;   //!
   TBranch        *b_ExitPX;   //!
   TBranch        *b_ExitPY;   //!
   TBranch        *b_ExitPZ;   //!
   TBranch        *b_TrackLength;   //!

   geninfo(TTree *tree=0);
   virtual ~geninfo();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef geninfo_cxx
geninfo::geninfo(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("muon_265_user.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("muon_265_user.root");
      }
      f->GetObject("geninfo",tree);

   }
   Init(tree);
}

geninfo::~geninfo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t geninfo::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t geninfo::LoadTree(Long64_t entry)
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

void geninfo::Init(TTree *tree)
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
   fChain->SetBranchAddress("nInitParticles", &nInitParticles, &b_nInitParticles);
   fChain->SetBranchAddress("InitPDGID", InitPDGID, &b_InitPDGID);
   fChain->SetBranchAddress("InitTRKID", InitTRKID, &b_InitTRKID);
   fChain->SetBranchAddress("InitX", InitX, &b_InitX);
   fChain->SetBranchAddress("InitY", InitY, &b_InitY);
   fChain->SetBranchAddress("InitZ", InitZ, &b_InitZ);
   fChain->SetBranchAddress("InitPX", InitPX, &b_InitPX);
   fChain->SetBranchAddress("InitPY", InitPY, &b_InitPY);
   fChain->SetBranchAddress("InitPZ", InitPZ, &b_InitPZ);
   fChain->SetBranchAddress("InitMass", InitMass, &b_InitMass);
   fChain->SetBranchAddress("InitTime", InitTime, &b_InitTime);
   fChain->SetBranchAddress("ExitX", ExitX, &b_ExitX);
   fChain->SetBranchAddress("ExitY", ExitY, &b_ExitY);
   fChain->SetBranchAddress("ExitZ", ExitZ, &b_ExitZ);
   fChain->SetBranchAddress("ExitT", ExitT, &b_ExitT);
   fChain->SetBranchAddress("ExitPX", ExitPX, &b_ExitPX);
   fChain->SetBranchAddress("ExitPY", ExitPY, &b_ExitPY);
   fChain->SetBranchAddress("ExitPZ", ExitPZ, &b_ExitPZ);
   fChain->SetBranchAddress("TrackLength", TrackLength, &b_TrackLength);
   Notify();
}

Bool_t geninfo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void geninfo::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t geninfo::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef geninfo_cxx
