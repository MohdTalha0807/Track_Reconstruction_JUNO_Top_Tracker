//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  2 22:46:04 2022 by ROOT version 6.24/02
// from TTree TT_Tracks/TT_Tracks
// found on file: output_canvas.root
//////////////////////////////////////////////////////////

#ifndef TT_Tracks_h
#define TT_Tracks_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using namespace std;

class TT_Tracks {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Num_TracksXZ_geninfo;
   Int_t           Num_TracksYZ_geninfo;
   vector<float>   *ThetaRhoXZ_gen;
   vector<float>   *ThetaRhoYZ_gen;
   Int_t           Num_TracksXZ_hough;
   Int_t           Num_TracksYZ_hough;
   vector<float>   *ThetaRhoXZ_hough;
   vector<float>   *ThetaRhoYZ_hough;

   // List of branches
   TBranch        *b_ntracksXZ_gen;   //!
   TBranch        *b_ntracksYZ_gen;   //!
   TBranch        *b_ThetaRhoXZ_gen;   //!
   TBranch        *b_ThetaRhoYZ_gen;   //!
   TBranch        *b_ntracksXZ_hough;   //!
   TBranch        *b_ntracksYZ_hough;   //!
   TBranch        *b_ThetaRhoXZ_hough;   //!
   TBranch        *b_ThetaRhoYZ_hough;   //!

   TT_Tracks(TTree *tree=0);
   virtual ~TT_Tracks();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TT_Tracks_cxx
TT_Tracks::TT_Tracks(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_canvas.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_canvas.root");
      }
      f->GetObject("TT_Tracks",tree);

   }
   Init(tree);
}

TT_Tracks::~TT_Tracks()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TT_Tracks::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TT_Tracks::LoadTree(Long64_t entry)
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

void TT_Tracks::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ThetaRhoXZ_gen = 0;
   ThetaRhoYZ_gen = 0;
   ThetaRhoXZ_hough = 0;
   ThetaRhoYZ_hough = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Num_TracksXZ_geninfo", &Num_TracksXZ_geninfo, &b_ntracksXZ_gen);
   fChain->SetBranchAddress("Num_TracksYZ_geninfo", &Num_TracksYZ_geninfo, &b_ntracksYZ_gen);
   fChain->SetBranchAddress("ThetaRhoXZ_gen", &ThetaRhoXZ_gen, &b_ThetaRhoXZ_gen);
   fChain->SetBranchAddress("ThetaRhoYZ_gen", &ThetaRhoYZ_gen, &b_ThetaRhoYZ_gen);
   fChain->SetBranchAddress("Num_TracksXZ_hough", &Num_TracksXZ_hough, &b_ntracksXZ_hough);
   fChain->SetBranchAddress("Num_TracksYZ_hough", &Num_TracksYZ_hough, &b_ntracksYZ_hough);
   fChain->SetBranchAddress("ThetaRhoXZ_hough", &ThetaRhoXZ_hough, &b_ThetaRhoXZ_hough);
   fChain->SetBranchAddress("ThetaRhoYZ_hough", &ThetaRhoYZ_hough, &b_ThetaRhoYZ_hough);
   Notify();
}

Bool_t TT_Tracks::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TT_Tracks::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TT_Tracks::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TT_Tracks_cxx
