/*

Compile with:
g++ rootread.cpp testTree.C -o rootread `root-config --cflags --libs`

*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstring>

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

Double_t Hough_rho(Double_t x, Double_t y, Double_t theta)
{
  Double_t pi= (TMath::ATan(1))*4;
  Double_t theta_degree = theta;
  Double_t theta_rad =(pi*theta_degree)/(180);
  Double_t x0 = x;
  Double_t y0 = y;
  Double_t rho = (x0*TMath::Cos(theta_rad)) + (y0*(TMath::Sin(theta_rad)));
  return rho;
}


int main() {

  TFile *fileIn1 = new TFile("muon_265_user.root","READ");
  TFile *fileIn2 = new TFile("muon_265_user.root","READ");
  TTree *treeIn1 = (TTree*)fileIn1->Get("TTDigit");  // creating a pointer of the current TTree
  TTree *treeIn2 = (TTree*)fileIn2->Get("geninfo"); 

































      
