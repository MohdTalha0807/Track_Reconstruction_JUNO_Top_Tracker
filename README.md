# Track_Reconstruction_JUNO_Top_Tracker

# Overview
We worked on the Development of a reconstruction algorithm for the
JUNO Top Tracker which is part of the veto system along with the Water Cherenkov
Detector to reject backgrounds coming from the muon interactions in the central detector.
Therefore, knowing precisely where the muons are coming from is essential
for the veto system to be efficient. Now the current track reconstruction algorithm is
based on χ2 fit of the data which becomes inefficient as the noise rate increases in the
detector. The master project involved the development of an efficient algorithm to implement
the Hough Transform method to perform a selection of the muon PMT triggers
as compared to those triggered by the noise (thanks to their random nature) and thus
using the topological difference between them for the track reconstruction. 

# Data Structure
The folder TT_tracks_Recon contains the source file (TT_Track_Recon.cc) along with all the include, Makefile and the input files. The program takes the input from the generated muon simulation root file (muon_265_user.root) by the JUNO official software.

# Setting
Prerequisites: Root, C++ (better if the latest version)
