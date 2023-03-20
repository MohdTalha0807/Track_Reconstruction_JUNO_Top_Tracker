#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <map>

using namespace std;

int main()
{
    map<int, vector<float>> mymap; 
    vector<float> v; 
    int tdm_ch, tdm;
    float cx, cy, cz;
    
   fstream myFile;
   myFile.open("ChannelPosition.txt", ios::in);
   
   if(myFile.is_open())
   {
     
     
     while( myFile>> tdm_ch >> cx >> cy >> cz)
     {
           v.push_back(cx);         
           v.push_back(cy);
           v.push_back(cz);
           mymap[tdm_ch]= v;
           tdm = tdm_ch; 
     }
     
     myFile.close();
     
   }

      vector<float> tmp = mymap[tdm]; 
       
      cout<<tmp[0]<<"\t"<<tmp[1]<<"\t"<<tmp[2]<<endl; 
 }


