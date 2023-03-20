
void rnd()
{
  TRandom2  *rand = new TRandom2();
  
  for(int i = 0; i < 5; i++)
  {  
    int r = rand->Integer(64000);
    cout<<r<<endl;

  } 
}
