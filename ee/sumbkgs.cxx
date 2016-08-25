#include <iostream>
#include <iomanip>
#include <cmath>

double total=0;
double toterr=0;

void printshit(char *source){
  ifstream infile;
  infile.open(source);
  double n=0;	      
  double Nexp=0;
  double err=0;
  double temp=0;
  while(infile>>n){
    //cout<<n<<endl;
    Nexp+=n;
    infile>>n;
    //cout<<n<<endl;
    err+=(n*n);
    //cout<<endl;
  }
  err=sqrt(err);
  
  cout<<endl;
  cout<<fixed; 
  cout<<setprecision(1);
  cout<<setw(16)<<source<<setw(9)<<" Nexp = "<<setw(7)<<Nexp<<" +/- "<<err<<endl;
  
  total+=Nexp;
  toterr+=(err*err);
  
  infile.close();
}

void sumbkgs(){
  printshit("foils.zr96.dat");
  printshit("ca48.dat");
  printshit("nd150.dat");
  printshit("exbg.dat");
  //printshit("radon.dat");
  cout<<endl;
  cout<<fixed; 
  cout<<setprecision(1);
  cout<<setw(16)<<" TOTAL "<<setw(9)<<" Nexp = "<<setw(7)<<total<<" +/- "<<sqrt(toterr)<<endl;
  cout<<endl;
  
}

