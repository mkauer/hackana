#include "h10.h"
#include "TH2.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <string>

Int_t h10::Cut(Int_t * results){
  //User supplied function for slim2 program
  //Called from h10::Loop for each event
  //returns 1 if event selected, 0 otherwise
  //For the cuts passed by this events, increment results[20] array (if interested)


  //Cut values for event preselection;

  //Selection code
   
   return 1;

}
