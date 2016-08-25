#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class cluster
{ public:
  Int_t ngroup;
  Int_t index[250];
  Float_t Esum;
  Int_t Iseed;
  Float_t probint;
  Float_t probext;
  Bool_t assigned;
};

class scintillator
{ public:
  Int_t icluster;
  cluster ncluster[250];
};
