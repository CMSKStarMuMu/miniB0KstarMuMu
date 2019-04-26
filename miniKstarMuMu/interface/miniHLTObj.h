#ifndef  miniHLTObj_h
#define  miniHLTObj_h

#include "TROOT.h"
#include "TMath.h"
#include <vector>
#include <string>


class miniHLTObj {
public:

  std::string pathName; // name of filter passed by the object
  Float_t pt;            // pt of the object passing the filter [GeV]
  Float_t eta;           // eta of the object passing the filter
  Float_t phi;           // phi of the object passing the filter
  Int_t   charge;       
  Int_t   pdgId;        
  
  miniHLTObj(){};
  virtual ~miniHLTObj(){};

//   ClassDef(miniHLTObj,1)

};

// typedef std::vector<miniHLTObj> miniHLTObjCollection;

#endif