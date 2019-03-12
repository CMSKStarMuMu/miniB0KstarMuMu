// #define G__DICTIONARY

#include "DataFormats/Common/interface/Wrapper.h"
#include "miniB0KstarMuMu/miniKstarMuMu/interface/miniHLTObj.h"

// namespace { struct dictionary {
//   miniHLTObjCollection dummy0;
//   edm::Wrapper<miniHLTObjCollection> dummy1;
//   };
// }
// 
namespace {
  struct miniB0KstarMuMu_miniKstarMuMu {
    miniHLTObj miniHLTObj_;
    std::vector<miniHLTObj> vmh;
    edm::Wrapper<std::vector<miniHLTObj>> wvmh;
  };
}
