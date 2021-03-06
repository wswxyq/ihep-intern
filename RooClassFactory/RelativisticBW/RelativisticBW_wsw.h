/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef RELATIVISTICBW_WSW
#define RELATIVISTICBW_WSW

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RelativisticBW_wsw : public RooAbsPdf {
public:
  RelativisticBW_wsw() {} ; 
  RelativisticBW_wsw(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _M,
	      RooAbsReal& _gamma);
  RelativisticBW_wsw(const RelativisticBW_wsw& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RelativisticBW_wsw(*this,newname); }
  inline virtual ~RelativisticBW_wsw() { }

protected:

  RooRealProxy x ;
  RooRealProxy M ;
  RooRealProxy gamma ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RelativisticBW_wsw,1) // Your description goes here...
};
 
#endif
