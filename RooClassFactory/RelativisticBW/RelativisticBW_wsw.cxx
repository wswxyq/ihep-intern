/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "RelativisticBW_wsw.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"
using namespace TMath;

ClassImp(RelativisticBW_wsw);

RelativisticBW_wsw::RelativisticBW_wsw(const char *name, const char *title,
                                       RooAbsReal& _x,
                                       RooAbsReal& _M,
                                       RooAbsReal& _gamma) :
    RooAbsPdf(name,title),
    x("x","x",this,_x),
    M("M","M",this,_M),
    gamma("gamma","gamma",this,_gamma)
{
}


RelativisticBW_wsw::RelativisticBW_wsw(const RelativisticBW_wsw& other, const char* name) :
    RooAbsPdf(other,name),
    x("x",this,other.x),
    M("M",this,other.M),
    gamma("gamma",this,other.gamma)
{
}



Double_t RelativisticBW_wsw::evaluate() const
{
    // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
    return 14.0*gamma*gamma*M*M/(22.0* (x*x-M*M)*(x*x-M*M)+x*x*x*x*(gamma*gamma)/(M*M) );
}



