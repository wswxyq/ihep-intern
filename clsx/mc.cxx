//  Example demostrating calculation of CLs for the on/off problem.
//  Glen Cowan, RHUL Physics, July 2013.
//  Shaowei Wu, modified on Oct 20, 2019
//  Based on part 6 of roostats macro HybridInstructional.C
//  To run, at the root prompt type: .x ws2file.C

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooPolynomial.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooGenericPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TH1.h>
#include <RooPlot.h>
#include <RooMsgService.h>

#include <TChain.h>
#include <TFile.h>

#include <RooStats/NumberCountingUtils.h>
#include <RooStats/ModelConfig.h>
#include <RooStats/HybridCalculator.h>
#include <RooStats/AsymptoticCalculator.h>
#include <RooStats/ToyMCSampler.h>
#include <RooStats/HypoTestPlot.h>
#include <RooStats/HypoTestInverter.h>
#include <RooStats/HypoTestInverterPlot.h>
#include <RooStats/HypoTestResult.h>
#include <RooStats/ProfileLikelihoodTestStat.h>
#include <RooStats/SimpleLikelihoodRatioTestStat.h>
#include <RooStats/RatioOfProfiledLikelihoodsTestStat.h>
#include <RooStats/MaxLikelihoodEstimateTestStat.h>

using namespace RooFit;
using namespace RooStats;
using namespace std;


void mc() {

}

