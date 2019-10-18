//  Example demostrating calculation of CLs for the on/off problem.
//  Glen Cowan, RHUL Physics, July 2013.
//  Based on part 6 of roostats macro HybridInstructional.C
//  To run, at the root prompt type: .x SimpleCLs.C

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooProdPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TH1.h>
#include <RooPlot.h>
#include <RooMsgService.h>

#include <RooStats/NumberCountingUtils.h>
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

void SimpleCLs();

int main() {

  SimpleCLs();
  return 0;

}

void SimpleCLs() {

  TCanvas* c1 = new TCanvas;

  // Define the model:  n ~ Poisson (s+b) and m ~ Poisson (tau*b), set tau = 1.
  RooWorkspace* w = new RooWorkspace("w");
  w->factory("Poisson::Pn(n[150,0,500], sum::splusb(s[0,0,100], b[100,0,300]))");
  w->factory("Poisson::Pm(m[100,0,500], prod::taub(tau[1.], b))");
  w->factory("PROD::model(Pn,Pm)");
  w->defineSet("poi", "s");

  // create a toy dataset with the observed data values n and m
  double n = 120.;
  double m = 100.;
  w->defineSet("obsNM","n,m");   
  w->var("n")->setVal(n);
  w->var("m")->setVal(m);
  RooDataSet* dataNM = new RooDataSet("dNM", "dNM", *w->set("obsNM"));
  dataNM->add(*w->set("obsNM"));

  // now we need new model configs, with PDF="model"
  ModelConfig b_modelNM("B_modelNM", w);
  b_modelNM.SetPdf(*w->pdf("model"));
  b_modelNM.SetObservables(*w->set("obsNM"));
  b_modelNM.SetParametersOfInterest(*w->set("poi"));
  w->var("s")->setVal(0.0);
  b_modelNM.SetSnapshot(*w->set("poi"));     // sets up b hypothesis as s = 0

  // create the alternate (s+b) ModelConfig with given value of s
  double s = 40.;
  ModelConfig sb_modelNM("S+B_modelNM", w);
  sb_modelNM.SetPdf(*w->pdf("model"));
  sb_modelNM.SetObservables(*w->set("obsNM"));
  sb_modelNM.SetParametersOfInterest(*w->set("poi"));
  RooRealVar* poi = (RooRealVar*) sb_modelNM.GetParametersOfInterest()->first();
  w->var("s")->setVal(s);
  sb_modelNM.SetSnapshot(*w->set("poi"));  // set up sb hypothesis with given s

  // test statistic \lambda(s) = -log L(s,\hat\hat{b})/L(\hat{s},\hat{b})
  ProfileLikelihoodTestStat profll(*sb_modelNM.GetPdf());

  // Set up Hybrid Calculator; b is alternate, sb is null
  HybridCalculator hc(*dataNM, b_modelNM, sb_modelNM);
  ToyMCSampler* toymcs = (ToyMCSampler*)hc.GetTestStatSampler();
  toymcs->SetNEventsPerToy(1);
  hc.SetToys(2000, 2000);                // # of toy exp for each hypothesis
  toymcs->SetTestStatistic(&profll);

  // Get the result and compute CLs
  HypoTestResult* result = hc.GetHypoTest();
  result->SetPValueIsRightTail(true);
  double psb = result->NullPValue();
  result->SetPValueIsRightTail(false);       // this affects pb
  cout << "results based on toy MC:" << endl;
  double pb = result->AlternatePValue();
  cout << "psb = " << psb << endl;
  cout << "pb  = " << pb << endl;
  double clb = 1. - pb;
  double clsb = psb;
  double cls = 9999.;
  if ( clb > 0 ) { cls = clsb/clb; }
  cout << "cls = " << cls << endl;
  cout << endl;

  // Make a plot
  c1->SetLogy();
  result->SetPValueIsRightTail(true);
  HypoTestPlot* plot = new HypoTestPlot(*result, 80);
  plot->Draw();
  c1->SaveAs("SimpleCLs.pdf");

  // Now compute using asymptotic formula; b is alt, sb is null
  AsymptoticCalculator ac(*dataNM, b_modelNM, sb_modelNM);
  ac.SetOneSided(false);     // KLUDGE -- should want one sided (true) for limit
  AsymptoticCalculator::SetPrintLevel(-1);
  HypoTestResult* asympResult = ac.GetHypoTest();

  asympResult->SetPValueIsRightTail(false);         // appears to do nothing?!
  double asymp_pb = 1. - asympResult->AlternatePValue(); // KLUDGE!!! Needs 1 -   
  asympResult->SetPValueIsRightTail(true);
  double asymp_psb = asympResult->NullPValue();

  cout << "Results based on asymptotic formulae:" << endl;
  cout << "psb = " << asymp_psb << endl;
  cout << "pb  = " << asymp_pb << endl;
  double asymp_clb = 1. - asymp_pb;
  double asymp_clsb = asymp_psb;
  double asymp_cls = 9999.;
  if ( asymp_clb > 0 ) {
    asymp_cls = asymp_clsb/asymp_clb;
  }
  cout << "cls = " << asymp_cls << endl;
  cout << endl;

  // create hypotest inverter passing the desired calculator (hc or ac)
  HypoTestInverter calc(ac);
  calc.SetVerbose(false);
  calc.SetConfidenceLevel(0.95);
  bool useCLs = true;
  calc.UseCLs(useCLs);
  if (useCLs) { profll.SetOneSided(true); }

  int npoints = 10;  // number of points to scan
  // min and max for scan (better to choose smaller intervals)
  double poimin = poi->getMin();
  double poimax = poi->getMax();
  cout << "Doing a fixed scan  in interval : " 
       << poimin << " , " << poimax << endl;
  calc.SetFixedScan(npoints, poimin, poimax);
  HypoTestInverterResult* r = calc.GetInterval();

  double upperLimit = r->UpperLimit();
  double ulError = r->UpperLimitEstimatedError();
  cout << "The computed upper limit is: " << upperLimit 
       << " +/- " << ulError << endl;
  
   // compute expected limit
   cout << "Expected upper limits using b-only model : " << endl;
   cout << "median limit = " << r->GetExpectedUpperLimit(0) << endl;
   cout << "med. limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
   cout << "med. limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
   cout << endl;

   // plot result of the scan 
   HypoTestInverterPlot* plot2 = 
      new HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
   TCanvas* c2 = new TCanvas("HypoTestInverter Scan"); 
   c2->SetLogy(false);
   plot2->Draw("2CL");
   c2->SaveAs("SimpleCLsLimit.pdf");

}