//  Example demostrating calculation of CLs for the on/off problem.
//  Glen Cowan, RHUL Physics, July 2013.
//  Based on part 6 of roostats macro HybridInstructional.C
//  To run, at the root prompt type: .x myCLs.C

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooPolynomial.h>
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

void myCLs();

int main() {

    myCLs();
    return 0;

}

void myCLs() {

    TCanvas* c1 = new TCanvas;

    // ------------------Define the model:  n ~ Poisson (s+b) and m ~ Poisson (tau*b), set tau = 1.-------------------------
    RooWorkspace* w = new RooWorkspace("w");
    //w->factory("Poisson::Pn(n[150,0,500], sum::splusb(s[0,0,100], b[100,0,300]))");
    //w->factory("Poisson::Pm(m[100,0,500], prod::taub(tau[1.], b))");
    //w->factory("PROD::model(Pn,Pm)");
    //w->defineSet("poi", "s");
    RooRealVar B_DTF_M("B_DTF_M", "B_DTF_M", 4200, 4600, "MeV");
    RooGenericPdf rtbw4312("rtbw4312", "14.0*5.3*5.3*4312.0*4312.0/(22.0*(B_DTF_M*B_DTF_M-4312.0*4312.0)*(B_DTF_M*B_DTF_M-4312.0*4312.0)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(5.3*5.3)/(4312.0*4312.0))", RooArgList(B_DTF_M));
    RooGenericPdf rtbw4440("rtbw4400", "14.0*25.2*25.2*4440.2*4440.2/(22.0*(B_DTF_M*B_DTF_M-4440.2*4440.2)*(B_DTF_M*B_DTF_M-4440.2*4440.2)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(25.2*25.2)/(4440.2*4440.2))", RooArgList(B_DTF_M));
    RooGenericPdf rtbw4457("rtbw4457", "14.0*5.5*5.5*4456.6*4456.6/(22.0*(B_DTF_M*B_DTF_M-4456.6*4456.6)*(B_DTF_M*B_DTF_M-4456.6*4456.6)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(5.5*5.5)/(4456.6*4456.6))", RooArgList(B_DTF_M));    
    RooGenericPdf rtbwx("rtbwx", "14.0*62.7*62.7*4394.7*4394.7/(22.0*(B_DTF_M*B_DTF_M-4394.7*4394.7)*(B_DTF_M*B_DTF_M-4394.7*4394.7)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(62.7*62.7)/(4394.7*4394.7))", RooArgList(B_DTF_M));
    //fraction of b-w
    RooRealVar signal_frac_4312("signal_frac_4312", "signal_frac_4312", 0.1, 0., 1.);
    RooRealVar signal_frac_4440("signal_frac_4440", "signal_frac_4440", 0.1, 0., 1.);
    RooRealVar signal_frac_4457("signal_frac_4457", "signal_frac_4457", 0.1, 0., 1.);
    RooRealVar signal_frac("signal_frac", "signal_frac", 0.0001, 0., 0.0002);
    //polynimial parameter
    //fit parameter first//////////////////////////
    RooRealVar x1("x1", "para1", 74.3657);
	RooRealVar x2("x2", "para2", -6.47583e-05);
	RooRealVar x3("x3", "para3", -3.33708e-06);

    RooAddPdf smodel("smodel", "smodel", RooArgList(rtbw4312, rtbw4440, rtbw4457, rtbwx),
                        RooArgList(signal_frac_4312, signal_frac_4440, signal_frac_4457));
    RooPolynomial bmodel("bmodel", "bmodel", B_DTF_M, RooArgList(x1, x2, x3));

    RooAddPdf model("model", "model", RooArgList(smodel, bmodel), RooArgList(signal_frac));

    //RooRealVar s("s", "s", 0, 0, 4658);
    //RooRealVar b("b", "b", 0, 4658, 4658);

    w->import(model);
    //w->import(s);
    w->defineSet("poi", "signal_frac");
    w->Print();


    
    cout<<"========== create a toy dataset with the observed data values n and m =========="<<endl;
    //double n = 120.;
    //double m = 100.;
    //w->defineSet("obsNM","n,m");   
    //w->var("n")->setVal(n);
    //w->var("m")->setVal(m);
    //RooDataSet* dataNM = new RooDataSet("dNM", "dNM", *w->set("obsNM"));
    //dataNM->add(*w->set("obsNM"));

    //set chain
    TChain *chain = new TChain("ReducedTree");
	//load file
    chain->Add("../data_selection/BDT_reduced.root");
    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(B_DTF_M), Import(*chain), 
		            Cut(""));
    w->defineSet("B_DTF_M", "B_DTF_M");
    w->import(*ds);
    w->Print();


    
    cout<<"========== now we need new model configs, with PDF=\"model\"=========="<<endl;
    ModelConfig b_modelNM("B_modelNM", w);
    b_modelNM.SetPdf(*w->pdf("model"));
    b_modelNM.SetObservables(*w->set("B_DTF_M"));
    b_modelNM.SetParametersOfInterest(*w->set("poi"));
    w->var("signal_frac")->setVal(0.0);
    b_modelNM.SetSnapshot(*w->set("poi"));     // sets up b hypothesis as s = 0

    
    cout<<"========== create the alternate (s+b) ModelConfig with given value of s =========="<<endl;
    double s_value = 0.001;
    ModelConfig sb_modelNM("S+B_modelNM", w);
    sb_modelNM.SetPdf(*w->pdf("model"));
    sb_modelNM.SetObservables(*w->set("B_DTF_M"));
    sb_modelNM.SetParametersOfInterest(*w->set("poi"));
    RooRealVar* poi = (RooRealVar*) sb_modelNM.GetParametersOfInterest()->first();
    w->var("signal_frac")->setVal(s_value);
    sb_modelNM.SetSnapshot(*w->set("poi"));  // set up sb hypothesis with given s

    
    cout<<"==========  test statistic \\lambda(s) = -log L(s,\\hat\\hat{b})/L(\\hat{s},\\hat{b})=========="<<endl;
    ProfileLikelihoodTestStat profll(*sb_modelNM.GetPdf());
    
    cout<<"========== Set up Hybrid Calculator; b is alternate, sb is null =========="<<endl;
    HybridCalculator hc(*ds, b_modelNM, sb_modelNM);
    ToyMCSampler* toymcs = (ToyMCSampler*)hc.GetTestStatSampler();
    toymcs->SetNEventsPerToy(1);
    hc.SetToys(2000, 2000);                // # of toy exp for each hypothesis
    toymcs->SetTestStatistic(&profll);

    b_modelNM.Print();
    sb_modelNM.Print();

    
    cout<<"========== Get the result and compute CLs =========="<<endl;
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
    
    

    cout<<"========== Make a plot =========="<<endl;
    c1->SetLogy();
    result->SetPValueIsRightTail(true);
    HypoTestPlot* plot = new HypoTestPlot(*result, 80);
    plot->Draw();
    c1->SaveAs("myCLs.pdf");
    
    // Now compute using asymptotic formula; b is alt, sb is null
    /*
    AsymptoticCalculator ac(*ds, b_modelNM, sb_modelNM);
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
    */

    
    cout<<"========== create hypotest inverter passing the desired calculator (hc or ac) =========="<<endl;
    HypoTestInverter calc(hc);
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

    cout<<"========== compute expected limit =========="<<endl;
    cout << "Expected upper limits using b-only model : " << endl;
    cout << "median limit = " << r->GetExpectedUpperLimit(0) << endl;
    cout << "med. limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << endl;
    cout << "med. limit (+1 sig) " << r->GetExpectedUpperLimit(1) << endl;
    cout << endl;

    cout<<"========== plot result of the scan =========="<<endl;
    HypoTestInverterPlot* plot2 = 
        new HypoTestInverterPlot("HTI_Result_Plot", "CLs upper limit", r);
    TCanvas* c2 = new TCanvas("HypoTestInverter Scan"); 
    c2->SetLogy(false);
    plot2->Draw("2CL");
    c2->SaveAs("myCLsLimit.pdf");
    
}