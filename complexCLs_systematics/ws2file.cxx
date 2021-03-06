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

void ws2file();

int main() {

    ws2file();
    return 0;

}

void ws2file() {

    TCanvas* c1 = new TCanvas;

    // ------------------Define the model:  n ~ Poisson (s+b) and m ~ Poisson (tau*b), set tau = 1.-------------------------
    RooWorkspace* w = new RooWorkspace("w");

    RooRealVar B_DTF_M("B_DTF_M", "B_DTF_M", 4200, 4600, "MeV");
    RooGenericPdf rtbw4312("rtbw4312", "14.0*5.3*5.3*4312.0*4312.0/(22.0*(B_DTF_M*B_DTF_M-4312.0*4312.0)*(B_DTF_M*B_DTF_M-4312.0*4312.0)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(5.3*5.3)/(4312.0*4312.0))", RooArgList(B_DTF_M));
    RooGenericPdf rtbw4440("rtbw4400", "14.0*25.2*25.2*4440.2*4440.2/(22.0*(B_DTF_M*B_DTF_M-4440.2*4440.2)*(B_DTF_M*B_DTF_M-4440.2*4440.2)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(25.2*25.2)/(4440.2*4440.2))", RooArgList(B_DTF_M));
    RooGenericPdf rtbw4457("rtbw4457", "14.0*5.5*5.5*4456.6*4456.6/(22.0*(B_DTF_M*B_DTF_M-4456.6*4456.6)*(B_DTF_M*B_DTF_M-4456.6*4456.6)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(5.5*5.5)/(4456.6*4456.6))", RooArgList(B_DTF_M));
    RooGenericPdf rtbwx("rtbwx", "14.0*62.7*62.7*4394.7*4394.7/(22.0*(B_DTF_M*B_DTF_M-4394.7*4394.7)*(B_DTF_M*B_DTF_M-4394.7*4394.7)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(62.7*62.7)/(4394.7*4394.7))", RooArgList(B_DTF_M));
    //fraction of b-w
    RooRealVar signal_frac_4312("signal_frac_4312", "signal_frac_4312", 8.8406e-03, 0., 1.);
    RooRealVar signal_frac_4440("signal_frac_4440", "signal_frac_4440", 8.4540e-01, 0., 1.);
    RooRealVar signal_frac_4457("signal_frac_4457", "signal_frac_4457", 2.1666e-04, 0., 1.);
    RooRealVar signal_frac("signal_frac", "signal_frac", 2.1983e+00, 0, 100);
    RooRealVar background_frac("background_frac", "background_frac", 2.4580e+03, 2000., 3000.);
    //polynimial parameter
    //fit parameter first//////////////////////////
    RooRealVar x1("x1", "para1", 7.3718e+01);
    RooRealVar x2("x2", "para2", -1.0000e-04);
    RooRealVar x3("x3", "para3", -3.3642e-06);

    RooAddPdf smodel("smodel", "smodel", RooArgList(rtbw4312, rtbw4440, rtbw4457, rtbwx),
                     RooArgList(signal_frac_4312, signal_frac_4440, signal_frac_4457));
    RooPolynomial bmodel("bmodel", "bmodel", B_DTF_M, RooArgList(x1, x2, x3));

    RooAddPdf model0("model0", "model0", RooArgList(smodel, bmodel), RooArgList(signal_frac, background_frac));


    w->import(model0);
    // add constraint of background here.
    w->factory( "Gaussian::constr_b(global_b[2000.0,3000.0],background_frac,20)" );

    w->factory("PROD::model(model0,constr_b)");

    //w->import(s);
    w->defineSet("poi", "signal_frac");

    //this is necessary for setting global observables
    w->var("global_b")->setConstant(true);
    w->Print();


    cout<<"========== import a dataset with the observed data values B_DTF_M =========="<<endl;

    //set chain
    TChain *chain = new TChain("ReducedTree");
    //load file
    chain->Add("BDT_reduced.root");
    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(B_DTF_M), Import(*chain),
                                  Cut(""));
    w->defineSet("B_DTF_M", "B_DTF_M");
    w->defineSet("globObs", "global_b");
    w->defineSet("nuisance", "background_frac");
    w->import(*ds);
    w->Print();

    cout<<"========== now we need new model configs, with PDF=\"model\"=========="<<endl;
    ModelConfig b_modelNM("B_modelNM", w);

    b_modelNM.SetGlobalObservables(*w->set("globObs"));

    b_modelNM.SetPdf(*w->pdf("model"));
    b_modelNM.SetObservables(*w->set("B_DTF_M"));
    b_modelNM.SetParametersOfInterest(*w->set("poi"));
    b_modelNM.SetNuisanceParameters(*w->set("nuisance"));

    w->var("signal_frac")->setVal(0.0);
    b_modelNM.SetSnapshot(*w->set("poi"));     // sets up b hypothesis as s = 0


    cout<<"========== create the alternate (s+b) ModelConfig with given value of s =========="<<endl;
    double s_value = 0.001;
    ModelConfig sb_modelNM("S+B_modelNM", w);

    sb_modelNM.SetGlobalObservables(*w->set("globObs"));

    sb_modelNM.SetPdf(*w->pdf("model"));
    sb_modelNM.SetObservables(*w->set("B_DTF_M"));
    sb_modelNM.SetParametersOfInterest(*w->set("poi"));
    sb_modelNM.SetNuisanceParameters(*w->set("nuisance"));

    RooRealVar* poi = (RooRealVar*) sb_modelNM.GetParametersOfInterest()->first();
    w->var("signal_frac")->setVal(s_value);
    sb_modelNM.SetSnapshot(*w->set("poi"));  // set up sb hypothesis with given s

    w->import(b_modelNM);
    w->import(sb_modelNM);

    w->writeToFile("ws2file.root");

}

