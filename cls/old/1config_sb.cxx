#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"

#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "RooBreitWigner.h"
#include "RooWorkspace.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumEventsTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "../RooClassFactory/RelativisticBW/RelativisticBW_wsw.cxx"

#include <cassert>
using namespace std;
using namespace RooFit ;
using namespace RooStats;


void config_sb(){
    RooWorkspace w("w");    

    w.importClassCode(RelativisticBW_wsw::Class(), kTRUE);
    //w.addClassDeclImportDir("../RooClassFactory/RelativisticBW/");
    //w.addClassImplImportDir("../RooClassFactory/RelativisticBW/");
    //set chain
    TChain *chain = new TChain("ReducedTree");
	//load file
    chain->Add("../data_selection/BDT_reduced.root");

    RooRealVar nsignal("nsignal", "nsignal", 0, 0, 4658);
    RooRealVar nbackground("nbackground", "nbackground", 4658, 0, 4658);

	RooRealVar B_DTF_M("B_DTF_M", "B_DTF_M", 4200, 4600, "MeV");

    //invariant mass
    RooRealVar M4312("M4312", "M4312", 4312.0);
    RooRealVar M4440("M4440", "M4440", 4440.2);
    RooRealVar M4457("M4457", "M4457", 4456.6);
    RooRealVar Mx("Mx", "Mx", 4394.7);
    
    //polynimial parameter
    RooRealVar x1("x1", "para1", 74.3657, -80, 80);
	RooRealVar x2("x2", "para2", -6.47583e-05, -10.e-05, -4.e-05);
	RooRealVar x3("x3", "para3", -3.33708e-06, -10.e-06, 2.0e-06);
    
    //width
    RooRealVar gamma4312("gamma4312", "gamma4312", 5.3);
    RooRealVar gamma4440("gamma4440", "gamma4440", 25.2);
    RooRealVar gamma4457("gamma4457", "gamma4457", 5.5);
    RooRealVar gammax("gammax", "gammax", 62.7);

    //breit wigner function
    RooBreitWigner rtbw4312("rtbw4312", "rtbw4312", B_DTF_M, M4312, gamma4312);
    RooBreitWigner rtbw4440("rtbw4440", "rtbw4440", B_DTF_M, M4440, gamma4440);
    RooBreitWigner rtbw4457("rtbw4457", "rtbw4457", B_DTF_M, M4457, gamma4457);
    RooBreitWigner rtbwx("rtbwx", "rtbwx", B_DTF_M, Mx, gammax);

    //fraction of b-w
    RooRealVar signal_frac_4312("signal_frac_4312", "signal_frac_4312", 0.1, 0., 1.);
    RooRealVar signal_frac_4440("signal_frac_4440", "signal_frac_4440", 0.1, 0., 1.);
    RooRealVar signal_frac_4457("signal_frac_4457", "signal_frac_4457", 0.1, 0., 1.);
    RooRealVar signal_frac("signal_frac", "signal_frac", 0.1, 0., 1.);

    RooAddPdf smodel("smodel", "smodel", RooArgList(rtbw4312, rtbw4440, rtbw4457, rtbwx),
                        RooArgList(signal_frac_4312, signal_frac_4440, signal_frac_4457));
    RooPolynomial bmodel("bmodel", "bmodel", B_DTF_M, RooArgList(x1, x2, x3));

    RooAddPdf model("model", "model", RooArgList(smodel, bmodel), RooArgList(signal_frac));



    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(B_DTF_M), Import(*chain), 
		                Cut(""));



    RooStats::ModelConfig mc("ModelConfig",&w);

    mc.SetPdf(smodel);
    mc.SetPdf(bmodel);
    mc.SetPdf(model);

    mc.SetParametersOfInterest(nsignal);
    mc.SetObservables(B_DTF_M);
    mc.SetNuisanceParameters(RooArgSet(x1, x2, x3, signal_frac_4312,
            signal_frac_4440, signal_frac_4457, signal_frac, nbackground, "nuisParams"));

    // import the ModelConfig in the workspace
    w.import(mc); 

    // import the data 
    // for performance reason we transform the data in a binned data set
    TH1 * h1 = ds->createHistogram("B_DTF_M");
    RooDataHist binData("obsData","obsData",B_DTF_M, h1); 
    w.import(binData);

    // write workspace in the file 
    w.writeToFile("CLs_new.root");


}