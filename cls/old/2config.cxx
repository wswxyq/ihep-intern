#include<stdio.h>
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
#include "../RooClassFactory/RelativisticBW/RelativisticBW_wsw.h"

#include <cassert>
using namespace RooStats;
using namespace std;
// macro to fit Higgs to gg spectrum

using namespace RooFit;

void config() {

    //set chain
    TChain *chain = new TChain("ReducedTree");
	//load file
    chain->Add("../data_selection/BDT_reduced.root");

    // make the RooFit model     
    RooWorkspace w("w");    
    w.factory("B_DTF_M[4200,4600]");  

    // invariant mass        
    w.factory("M4312[4312]");
    w.factory("M4440[4440]");
    w.factory("M4457[4457]");
    w.factory("Mx[4394.7]");

    //width
    w.factory("gamma4312[5.3]");
    w.factory("gamma4440[25.2]");
    w.factory("gamma4457[5.5]");
    w.factory("gammax[62.7]");

    w.factory("expr::rtbw4312('14.0*gamma4312*gamma4312*M4312*M4312/\
                (22.0*(B_DTF_M*B_DTF_M-M4312*M4312)*(B_DTF_M*B_DTF_M-M4312*M4312)\
                +B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(gamma4312*gamma4312)/(M4312*M4312))', B_DTF_M, M4312, gamma4312)");
    //cout<<"llll\n"<<endl;
    w.factory("expr::rtbw4440('14.0*gamma4440*gamma4440*M4440*M4440/(22.0*(B_DTF_M*B_DTF_M-M4440*M4440)*(B_DTF_M*B_DTF_M-M4440*M4440)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(gamma4440*gamma4440)/(M4440*M4440))', B_DTF_M, M4440, gamma4440)");
    w.factory("expr::rtbw4457('14.0*gamma4457*gamma4457*M4457*M4457/(22.0*(B_DTF_M*B_DTF_M-M4457*M4457)*(B_DTF_M*B_DTF_M-M4457*M4457)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(gamma4457*gamma4457)/(M4457*M4457))', B_DTF_M, M4457, gamma4457)");
    w.factory("expr::rtbwx('14.0*gammax*gammax*Mx*Mx/(22.0*(B_DTF_M*B_DTF_M-Mx*Mx)*(B_DTF_M*B_DTF_M-Mx*Mx)+B_DTF_M*B_DTF_M*B_DTF_M*B_DTF_M*(gammax*gammax)/(Mx*Mx))', B_DTF_M, Mx, gammax)");
    

    w.factory("signal_frac_4312[0.1, 0., 1.]");
    w.factory("signal_frac_4440[0.1, 0., 1.]");
    w.factory("signal_frac_4457[0.1, 0., 1.]");
    w.factory("signal_frac_x[0.1, 0., 1.]");


    w.factory("nbackground[4658, 0, 4658]");    
    w.factory("nsignal[0, 0, 4658]");

    //w.factory("Exponential::z1(x, a1[-1,-10,0])");    
    /*w.var("nbackground")->setVal(nevt);
    w.var("nbackground")->setMin(0.1*nevt);
    w.var("nbackground")->setMax(10*nevt);*/

    // create bmodel
    w.factory("x1[ 74.3657, -80, 80]");
    w.factory("x2[-6.47583e-05, -10.e-05, -2.e-05]");
    w.factory("x3[-3.33708e-06, -10.e-06, 4.0e-06]");
    w.factory("RooPolynomial::bmodel(x1, x2, x3)");
    //w.factory("expr::z('-(a1*x/100 + a2*(x/100)^2)', a1, a2, x)");

    // signal model   
    //w.factory("mass[%f, %f, %f]' % (massguess, massmin, massmax))
    RooAbsPdf * rtbw4312 = w.pdf("rtbw4312");
    RooAbsPdf * rtbw4440 = w.pdf("rtbw4440");
    RooAbsPdf * rtbw4457 = w.pdf("rtbw4457");
    RooAbsPdf * rtbwx = w.pdf("rtbwx");
    cout<<"1\n"<<endl;
    w.factory("SUM::smodel(signal_frac_4312*rtbw4312,signal_frac_4440*rtbw4440,signal_frac_4457*rtbw4457,rtbwx)");
    cout<<"2\n"<<endl;
    w.factory("SUM::model(nbackground*bmodel,nsignal*smodel)");
    RooAbsPdf * smodel = w.pdf("smodel");
    RooAbsPdf * bmodel = w.pdf("bmodel");
    RooAbsPdf * model = w.pdf("model");
    

    // create RooDataSet
    RooDataSet data("data","data",*w.var("B_DTF_M"),Import(*chain) );
    //RooDataSet ds("ds", "ds", RooArgSet(x), Import(*chain), 
	//	                Cut(""));
    //new
    RooStats::ModelConfig mc("ModelConfig",&w);
    //mc.SetPdf(*w.pdf("model"));
    //mc.SetPdf(*w.pdf("bmodel"));
    //mc.SetPdf(*bmodel);
    //mc.SetPdf(*smodel);
    //mc.SetPdf(*model);

    mc.SetParametersOfInterest(*w.var("nsignal"));
    mc.SetObservables(*w.var("B_DTF_M"));
    // define set of nuisance parameters
    w.defineSet("nuisParams","x1,x2,x3,nbackground,signal_frac_4312,signal_frac_4440,signal_frac_4457,signal_frac_x");

    mc.SetNuisanceParameters(*w.set("nuisParams"));


    //w.var("mass")->setConstant(true);
    //w.var("width")->setConstant(true);



    // import the ModelConfig in the workspace
    w.import(mc); 

    // import the data 
    // for performance reason we transform the data in a binned data set
    TH1 * h1 = data.createHistogram("B_DTF_M");
    RooDataHist binData("obsData","obsData",*w.var("B_DTF_M"), h1); 
    w.import(binData);
    //ClassDef(MyTest, 1)   // 
    // write workspace in the file 
    w.writeToFile("CLs_Model.root");     

}

