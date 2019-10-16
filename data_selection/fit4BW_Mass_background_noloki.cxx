#include<stdio.h>
#include<RooFit.h>
#include<TChain.h>
#include<TProof.h>
#include<TCut.h>
#include<TH1F.h>
#include<RooRealVar.h>
#include<RooPolynomial.h>
#include<RooGaussian.h>
#include<RooAddPdf.h>
#include<RooDataSet.h>
#include<RooDataHist.h>
#include<RooHist.h>
#include<TCanvas.h>
#include<RooPlot.h>
#include<TMath.h>
#include<RooFitResult.h>
#include "../RooClassFactory/RelativisticBW/RelativisticBW_wsw.cxx"
using namespace std;
using namespace RooFit ;


void fit4BW_Mass_background_noloki()
{

    //TProof::Open("");

	//set chain
    TChain *chain = new TChain("ReducedTree");
	//load file
    chain->Add("/scratchfs/others/wusw/BDT_reduced_noloki.root");
    //chain->SetProof();
 
    //RooRealVar B_DTF_M("B_DTF_M", "B_DTF_M", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar B_BDT("B_BDT", "B_BDT", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar B_LOKI_FDS("B_LOKI_FDS", "B_LOKI_FDS", -RooNumber::infinity(), RooNumber::infinity());

	RooRealVar x("B_DTF_M", "B_DTF_M", 4150, 4600, "MeV");

    RooRealVar x1("x1", "para1", 3.1740e+01, -1.e2, 1.e2);
	RooRealVar x2("x2", "para2", 1.2170e-02, -1.e2, 1.e2);
	RooRealVar x3("x3", "para3", -4.02805e-06, -1.e-01, 1.e-01);
	RooRealVar x4("x4", "para4", 1.1642e-09, -1.e-6, 1.e-6);


    RooRealVar sigma("sigma", "sigma", 200, 0, 400);
    RooRealVar sigma1("sigma1", "sigma1", 100, 0, 200);
    RooRealVar sigma2("sigma2", "sigma2", 150, 0, 200);
    RooRealVar mean("mean", "mean", 4.2000e+03, 4150, 4250);
    RooRealVar mean1("mean1", "mean1", 4.220e+03, 4150, 4300);
    RooRealVar mean2("mean2", "mean2", 4.500e+03, 4450, 4600);

    RooRealVar gauss_frac("gauss_frac", "gauss_frac", 7.8629e-01, 0., 1.);
    RooRealVar gauss_frac1("gauss_frac1", "gauss_frac1", 0.1, 0., 0.5);
    RooRealVar gauss_frac2("gauss_frac2", "gauss_frac2", 0.1, 0., 0.5);


    RooPolynomial poly("poly", "poly", x, RooArgList(x1, x2, x3, x4));
    RooGaussian gauss("gauss", "gauss", x, mean, sigma);
    RooGaussian gauss1("gauss1", "gauss1", x, mean1, sigma1);
    RooGaussian gauss2("gauss2", "gauss2", x, mean2, sigma2);

    RooAddPdf background("background", "background", RooArgList(gauss, gauss1, gauss2, poly), 
                RooArgList(gauss_frac, gauss_frac1, gauss_frac2));

    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(x, B_BDT, B_LOKI_FDS), Import(*chain), 
		                Cut(""));
    RooDataHist *dh=ds->binnedClone("dh", "dh");
                    

    auto result= background.fitTo(*dh, RooFit::NumCPU(60), RooFit::Save(kTRUE), RooFit::Minos(kTRUE));

	RooPlot* xframe = x.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = x.frame(Title(" ")) ;

	xframe_2->SetYTitle("pull distribution");

	ds->plotOn(xframe);
	background.plotOn(xframe);
	RooHist* hpull = xframe->pullHist() ;
    background.plotOn(xframe, Components(gauss),LineColor(kRed),LineStyle(kDashed));
    background.plotOn(xframe, Components(gauss1),LineColor(kGreen),LineStyle(kDashed));
    background.plotOn(xframe, Components(gauss2),LineColor(kYellow),LineStyle(kDashed));
    background.plotOn(xframe, Components(poly),LineColor(kOrange),LineStyle(kDashed));

    xframe_2->addPlotable(hpull, "P") ;
    
    TCanvas* c = new TCanvas("total_plot","total_plot", 1200, 1000) ;
	c->Divide(1,2) ;
    c->cd(1);
    xframe->Draw() ;
    c->cd(2);
	xframe_2->Draw() ;


    c->Print("./fit_background_noloki.pdf");

	std::ofstream myfile;
    myfile.open ("./fit_background_noloki.txt");
	myfile<<"number of entries in dataset: "<<ds->numEntries()<<endl;
	myfile<<"x1: "<<x1.getVal()<<endl;
	myfile<<"x2: "<<x2.getVal()<<endl;
	myfile<<"x3: "<<x3.getVal()<<endl;
    myfile.close();
    //x1.Print();
    //x2.Print();
    //x3.Print();
    result->Print("v");
    
    // Access list of final fit parameter values
    cout << "final value of floating parameters" << endl ;
    result->floatParsFinal().Print("s");
    ds->Print();


}
