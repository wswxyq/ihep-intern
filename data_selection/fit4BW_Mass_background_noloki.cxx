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

    RooRealVar x1("x1", "para1", 3.1740e+01, 1.e+01, 5.e+01);
	RooRealVar x2("x2", "para2", 1.2170e-02, -1.e0, 1.e0);
	RooRealVar x3("x3", "para3", -4.02805e-06, -1.e-01, 1.e-01);
	RooRealVar x4("x4", "para4", 1.1642e-09, -1.e-6, 1.e-6);


    RooRealVar sigma("sigma", "sigma", 200, 0, 400);
    RooRealVar mean("mean", "mean", 4.2000e+03, 4150, 4250);

    RooRealVar gauss_frac("gauss_frac", "gauss_frac", 7.8629e-01, 0., 1.);


    RooPolynomial poly("poly", "poly", x, RooArgList(x1, x2, x3, x4));
    RooGaussian gauss("gauss", "gauss", x, mean, sigma);

    RooAddPdf background("background", "background", RooArgList(gauss, poly), RooArgList(gauss_frac));

    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(x, B_BDT, B_LOKI_FDS), Import(*chain), 
		                Cut(""));
    RooDataHist *dh=ds->binnedClone("dh", "dh");
                    

    auto result= background.fitTo(*dh, RooFit::NumCPU(60), RooFit::Save(kTRUE), RooFit::Minos(kTRUE));

	RooPlot* xframe = x.frame(Title("Event p.d.f.")) ;

	ds->plotOn(xframe);
	background.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;

    
    TCanvas* c = new TCanvas("total_plot","total_plot", 600, 600) ;

    xframe->Draw() ;
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
