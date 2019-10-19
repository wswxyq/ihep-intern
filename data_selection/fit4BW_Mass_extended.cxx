#include<stdio.h>
#include<RooFit.h>
#include<TChain.h>
#include<TProof.h>
#include<TCut.h>
#include<TH1F.h>
#include<RooRealVar.h>
#include<RooPolynomial.h>
#include<RooAddPdf.h>
#include<RooDataSet.h>
#include<RooGaussian.h>
#include<RooHist.h>
#include<TCanvas.h>
#include<RooPlot.h>
#include<TMath.h>
#include<RooFitResult.h>
#include "../RooClassFactory/RelativisticBW/RelativisticBW_wsw.cxx"
using namespace std;
using namespace RooFit;


void fit4BW_Mass_extended()
{

    //TProof::Open("");

	//set chain
    TChain *chain = new TChain("ReducedTree");
	//load file
    chain->Add("/scratchfs/others/wusw/BDT_reduced.root");
    //chain->SetProof();
 
    //RooRealVar B_DTF_M("B_DTF_M", "B_DTF_M", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar B_BDT("B_BDT", "B_BDT", -RooNumber::infinity(), RooNumber::infinity());
    RooRealVar B_LOKI_FDS("B_LOKI_FDS", "B_LOKI_FDS", -RooNumber::infinity(), RooNumber::infinity());



	/*
    //...mass range
    const Double_t MassMin = 4000;
    const Double_t MassMax = 6000;
    const Int_t BinNum = 200;

    //...mass resolution
    const Double_t MsigmaAve = 3. ;
    const Double_t MsigmaMin = -10. ;
    const Double_t MsigmaMax = 300. ;

    TH1F *h10 = new TH1F("h10", "M(J/#psip)", BinNum, MassMin, MassMax);

    chain->Project("h10", "B_DTF_M", totCuts );

    h10->Draw("E");
	*/

	
	RooRealVar x("B_DTF_M", "B_DTF_M", 4200, 4600, "MeV");
    //ref: page 131 LHCb-ANA-2018-043
    /*
    Discovery of narrow Pc(4312)+ → J/ψ p state in Λ0b → J/ψpK− decays, 
    and observation of two-peak structure of
    the Pc(4450)+
    */

    RooRealVar M4312("M4312", "M4312", 4312.0);
    RooRealVar M4440("M4440", "M4440", 4440.2);
    RooRealVar M4457("M4457", "M4457", 4456.6);
    RooRealVar Mx("Mx", "Mx", 4394.7);
    
    RooRealVar x1("x1", "para1", 74.3657, -80, 150);
	RooRealVar x2("x2", "para2", -6.47583e-05, -10.e-05, -4.e-05);
	RooRealVar x3("x3", "para3", -3.33708e-06, -10.e-06, 2.0e-06);
    

    RooRealVar gamma4312("gamma4312", "gamma4312", 5.3);
    RooRealVar gamma4440("gamma4440", "gamma4440", 25.2);
    RooRealVar gamma4457("gamma4457", "gamma4457", 5.5);
    RooRealVar gammax("gammax", "gammax", 62.7);

    RelativisticBW_wsw rtbw4312("rtbw4312", "rtbw4312", x, M4312, gamma4312);
    RelativisticBW_wsw rtbw4440("rtbw4440", "rtbw4440", x, M4440, gamma4440);
    RelativisticBW_wsw rtbw4457("rtbw4457", "rtbw4457", x, M4457, gamma4457);
    RelativisticBW_wsw rtbwx("rtbwx", "rtbwx", x, Mx, gammax);

    RooRealVar signal_frac_4312("signal_frac_4312", "signal_frac_4312", 0.1, 0., 1.);
    RooRealVar signal_frac_4440("signal_frac_4440", "signal_frac_4440", 0.1, 0., 1.);
    RooRealVar signal_frac_4457("signal_frac_4457", "signal_frac_4457", 0.1, 0., 1.);

    RooRealVar signal_num("signal_num", "signal_num", 0.1, 0., 100.);
    RooRealVar background_num("background_num", "background_num", 2458., 0., 3000.);

    RooAddPdf signal("signal", "signal", RooArgList(rtbw4312, rtbw4440, rtbw4457, rtbwx),
                        RooArgList(signal_frac_4312, signal_frac_4440, signal_frac_4457));
    RooPolynomial background("background", "background", x, RooArgList(x1, x2, x3));

    RooAddPdf event("event", "event", RooArgList(signal, background), RooArgList(signal_num, background_num));



    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(x, B_BDT, B_LOKI_FDS), Import(*chain), 
		                Cut(""));

    auto result= event.fitTo(*ds, RooFit::NumCPU(60), RooFit::Save(kTRUE), RooFit::Minos(kTRUE));

	RooPlot* xframe = x.frame(Title("Event p.d.f.")) ;
	RooPlot* xframe_2 = x.frame(Title(" ")) ;

	xframe_2->SetYTitle("pull distribution");


	ds->plotOn(xframe);
	event.plotOn(xframe) ;
	RooHist* hpull = xframe->pullHist() ;
	event.plotOn(xframe,Components(rtbw4312),LineColor(kPink),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(rtbw4440),LineColor(kGreen),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(rtbw4457),LineColor(kYellow),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(rtbwx),LineColor(kOrange),LineStyle(kDashed)) ;
	event.plotOn(xframe,Components(background),LineColor(kCyan),LineStyle(kDashed)) ;
    xframe_2->addPlotable(hpull, "P") ;

    TCanvas* c = new TCanvas("total_plot","total_plot", 1200, 1000) ;
	c->Divide(1,2) ;
	c->cd(1) ; 
    xframe->Draw() ;

	c->cd(2) ; 
	//gPad->SetLeftMargin(0.15) ; 
	//xframe_2->GetYaxis()->SetTitleOffset(1.6) ; 
	//xframe_2->GetYaxis()->SetRangeUser(-5., 5.);
	xframe_2->Draw() ;
    c->Print("./fit_all_event_extended.pdf");



	std::ofstream myfile;
    myfile.open ("./fit_all_event_extended.txt");
	myfile<<"number of entries in dataset: "<<ds->numEntries()<<endl;
	myfile<<"x1: "<<x1.getVal()<<endl;
	myfile<<"x2: "<<x2.getVal()<<endl;
	myfile<<"x3: "<<x3.getVal()<<endl;
	myfile<<"signal_frac_4312: "<<signal_frac_4312.getVal()<<endl;
	myfile<<"signal_frac_4440: "<<signal_frac_4440.getVal()<<endl;
	myfile<<"signal_frac_4457: "<<signal_frac_4457.getVal()<<endl;
	myfile<<"signal_num: "<<signal_num.getVal()<<endl;
	myfile<<"background_num: "<<background_num.getVal()<<endl;
    myfile.close();
    x1.Print();
    x2.Print();
    x3.Print();
    signal_frac_4312.Print();
    signal_frac_4440.Print();
    signal_frac_4457.Print();
    signal_num.Print();
    background_num.Print();

    result->Print("v");
    
}
