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
#include<RooHist.h>
#include<TCanvas.h>
#include<RooPlot.h>
#include<TMath.h>
#include<RooFitResult.h>
#include "../RooClassFactory/RelativisticBW/RelativisticBW_wsw.cxx"
using namespace std;
using namespace RooFit ;


void fit4BW_Mass_background()
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


    RooRealVar x1("x1", "para1", 74.3657, -100, 100);
    RooRealVar x2("x2", "para2", -6.47583e-05, -100, 100);
    RooRealVar x3("x3", "para3", -3.33708e-06, -10, 10);
    RooRealVar x4("x4", "para4", 0.48, -10., 10.);
    RooRealVar x5("x5", "para5", 0.48, -10., 10.);

    /*
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
    RooRealVar signal_frac("signal_frac", "signal_frac", 0.1, 0., 1.);

    RooAddPdf signal("signal", "signal", RooArgList(rtbw4312, rtbw4440, rtbw4457, rtbwx),
                        RooArgList(signal_frac_4312, signal_frac_4440, signal_frac_4457));
    */

    RooPolynomial background("background", "background", x, RooArgList(x1, x2, x3), 1);

    //RooAddPdf event("event", "event", RooArgList(signal, background), RooArgList(signal_frac));



    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(x, B_BDT, B_LOKI_FDS), Import(*chain),
                                  Cut(""));

    auto result= background.fitTo(*ds, RooFit::NumCPU(64), RooFit::Save(kTRUE), RooFit::Minos(kTRUE));

    RooPlot* xframe = x.frame(Title("Event p.d.f.")) ;

    ds->plotOn(xframe, Binning(20));
    background.plotOn(xframe) ;
    RooHist* hpull = xframe->pullHist() ;


    TCanvas* c = new TCanvas("total_plot","total_plot", 600, 600) ;

    xframe->Draw() ;
    c->Print("./fit_background.pdf");

    result->Print("v");

    // Access list of final fit parameter values
    cout << "final value of floating parameters" << endl ;
    result->floatParsFinal().Print("s");
    ds->Print();


}
