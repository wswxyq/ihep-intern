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
#include<../RooClassFactory/RelativisticBW/RelativisticBW_wsw.cxx>
using namespace std;
using namespace RooFit ;




// $Id: $
void loadFilesToChain (TChain *chain, const std::string &filename){
    string line;
    ifstream source(filename.c_str(), ios::in);
    while ( getline (source, line) )
    {
        // you may want to clean line here...
        chain->AddFile( line.c_str() );
    }
}

void fit4BW_Mass_small()
{

    TProof::Open("");

	//set chain
    TChain *chain = new TChain("Pc2JpsipTuple/DecayTree");
	//load file
    //loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2011_2018_07_18.txt");
    //loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2012_2018_07_18.txt");
    //loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2015_2018_07_18.txt");
    //loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2016_2018_07_18.txt");
    //loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2017_2018_07_18.txt");
    chain->Add("/afs/ihep.ac.cn/users/j/jibo/gangadir/workspace/jhe/LocalXML/259/183/output/Tuple.root");
    chain->SetProof();
 
    TCut BMCut("B_DTF_M>0");
    TCut BDTCut("B_BDT>0.2");

    TCut TauCut("B_LOKI_FDS>49");

    TCut totCuts = 
        BMCut
        && BDTCut 
        && TauCut
             ;
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
	//define relativistic Breit Wigner distribution

	
	RooRealVar x("x", "mass spectrum", 4200, 4600, "MeV");
    //ref: page 131 LHCb-ANA-2018-043
    /*
    Discovery of narrow Pc(4312)+ → J/ψ p state in Λ0b → J/ψpK− decays, 
    and observation of two-peak structure of
    the Pc(4450)+
    */

    RooRealVar M4312("M4312", "M4312", 4312.1);
    RooRealVar M4440("M4440", "M4440", 4440);
    RooRealVar M4457("M4457", "M4457", 4457);
    RooRealVar Mx("Mx", "Mx", 4380, 4300, 4460);

    RooRealVar x1("x1", "para1", -1., 1.);
	RooRealVar x2("x2", "para2", -1., 1.);
	RooRealVar x3("x3", "para3", -1., 1.);

    RooRealVar gamma4312("gamma4312", "gamma4312", 10.0, 5.0, 15.0);
    RooRealVar gamma4440("gamma4440", "gamma4440", 20.0, 15.0, 40.0);
    RooRealVar gamma4457("gamma4457", "gamma4457", 10.0, 8.0, 12.0);
    RooRealVar gammax("gammax", "gammax", 100.0, 0.0, 200.0);

    RelativisticBW_wsw rtbw4312("rtbw4312", "rtbw4312", x, M4312, gamma4312);
    RelativisticBW_wsw rtbw4440("rtbw4440", "rtbw4440", x, M4440, gamma4440);
    RelativisticBW_wsw rtbw4457("rtbw4457", "rtbw4457", x, M4457, gamma4457);
    RelativisticBW_wsw rtbwx("rtbwx", "rtbwx", x, Mx, gammax);

    RooRealVar signal_frac_4312("signal_frac_4312", "signal_frac_4312", 0.1, 0., 1.);
    RooRealVar signal_frac_4440("signal_frac_4440", "signal_frac_4440", 0.1, 0., 1.);
    RooRealVar signal_frac_4457("signal_frac_4457", "signal_frac_4457", 0.1, 0., 1.0);
    RooRealVar signal_frac("signal_frac", "signal_frac", 0.1, 0., 1.);

    RooAddPdf signal("signal", "signal", RooArgList(rtbw4312, rtbw4440, rtbw4457, rtbwx),
                        RooArgList(signal_frac_4312, signal_frac_4440, signal_frac_4457));
    RooPolynomial background("background", "background", x, RooArgList(x1, x2, x3));

    RooAddPdf event("event", "event", RooArgList(signal, background), RooArgList(signal_frac));

    RooDataSet *ds=new RooDataSet("ds", "ds", RooArgSet(x), Import(*chain), 
		                Cut(totCuts));

    auto result= event.fitTo(*ds, RooFit::NumCPU(64), RooFit::Save(kTRUE), RooFit::Minos(kTRUE));



        
    
}
