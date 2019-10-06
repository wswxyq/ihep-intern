#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooDataSet.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TCut.h"
#include "TProof.h"

using namespace std;

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



int cut0()
{

    //TProof::Open("");

	//set chain
    TChain *chain = new TChain("Pc2JpsipTuple/DecayTree");
	//load file
    loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2011_2018_07_18.txt");
    loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2012_2018_07_18.txt");
    loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2015_2018_07_18.txt");
    loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2016_2018_07_18.txt");
    loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2017_2018_07_18.txt");

    //chain->SetProof();
 
    TCut BMCut("B_DTF_M>0");
    TCut BDTCut("B_BDT>0.2");

    TCut TauCut("B_LOKI_FDS>49");

    TCut totCuts = 
        BMCut
        && BDTCut 
        && TauCut
             ;


    TFile *newfile = new TFile("/scratchfs/others/wusw/new.root","recreate");


	//creat newtree
	chain->SetBranchStatus("*",0);
    chain->SetBranchStatus("B_DTF_M",1);
    chain->SetBranchStatus("B_BDT",1);
    chain->SetBranchStatus("B_LOKI_FDS",1);
	TTree *newtree = (TTree*)chain->CopyTree(totCuts);

	newfile->cd();
    newtree->Write("ReducedTree");
	return 0;
}
