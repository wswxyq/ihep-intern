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

void drawMass()
{

  TProof::Open("");

  TChain *chain = new TChain("Pc2JpsipTuple/DecayTree");

  loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2011_2018_07_18.txt");
  loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2012_2018_07_18.txt");
  loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2015_2018_07_18.txt");
  loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2016_2018_07_18.txt");
  loadFilesToChain(chain, "/afs/ihep.ac.cn/users/j/jibo/public/GangaScripts/Pc2JpsiP_BDTD_2017_2018_07_18.txt");

  chain->SetProof();
 
  TCut BMCut("B_DTF_M>0");
  TCut BDTCut("B_BDT>0.2");

  TCut TauCut("B_LOKI_FDS>49");

  TCut totCuts = 
    BMCut
    && BDTCut 
    && TauCut
       ;
 
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

    
  
}
