#include<stdio.h>
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "TCanvas.h"
using namespace std;
void Cal_CLs(){
    //Open the ROOT file
    TFile* f = TFile::Open("CLs_new.root") ;
    cout<<"Retrieve the workspace"<<endl;
    RooWorkspace* w = (RooWorkspace*) f->Get("w") ;
    w->Print();
    cout<<"Retrieve the ModelConfig for the S+B hypothesis"<<endl;
    
    RooAbsData* data = w->data("obsData") ;
    RooStats::ModelConfig* sbModel = (RooStats::ModelConfig*) w->obj("model") ;
    cout<<"Construct a ModelConfig for the B-only hypothesis"<<endl;
    RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) sbModel->Clone("bmodel") ;
    //RooStats::ModelConfig* bModel = (RooStats::ModelConfig*) w->obj("bmodel") ;
    cout<<"Set value POI parameter to zero"<<endl;
    RooRealVar* poi = (RooRealVar*) bModel->GetParametersOfInterest()->first();
    poi-> setVal(0) ;
    cout<<"Configure bModel to encode current poi=0 scenario as its hypothesis"<<endl;
    bModel->SetSnapshot( *poi  );
    cout<<"Construct an hypothesis p-value calculator"<<endl;
    RooStats::AsymptoticCalculator  asympCalc(*data, *bModel, *sbModel);
    cout<<"Configure calculator for a limit (=one-sided interval)"<<endl;
    asympCalc.SetOneSided(true);
    cout<<"Construct an hypothesis test inverter"<<endl;
    RooStats::HypoTestInverter inverter(asympCalc);
    cout<<"Statistical configuration of hypothesis test inverter"<<endl;
    inverter.SetConfidenceLevel(0.90);
    inverter.UseCLs(true);
    cout<<"Technical configuration of hypothesis test inverter"<<endl;
    inverter.SetVerbose(false);
    inverter.SetFixedScan(50,0.0,6.0); 
    cout<<" set number of points , xmin and xmax"<<endl;
    cout<<"Perform calculation of limit"<<endl;
    RooStats::HypoTestInverterResult* result =  inverter.GetInterval();
    cout<<"Print observed limit"<<endl;
    cout << 100*inverter.ConfidenceLevel() << "%  upper limit : " << result->UpperLimit() << endl;
    cout<<"compute expected limit"<<endl;
    std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << result->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << result->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << result->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << result->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << result->GetExpectedUpperLimit(2) << std::endl;  
    cout<<"Use the visualization tool of the PLC to show how the interval was calculated"<<endl;
    RooStats::HypoTestInverterPlot* plot = new RooStats::HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",result);
    
    //plot->Draw("CLb 2CL");  // plot also CLb and CLs+b
    //c1->Draw() ;


    //return 0;


}