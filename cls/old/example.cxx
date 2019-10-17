using namespace RooFit;

void fitHgg() {

   // read from the file and create a ROOT tree

   TTree tree("tree","tree");
   int nevt = tree.ReadFile("Hgg.txt","x");
   if (nevt <= 0) {       Error("fitHgg","Error reading data from input file ");       return;    }    std::cout << "Read " << nevt << " from the file " << std::endl;    // make the RooFit model     RooWorkspace w("w");    w.factory("x[110,160]");  // invariant mass        w.factory("nbackground[10000, 0, 10000]");    //w.factory("Exponential::z1(x, a1[-1,-10,0])");    w.var("nbackground")->setVal(nevt);
   w.var("nbackground")->setMin(0.1*nevt);
   w.var("nbackground")->setMax(10*nevt);

   // create exponential model as two components
   w.factory("a1[ 7.5, -500, 500]");
   w.factory("a2[-1.5, -500, 500]");
   w.factory("expr::z('-(a1*x/100 + a2*(x/100)^2)', a1, a2, x)");
   w.factory("Exponential::bmodel(z, 1)");

   // signal model   
   w.factory("nsignal[100, 0.0, 1000.0]");
   //w.factory("mass[%f, %f, %f]' % (massguess, massmin, massmax))
   w.factory("mass[130, 110, 150]");
   w.factory("width[1, 0.5, 5]");
   w.factory("Gaussian::smodel(x, mass, width)");
   RooAbsPdf * smodel = w.pdf("smodel");

   w.factory("SUM::model(nbackground*bmodel, nsignal*smodel)");
   RooAbsPdf * model = w.pdf("model");

   // create RooDataSet
   RooDataSet data("data","data",*w.var("x"),Import(tree) );
   
   RooFitResult * r = model->fitTo(data, Minimizer("Minuit2"),Save(true), Offset(true));

   // plot data and function

   RooPlot * plot = w.var("x")->frame();
   data.plotOn(plot);
   model->plotOn(plot);
   model->plotOn(plot, Components("bmodel"),LineStyle(kDashed));
   model->plotOn(plot, Components("smodel"),LineColor(kRed));

   plot->Draw();

}