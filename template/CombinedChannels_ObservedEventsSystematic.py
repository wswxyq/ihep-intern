#!/usr/bin/python
'''Simultanous fit to two different channels with same resonance but different background pdfs.
The total number of events is extracted by using the reconstruction efficiency in both channels and their branching ratios.
The efficiency and branching ratio have a systematic error.'''
from ROOT import *
import UpperLimit
from rootpy.interactive import wait

eff1 = 0.3
Br1 = 0.15
eff2 = 0.195
Br2 = 0.145

w = RooWorkspace('w')
w.factory('Exponential::bkg1_pdf(x[0,7], a1[-0.5,-2, 2])')
w.factory('Gaussian::sig_pdf(x,mass[2],sigma[0.3])')
#redefine nsig as product of efficiency and branching ratio (different for both channels) 
#and parameter of interest Nobs
w.factory('prod::effBr1(eff1[{0}], Br1[{1}])'.format(eff1, Br1))
w.factory('prod::nsig1(Nobs[280, 0, 1200], effBr1)')
w.factory('SUM::model1(nsig1*sig_pdf, nbkg1[1000,0,10000]*bkg1_pdf)')
#gaussian distribtuion of efficiency and branching ratio
w.factory('Gaussian::eff1Const(eff1g[{0}, 0, 1], eff1, eff1sig[0.01])'.format(eff1))
w.factory('Gaussian::br1Const(br1g[{0}, 0, 1], eff1, br1sig[0.05])'.format(Br1))
w.factory('PROD::constraintmodel1(model1, eff1Const, br1Const)') 
w.var('eff1g').setConstant(True)#do not vary in fit, treat as global variable later
w.var('br1g').setConstant(True)

w.factory('Chebychev::bkg2_pdf(x, {a2[-0.2,-2, 2]})')
w.factory('prod::effBr2(eff2[{0}], Br2[{1}])'.format(eff2, Br2))
w.factory('prod::nsig2(Nobs, effBr2)')
#gaussian distribtuion of efficiency and branching ratio
w.factory('Gaussian::eff2Const(eff2g[{0}, 0, 1], eff2, eff12ig[0.01])'.format(eff2))
w.factory('Gaussian::br2Const(br2g[{0}, 0, 1], eff2, br2sig[0.05])'.format(Br2))
w.factory('SUM::model2(nsig2*sig_pdf, nbkg2[200,0,10000]*bkg2_pdf)')
w.factory('PROD::constraintmodel2(model2, eff2Const, br2Const)') 
w.var('eff2g').setConstant(True)#do not vary in fit, treat as global variable later
w.var('br2g').setConstant(True)

#create discrete observable to label channels
w.factory('index[channel1,channel2]')
#create joint pdf
w.factory('SIMUL::jointModel(index,channel1=constraintmodel1,channel2=constraintmodel2)')
w.Print()

pdf = w.pdf('jointModel')
x = w.var('x')
index = w.cat('index')

#generate events with joint model
x.setBins(50)
data = pdf.generate(RooArgSet(x, index), RooFit.AllBinned())
data.SetName('data')
getattr(w,'import')(data)

#make example plots and fit both channels simultanously
plot1 = x.frame()
plot2 = x.frame()
plot3 = x.frame()
data.plotOn(plot1, RooFit.Cut('index==index::channel1'))
data.plotOn(plot2, RooFit.Cut('index==index::channel2'))
data.plotOn(plot3)

fitresult = pdf.fitTo(data, RooFit.Save(True))
fitresult.Print()
pdf.plotOn(plot1, RooFit.ProjWData(data), RooFit.Slice(w.cat('index'), 'channel1'))
pdf.plotOn(plot1, RooFit.ProjWData(data), RooFit.Components('bkg1_pdf'), RooFit.Slice(w.cat('index'), 'channel1'), RooFit.LineStyle(kDashed))
pdf.plotOn(plot2, RooFit.ProjWData(data), RooFit.Slice(w.cat('index'), 'channel2'))
pdf.plotOn(plot2, RooFit.ProjWData(data), RooFit.Components('bkg2_pdf'), RooFit.Slice(w.cat('index'), 'channel2'), RooFit.LineStyle(kDashed))
pdf.plotOn(plot3, RooFit.ProjWData(data))


canvas = TCanvas('simultanous fit',  'simultanous fit')
canvas.Divide(2, 2)
canvas.cd(1)
plot1.Draw()
canvas.cd(2)
plot2.Draw()
canvas.cd(3)
plot3.Draw()

#create ModelConfig
mc = RooStats.ModelConfig('ModelConfig', w)
mc.SetPdf(pdf)
mc.SetParametersOfInterest(RooArgSet(w.var('Nobs')))
mc.SetObservables(RooArgSet(w.var('x'), w.cat('index')))
w.defineSet('nuisParams','a1,nbkg1,nbkg2,a2')
w.defineSet('globalParams', 'eff1g,eff2g,br1g,br2g')
mc.SetNuisanceParameters(getattr(w,'set')('nuisParams'))
mc.SetGlobalObservables(getattr(w,'set')('globalParams'))
mc.SetSnapshot(RooArgSet(w.var('Nobs')))
getattr(w,'import')(mc)
w.Print()

#UpperLimit.get_significance(w, 'test', 'ModelConfig', True)
UpperLimit.restrict_nuissance_parameters(w, ['a1', 'nbkg1', 'nbkg2', 'a2'])
UpperLimit.get_upperlimit(w, 'test', 'ModelConfig', False, True, 15)

wait()