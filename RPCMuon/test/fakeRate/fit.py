#!/usr/bin/env python

import sys, os
from ROOT import *
if os.path.exists("rootlogon.C"):
    gROOT.ProcessLine(".x rootlogon.C")

objs = []

def fit(hA, hB, c = None):
    nA = hA.Integral()
    nB = hB.Integral()
    nTotal = nA+nB

    massMin = hA.GetXaxis().GetXmin()
    massMax = hA.GetXaxis().GetXmax()

    mass = RooRealVar("mass", "mass", massMin, massMax, "GeV/c^{2}")
    mass.setBinning(RooBinning(hA.GetNbinsX(), massMin, massMax))
    hDataA = RooDataHist("hDataA", "mass", RooArgList(mass), hA)
    hDataB = RooDataHist("hDataB", "mass", RooArgList(mass), hB)

    ws = RooWorkspace("ws")
    getattr(ws, 'import')(mass)
    getattr(ws, 'import')(hDataA)
    getattr(ws, 'import')(hDataB)

    ws.factory("m0[%f, %f, %f]" % ((massMax+massMin)/2, massMin, massMax))
    ws.factory("Voigtian::sigA(mass, m0, w0[1e-2, 1e-5, 1e-1], sigmaA[1e-2, 1e-3, 1e-1])")
    ws.factory("Voigtian::sigB(mass, m0, w0, sigmaA)")
    #ws.factory("Voigtian::sigB(mass, m0, w0, sigmaB[1e-2, 1e-3, 1e+0])")
    ws.factory("Chebychev::bkgA(mass, {p0A[0, -10, 10], p1A[0, -10, 10]})")
    ws.factory("Chebychev::bkgB(mass, {p0B[0, -10, 10], p1B[0, -10, 10]})")
    ws.factory("ratio[0.001, 0, 1]")
    ws.factory("EXPR::nSigA('nSig*ratio', nSig[%f, 0, %f], ratio)" % (0.5*nTotal, 1.1*nTotal))
    ws.factory("EXPR::nSigB('nSig*(1-ratio)', nSig, ratio)")
    ws.factory("SUM::pdfA(nSigA*sigA, nBkgA[%f, 0, %f]*bkgA)" % (0.5*nA, 1.1*nA))
    ws.factory("SUM::pdfB(nSigB*sigB, nBkgB[%f, 0, %f]*bkgB)" % (0.5*nB, 1.1*nB))
    ws.factory("weight[1, 0, 1e9]")
    ws.factory("index[A,B]")
    ws.exportToCint()

    ws.index = ws.cat('index')
    weight = ws.var('weight')
    ws.pdfA = ws.pdf('pdfA')
    ws.pdfB = ws.pdf('pdfB')
 
    simPdf = RooSimultaneous("simPdf", "simPdf", ws.index);
    simPdf.addPdf(ws.pdfA, "A");
    simPdf.addPdf(ws.pdfB, "B");

    dataA = RooDataSet("dataA", "mass", RooArgSet(mass, weight), RooFit.WeightVar("weight"))
    dataB = RooDataSet("dataB", "mass", RooArgSet(mass, weight), RooFit.WeightVar("weight"))
    for i in range(hDataA.numEntries()):
        dataA.add(hDataA.get(i), hDataA.weight())
        dataB.add(hDataB.get(i), hDataB.weight())
    dataAB = RooDataSet("dataAB", "mass", RooArgSet(mass, weight), RooFit.Index(ws.index), 
                        RooFit.Import("A", dataA), RooFit.Import("B", dataB), RooFit.WeightVar(weight))

    result = simPdf.fitTo(dataAB, RooFit.Save(), RooFit.Extended())
  
    if c != None:
        frameA = mass.frame()
        dataAB.plotOn(frameA, RooFit.Cut("index==index::A"))
        proj = RooFit.ProjWData(RooArgSet(ws.index), dataAB)
        slice = RooFit.Slice(ws.index, "A")
        simPdf.plotOn(frameA, slice, proj, RooFit.LineColor(kBlue))
        simPdf.plotOn(frameA, slice, proj, RooFit.LineColor(kRed), RooFit.Components("bkgA"), RooFit.LineStyle(kDashed));

        frameB = mass.frame()
        slice = RooFit.Slice(ws.index, "B")
        dataAB.plotOn(frameB, RooFit.Cut("index==index::B"))
        simPdf.plotOn(frameB, slice, proj, RooFit.LineColor(kBlue))
        simPdf.plotOn(frameB, slice, proj, RooFit.LineColor(kRed), RooFit.Components("bkgB"), RooFit.LineStyle(kDashed))

        c.Divide(2,1)
        c.cd(1)
        frameA.Draw();
        c.cd(2)
        frameB.Draw();

    #result.Print("v");

    #cA->Print(Form("cA_%s_eta_RPC_%s.png", category, region));
    #cB->Print(Form("cB_%s_eta_RPC_%s.png", category, region));

    ratio = result.floatParsFinal().find('ratio')
    return ratio;

histFile = TFile("hist.root")
fitFile  = TFile("fit.root", "RECREATE")
for catName in [x.GetName() for x in histFile.GetListOfKeys()]:
    catDir = histFile.GetDirectory(catName)
    if catDir == None: continue
    outCatDir = fitFile.mkdir(catName)

    for varName in [x.GetName() for x in catDir.GetListOfKeys()]:
        varDir = catDir.GetDirectory(varName)
        if varDir == None: continue
        hFrame = varDir.Get("hFrame")
        if hFrame == None: continue
        outVarDir = outCatDir.mkdir(varName)
        outVarDir.cd()
        hFrame = hFrame.Clone()

        hFrame.SetMinimum(0)
        hFrame.SetMaximum(1)
        c = TCanvas("c_%s_%s" % (catName, varName), "%s %s" % (catName, varName), 500, 500)
        grp = TGraphErrors()
        grp.SetName("ratio")
        grp.SetTitle(catName)
        for bin in range(hFrame.GetNbinsX()):
            binName = 'bin_%d' % bin
            binDir = varDir.GetDirectory(binName)
            if binDir == None: continue

            hM_pass = binDir.Get("hM_pass")
            hM_fail = binDir.Get("hM_fail")
            cFitCanvas = TCanvas("cFit_%s_%s_%s" % (catName, varName, binName), "c %s %s %s" % (catName, varName, binName), 800, 400)
            ratio = fit(hM_pass, hM_fail, cFitCanvas)
            cFitCanvas.Write()
            objs.append(cFitCanvas)

            x    = hFrame.GetBinCenter(bin+1)
            dx   = hFrame.GetBinWidth(bin+1)/2
            y    = 100*ratio.getVal()
            dy   = abs(100*ratio.getError())
            grp.SetPoint(bin, x, y)
            grp.SetPointError(bin, dx, dy)
            if y > hFrame.GetMaximum(): hFrame.SetMaximum(y*1.1)

        c.cd()
        hFrame.Draw()
        grp.Draw("P")

        c.Write()
        hFrame.Write()
        grp.Write()
        objs.extend([c, hFrame, grp])

