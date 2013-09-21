#!/usr/bin/env python

import sys, os
from ROOT import *
from urllib import urlretrieve
if not os.path.exists('rootlogon.C'):
    urlretrieve('https://raw.github.com/cms-top-kr/tools/master/rootlogon.C', 'rootlogon.C')
gROOT.ProcessLine(".x rootlogon.C")

inputFileName = "hist.root"
outputFileName = "fit.root"
imageDirName = "image"

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

    ws.factory("m0A[91.2, 89, 92]")
    ws.factory("m0B[91.2, 89, 92]")
    ws.factory("Voigtian::sigA(mass, m0A, w0[2.495], sigmaA[1e-1, 1e-1, 3])")
    ws.factory("Voigtian::sigB(mass, m0B, w0, sigmaB[1e-1, 1e-1, 3])")
    ws.factory("Chebychev::bkgA(mass, {p0A[0, -5, 5], p1A[0, -5, 5]})")
    ws.factory("Chebychev::bkgB(mass, {p0B[0, -5, 5], p1B[0, -5, 5]})")
    ws.factory("efficiency[0.9, 0, 1]")
    ws.factory("EXPR::nSigA('nSig*efficiency', nSig[%f, 0, %f], efficiency)" % (0.5*nTotal, 1.1*nTotal))
    ws.factory("EXPR::nSigB('nSig*(1-efficiency)', nSig, efficiency)")
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

    simPdf.fitTo(dataAB)
    simPdf.fitTo(dataAB, RooFit.Extended())
    result = simPdf.fitTo(dataAB, RooFit.Save(), RooFit.Extended(), RooFit.Minos())
  
    if c != None:
        frameA = mass.frame()
        dataAB.plotOn(frameA, RooFit.Cut("index==index::A"))
        proj = RooFit.ProjWData(RooArgSet(ws.index), dataAB)
        slice = RooFit.Slice(ws.index, "A")
        simPdf.plotOn(frameA, slice, proj, RooFit.LineColor(kGreen))
        simPdf.plotOn(frameA, slice, proj, RooFit.LineColor(kGreen), RooFit.Components("bkgA"), RooFit.LineStyle(kDashed));

        frameB = mass.frame()
        slice = RooFit.Slice(ws.index, "B")
        dataAB.plotOn(frameB, RooFit.Cut("index==index::B"))
        simPdf.plotOn(frameB, slice, proj, RooFit.LineColor(kRed))
        simPdf.plotOn(frameB, slice, proj, RooFit.LineColor(kRed), RooFit.Components("bkgB"), RooFit.LineStyle(kDashed))

        c.Divide(2,1)
        c.cd(1)
        frameA.Draw();
        c.cd(2)
        frameB.Draw();

    efficiency = result.floatParsFinal().find('efficiency')
    return efficiency;

summaryText = ""
if not os.path.isdir(imageDirName) :os.mkdir(imageDirName)
histFile = TFile(inputFileName)
fitFile  = TFile(outputFileName, "RECREATE")
for catName in [x.GetName() for x in histFile.GetListOfKeys()]:
    catDir = histFile.GetDirectory(catName)
    if catDir == None: continue
    outCatDir = fitFile.mkdir(catName)
    if not os.path.isdir("%s/%s" % (imageDirName, catName)): os.mkdir("%s/%s" % (imageDirName, catName))

    for varName in [x.GetName() for x in catDir.GetListOfKeys()]:
        varDir = catDir.GetDirectory(varName)
        if varDir == None: continue
        hFrame = varDir.Get("hFrame")
        if hFrame == None: continue
        outVarDir = outCatDir.mkdir(varName)
        outVarDir.cd()
        hFrame = hFrame.Clone()

        if not os.path.isdir("%s/%s/%s" % (imageDirName, catName, varName)): os.mkdir("%s/%s/%s" % (imageDirName, catName, varName))

        tab_header = []
        tab_values = []

        hFrame.SetMinimum(0)
        hFrame.SetMaximum(1)
        hFrame.GetYaxis().SetTitle("Fake rate (%)")
        c = TCanvas("c_%s_%s" % (catName, varName), "%s %s" % (catName, varName), 500, 500)
        grp = TGraphAsymmErrors()
        grp.SetName("efficiency")
        grp.SetTitle(catName)
        for bin in range(hFrame.GetNbinsX()):
            binName = 'bin_%d' % bin
            binDir = varDir.GetDirectory(binName)
            if binDir == None: continue

            hM_pass = binDir.Get("hM_pass")
            hM_fail = binDir.Get("hM_fail")
            cFitCanvas = TCanvas("cFit_%s_%s_%s" % (catName, varName, binName), "c %s %s %s" % (catName, varName, binName), 800, 400)
            efficiency = fit(hM_pass, hM_fail, cFitCanvas)
            cFitCanvas.Write()
            objs.append(cFitCanvas)
            cFitCanvas.Print("%s/%s/%s/%s.png" % (imageDirName, catName, varName, cFitCanvas.GetName()))
            cFitCanvas.Print("%s/%s/%s/%s.pdf" % (imageDirName, catName, varName, cFitCanvas.GetName()))

            x    = hFrame.GetBinCenter(bin+1)
            dx   = hFrame.GetBinWidth(bin+1)/2
            y    = 100*efficiency.getVal()
            dyLo = abs(100*efficiency.getErrorLo())
            dyHi = abs(100*efficiency.getErrorHi())
            grp.SetPoint(bin, x, y)
            grp.SetPointError(bin, dx, dx, dyLo, dyHi)
            if y > hFrame.GetMaximum(): hFrame.SetMaximum(y*1.1)

            tab_header.append("$%.2f-%.2f$" % (x-dx, x+dx))
            tab_values.append("$%.1f^{+%.1f}_{-%.1f}$" % (y, dyLo, dyHi))
        summaryText += ("    %s &" % varName ) + " & ".join(tab_header) + "\\\\" + "\n"
        summaryText += ("    %s &" % catName ) + " & ".join(tab_values) + "\\\\" + "\n"

        c.cd()
        hFrame.Draw()
        grp.Draw("P")

        c.Write()
        hFrame.Write()
        grp.Write()
        objs.extend([c, hFrame, grp])

        c.Print("%s/%s/%s/%s.png" % (imageDirName, catName, varName, c.GetName()))
        c.Print("%s/%s/%s/%s.pdf" % (imageDirName, catName, varName, c.GetName()))

print summaryText
