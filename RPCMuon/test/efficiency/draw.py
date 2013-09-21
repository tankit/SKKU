#!/usr/bin/env python

import sys, os
from ROOT import *
from urllib import urlretrieve
if not os.path.exists('rootlogon.C'):
    urlretrieve('https://raw.github.com/cms-top-kr/tools/master/rootlogon.C', 'rootlogon.C')
gROOT.ProcessLine(".x rootlogon.C")
gROOT.ForceStyle()

objects = []
f = TFile("tnpHist.root")

muonIdNames = ["RPCMuLoose", "RPCMuMedium", "RPCMuTight"]
varNames = ["pt", "abseta"]
colors = [kRed, kRed+1, kMagenta, kBlue, kAzure, kBlack]

modDir = f.GetDirectory("muonEffs")
for varName in varNames:
    c = TCanvas("c_%s" % varName, varName, 500, 500)
    leg = TLegend(0.6, 0.25, 0.92, 0.50)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    hFrame = None
    for i, muonIdName in enumerate(muonIdNames):
        varDir = modDir.GetDirectory("%s_%s" % (muonIdName, varName))
        if varDir == None: continue
        if hFrame == None:
            hFrame = varDir.Get("hFrame")
            hFrame.GetYaxis().SetTitle("Efficiency")
            hFrame.GetYaxis().SetRangeUser(0.50, 1.00)
            hFrame.Draw()

            tab_header = [varName]
            for b in range(hFrame.GetNbinsX()):
                binLo = hFrame.GetXaxis().GetBinLowEdge(b+1)
                binHi = binLo+hFrame.GetXaxis().GetBinWidth(b+1)
                tab_header.append("$%.2f-%.2f$" % (binLo, binHi))
            print "\hline\hline"
            print " & ".join(tab_header), "\\\\"
            print "\hline"

            objects.append(hFrame)
        hEffic = varDir.Get("hEfficiency")
        if hEffic == None: continue

        color = colors[i]
        hEffic.SetLineColor(color)
        hEffic.SetMarkerColor(color)
        hEffic.Draw("P")

        tab_content = []
        for b in range(hEffic.GetN()):
            y = hEffic.GetY()[b]
            eyLo = hEffic.GetEYlow()[b]
            eyHi = hEffic.GetEYhigh()[b]
            tab_content.append("$%.3f^{+%.3f}_{-%.3f}$" % (y, eyHi, eyLo))
        print muonIdName, "&",  " & ".join(tab_content), "\\\\"

        leg.AddEntry(hEffic, muonIdName, "lp")

        objects.append(hEffic)
    leg.Draw()
    c.Print("efficiency_%s.png" % (varName))
    c.Print("efficiency_%s.pdf" % (varName))

    objects.extend([c, leg])
