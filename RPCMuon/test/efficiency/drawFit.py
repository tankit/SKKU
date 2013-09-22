#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")
gROOT.ForceStyle()

imgPrefix = "20130922_efficiencyGain_"

varNames = ["Pt", "AbsEta"]
categories = {
#    "looseRPCMuons" :("RPCMuLoose" , kRed),
#    "mediumRPCMuons":("RPCMuMedium", kRed+1), 
#    "tightRPCMuons" :("RPCMuTight" , kMagenta),
    "tightMuons"     :("TightMuons" , kBlack),
    "looseRPCInclusive" :("Tight+RPCMuLoose" , kBlue),
    "mediumRPCInclusive":("Tight+RPCMuMedium", kAzure), 
    "tightRPCInclusive" :("Tight+RPCMuTight" , kGreen),
}

objs = []

summaryTableTex = ""
summaryTableTxt = ""

f = TFile("fit.root")
for varName in varNames:
    hFrame = None
    c = TCanvas("c%s" % varName, varName, 500, 500)
    leg = TLegend(0.6, 0.25, 0.92, 0.5)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for catName in categories:
        objDir = f.GetDirectory("%s/%s" % (catName, varName))
        if objDir == None: continue
        catTitle, color = categories[catName]

        if hFrame == None:
            hFrame = objDir.Get("hFrame")
            hFrame.SetMinimum(0)
            hFrame.SetMaximum(105)
            hFrame.Draw()
            summaryTableTex += "    %s " % varName
            summaryTableTxt += "| *%s* |" % varName
            for b in range(hFrame.GetNbinsX()):
                xLo = hFrame.GetXaxis().GetBinLowEdge(b+1)
                xHi = xLo + hFrame.GetXaxis().GetBinWidth(b+1)
                summaryTableTex += " & $%.2f-%.2f$" % (xLo, xHi)
                summaryTableTxt += " %.2f-%.2f |" % (xLo, xHi)
            summaryTableTex += "\\\\\n"
            summaryTableTxt += "| \n"
       
        grp = objDir.Get("efficiency")
        if grp == None: continue
        summaryTableTex += "    %s " % catTitle
        summaryTableTxt += "| *%s* |" % catTitle
        for b in range(grp.GetN()):
            y = grp.GetY()[b]
            ey = grp.GetEY()[b]
            summaryTableTex += " & $%.2f\pm%.2f$" % (y, ey)
            summaryTableTxt += " %.2f&plusmn;%.2f |" % (y, ey)
        summaryTableTxt += "\n"
        summaryTableTex += "\\\\\n"

        grp.SetLineColor(color)
        grp.SetMarkerColor(color)

        leg.AddEntry(grp, catTitle, "lp")
        grp.Draw("p")
        objs.append(grp)

    leg.Draw()

    objs.extend([hFrame, leg, c])

    c.Print("%s_%s.png" % (imgPrefix, varName))
    c.Print("%s_%s.pdf" % (imgPrefix, varName))

print "="*20
print summaryTableTex
print "="*20
print summaryTableTxt
print "="*20
