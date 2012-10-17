#!/usr/bin/env python

from ROOT import *
import os

gROOT.ProcessLine(".x rootlogon.C")
gObjects = []
gPrefix = "20121011"
gCategory = "pt45_nj4_nb2"
#gStackColors = [kGreen+1, kAzure-2, kBlue-2, kMagenta, kYellow, kRed]
gStackColors = [kGreen+1, kMagenta, kRed, kYellow, kBlue-2, kAzure-2]

def styleLegend(leg):
    leg.SetFillColor(kWhite)
    leg.SetBorderSize(1)
    return leg

def stackPlot(file, name, paths):
    l = styleLegend(TLegend(0.7, 0.65, 0.9, 0.9))
    hs = []
    for path in paths:
        h = file.Get(path)
        if not h: continue
        h = h.Clone()
        i = len(hs)
        h.SetFillColor(gStackColors[i])
        l.AddEntry(h, h.GetTitle(), "f")
        hs.append(h)

    hStack = THStack(name, "%s;%s;%s" % (name, hs[0].GetXaxis().GetTitle(), hs[0].GetYaxis().GetTitle()))
    hs.reverse()
    for h in hs:
        hStack.Add(h)

    #c = TCanvas("c_"+name, name, 500, 500)
    #hStack.Draw()
    #l.Draw()

    if "SEvt" in name : dir = "SEvt"
    else: dir = "BEvt"
    #hh = file.Get(dir+"/hM_JJ").Clone()
    #hh.Sumw2()
    #hh.Draw("same")
    #gObjects.append(hh)

    #c.Print("%s_%s.png" % (gPrefix, name))

    sum = 0
    for h in hs: sum += h.Integral()
    for h in hs:
        print "^ %s | %.1f(%.1f%%) |" % (h.GetTitle(), h.Integral(), 100*h.Integral()/sum)

    #gObjects.extend([c, hStack, l])
    gObjects.extend([hStack, l])

    return (hStack, l)

def cmpPlot(file, plot1, plot2):
    title1, path1 = plot1
    title2, path2 = plot2

    h1 = file.Get(path1)
    h2 = file.Get(path2)

    if not h1 or not h2: return

    h1 = h1.Clone()
    h2 = h2.Clone()

    name = os.path.basename(path1)
    c = TCanvas("c_"+name, name, 500, 500)
    c.SetLogy()

    h1.SetLineColor(kBlue)
    h2.SetLineColor(kRed)

    l = styleLegend(TLegend(0.7, 0.75, 0.9, 0.9))
    l.AddEntry(h1, title1, "l")
    l.AddEntry(h2, title2, "l")

    h1.GetYaxis().SetTitle("Normalized")
    h2.GetYaxis().SetTitle("Normalized")

    yMax1 = h1.GetMaximum()
    yMax2 = h2.GetMaximum()
    if h1.Integral() > 0 : yMax1 /= h1.Integral()
    if h2.Integral() > 0 : yMax2 /= h2.Integral()
    if yMax2 > yMax1:
        h1, h2 = h2, h1
        yMax1, yMax2 = yMax2, yMax1

    if yMax1 > 0:
        h1.DrawNormalized()
        if yMax2 > 0: h2.DrawNormalized("same")
        else: h2.Draw("same")
    else:
        h1.Draw()
        
    l.Draw()

    #c.Print("%s_cmp_%s.png" % (gPrefix, name))

    gObjects.extend([c, l, h1, h2])

def getIntegral(hStack):
    sum = 0
    for h in hStack.GetHists():
        sum += h.Integral()
    return sum

if __name__ == "__main__":
    f = TFile.Open("hist.root")
    cmpPlot(f, ["Same evet", "SEvt/hMMC_JJ"], ["Bi event", "BEvt/hMMC_JJ"])
    cmpPlot(f, ["Same evet", "SEvt/hMMC_JK"], ["Bi event", "BEvt/hMMC_JK"])
    cmpPlot(f, ["Same evet", "SEvt/hMMC_KK"], ["Bi event", "BEvt/hMMC_KK"])
    cmpPlot(f, ["Same evet", "SEvt/hMMC_HBJ"], ["Bi event", "BEvt/hMMC_HBJ"])
    cmpPlot(f, ["Same evet", "SEvt/hMMC_LBJ"], ["Bi event", "BEvt/hMMC_LBJ"])
    cmpPlot(f, ["Same evet", "SEvt/hMMC_BX"], ["Bi event", "BEvt/hMMC_BX"])

    #stackPlot(f, "SEvt", ["SEvt/hMMC_JJ", "SEvt/hMMC_JK", "SEvt/hMMC_KK", "SEvt/hMMC_BJ", "SEvt/hMMC_BX"])
    #stackPlot(f, "BEvt", ["BEvt/hMMC_JJ", "BEvt/hMMC_JK", "BEvt/hMMC_KK", "BEvt/hMMC_BJ", "BEvt/hMMC_BX"])

    fOriginal = TFile.Open("original/best_TTJets_mass172_5_MuHad_loosePatJetsPF_TriCentral_looseJetCut_use2Jets.root")
    hSEvtOrig_M_JJ = fOriginal.Get("%s/hMassW_same" % gCategory)
    hBEvtOrig_M_JJ = fOriginal.Get("%s/hMassW_bi" % gCategory)
    hSEvtOrig_M_JJ.SetMarkerStyle(20)
    hBEvtOrig_M_JJ.SetMarkerStyle(20)
    hSEvtOrig_M_JJ.SetMarkerSize(1)
    hBEvtOrig_M_JJ.SetMarkerSize(1)
    hSEvtOrig_M_JJ.SetTitle("Same event;Dijet mass (GeV/c^{2});")
    hBEvtOrig_M_JJ.SetTitle("Bi event;Dijet mass (GeV/c^{2});")
    #hSEvt.SetTitle(";;")
    #hBEvt.SetTitle(";;")

    hSEvt, legSEvt = stackPlot(f, "SEvt", ["SEvt/hMMC_JJ", "SEvt/hMMC_HBJ", "SEvt/hMMC_LBJ", "SEvt/hMMC_BX", "SEvt/hMMC_KK", "SEvt/hMMC_JK"])
    hSEvtOrig_M_JJ.Scale(getIntegral(hSEvt)/hBEvtOrig_M_JJ.Integral())
    cSEvtCmp = TCanvas("cSEvtCmp", "cSEvtCmp", 500, 500)
    hSEvtOrig_M_JJ.SetMaximum(hSEvtOrig_M_JJ.GetMaximum()*1.2)
    hSEvtOrig_M_JJ.SetMinimum(0)
    hSEvtOrig_M_JJ.Draw("")
    hSEvt.Draw("same")
    hSEvtOrig_M_JJ.Draw("same")
    legSEvt.AddEntry(hSEvtOrig_M_JJ, "Original", "p")
    legSEvt.Draw()
    cSEvtCmp.Print("%s_sevt_%s.png" % (gPrefix, gCategory))

    hBEvt, legBEvt = stackPlot(f, "BEvt", ["BEvt/hMMC_JJ", "BEvt/hMMC_HBJ", "BEvt/hMMC_LBJ", "BEvt/hMMC_BX", "BEvt/hMMC_KK", "BEvt/hMMC_JK"])
    hBEvtOrig_M_JJ.Scale(getIntegral(hBEvt)/hBEvtOrig_M_JJ.Integral())
    cBEvtCmp = TCanvas("cBEvtCmp", "cBEvtCmp", 500, 500)
    hBEvtOrig_M_JJ.SetMaximum(hBEvtOrig_M_JJ.GetMaximum()*1.2)
    hBEvtOrig_M_JJ.SetMinimum(0)
    hBEvtOrig_M_JJ.Draw("")
    hBEvt.Draw("same")
    hBEvtOrig_M_JJ.Draw("same")
    legBEvt.AddEntry(hBEvtOrig_M_JJ, "Original", "p")
    legBEvt.Draw()
    cBEvtCmp.Print("%s_bevt_%s.png" % (gPrefix, gCategory))
