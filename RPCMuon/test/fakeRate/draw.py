#!/usr/bin/env python

from ROOT import *

gROOT.ProcessLine(".x rootlogon.C")

f_pi = TFile("fit_Kshort.root")
#muonIds = ["RPCMuLoose", "RPCMuMedium", "RPCMuTight", "TMOneStationLoose", "TMOneStationTight",]# "TMTwoStationTest"]
muonIds = ["RPCMuLoose", "RPCMuTight", "TMOneStationLoose", "TMOneStationTight", "globalMuonTight"]# "TMTwoStationTest"]
colors = [kBlack, kBlue, kRed, kMagenta, kGreen+2]
legMuonId = TLegend(0.20, 0.65, 0.50, 0.90)
legMuonId.SetBorderSize(0)
legMuonId.SetFillStyle(0)

hFramePt  = f_pi.Get("%s/Pt/hFrame" % muonIds[0])
hFrameEta = f_pi.Get("%s/AbsEta/hFrame" % muonIds[0])
hFramePt.SetMaximum(2)
hFrameEta.SetMaximum(2)
cPt  = TCanvas("cPt", "cPt", 500, 500)
hFramePt.Draw()
cEta = TCanvas("cEta", "cEta", 500, 500)
hFrameEta.Draw()
hists = []
for muonId in muonIds:
    grpPt  = f_pi.Get("%s/Pt/ratio" % muonId)
    grpEta = f_pi.Get("%s/AbsEta/ratio" % muonId)

    grpPt.SetLineColor(colors[muonIds.index(muonId)])
    grpPt.SetMarkerColor(colors[muonIds.index(muonId)])
    grpEta.SetLineColor(colors[muonIds.index(muonId)])
    grpEta.SetMarkerColor(colors[muonIds.index(muonId)])

    legMuonId.AddEntry(grpPt, muonId, "lp")

    cPt.cd()
    grpPt.Draw("p")

    print "*Pion fake rate of %s*" % muonId
    print "| p_{T} bin | Fake rate (%) |"
    for i in range(grpPt.GetN()):
        x = grpPt.GetX()[i]
        ex = grpPt.GetEX()[i]
        y = grpPt.GetY()[i]
        ey = grpPt.GetEY()[i]
        print "| %g-%g | %.2f +- %.2f |" % (x-ex, x+ex, y, ey)

    cEta.cd()
    grpEta.Draw("p")

    print "| abseta bin | Fake rate (%) |"
    for i in range(grpEta.GetN()):
        x = grpEta.GetX()[i]
        ex = grpEta.GetEX()[i]
        y = grpEta.GetY()[i]
        ey = grpEta.GetEY()[i]
        print "| %g-%g | %.2f +- %.2f |" % (x-ex, x+ex, y, ey)

cPt.cd()
legMuonId.Draw()
cEta.cd()
legMuonId.Draw()

cPt.Update()
cEta.Update()

cPt.Print("cFakerate_Pt.png")
cEta.Print("cFakerate_Eta.png")

## Pion fakerate in differential charge
f_pi1 = TFile("fit_Kshort_1.root")
f_pi2 = TFile("fit_Kshort_2.root")
muonId = "RPCMuLoose"

legPiCharge = TLegend(0.20, 0.65, 0.50, 0.90)
legPiCharge.SetBorderSize(0)
legPiCharge.SetFillStyle(0)

cPt_pi = TCanvas("cPt_pi", "cPt_pi", 500, 500)
hFramePt.Draw()
grpPt_pi  = f_pi.Get("%s/Pt/ratio" % muonId)
grpPt_pi1 = f_pi1.Get("%s/Pt/ratio" % muonId)
grpPt_pi2 = f_pi2.Get("%s/Pt/ratio" % muonId)
grpPt_pi.SetLineColor(kBlack)
grpPt_pi1.SetLineColor(kRed)
grpPt_pi2.SetLineColor(kBlue)
grpPt_pi.SetMarkerColor(kBlack)
grpPt_pi1.SetMarkerColor(kRed)
grpPt_pi2.SetMarkerColor(kBlue)
grpPt_pi.Draw("p")
grpPt_pi1.Draw("p")
grpPt_pi2.Draw("p")
legPiCharge.AddEntry(grpPt_pi1, "#pi^{+}", "lp")
legPiCharge.AddEntry(grpPt_pi, "#pi^{#pm}", "lp")
legPiCharge.AddEntry(grpPt_pi2, "#pi^{-}", "lp")
legPiCharge.Draw()

cEta_pi = TCanvas("cEta_pi", "cEta_pi", 500, 500)
hFrameEta.Draw()
grpEta_pi  = f_pi.Get("%s/AbsEta/ratio" % muonId)
grpEta_pi1 = f_pi1.Get("%s/AbsEta/ratio" % muonId)
grpEta_pi2 = f_pi2.Get("%s/AbsEta/ratio" % muonId)
grpEta_pi.SetLineColor(kBlack)
grpEta_pi1.SetLineColor(kRed)
grpEta_pi2.SetLineColor(kBlue)
grpEta_pi.SetMarkerColor(kBlack)
grpEta_pi1.SetMarkerColor(kRed)
grpEta_pi2.SetMarkerColor(kBlue)
grpEta_pi.Draw("p")
grpEta_pi1.Draw("p")
grpEta_pi2.Draw("p")
legPiCharge.Draw()

## Kaons
f_K  = TFile("fit_Phi.root")
f_K1 = TFile("fit_Phi_1.root")
f_K2 = TFile("fit_Phi_2.root")

legKCharge = TLegend(0.20, 0.65, 0.50, 0.90)
legKCharge.SetBorderSize(0)
legKCharge.SetFillStyle(0)

cPt_K = TCanvas("cPt_K", "cPt_K", 500, 500)
hFramePt.Draw()
grpPt_K  = f_K.Get("%s/Pt/ratio" % muonId)
grpPt_K1 = f_K1.Get("%s/Pt/ratio" % muonId)
grpPt_K2 = f_K2.Get("%s/Pt/ratio" % muonId)
grpPt_K.SetLineColor(kBlack)
grpPt_K1.SetLineColor(kRed)
grpPt_K2.SetLineColor(kBlue)
grpPt_K.SetMarkerColor(kBlack)
grpPt_K1.SetMarkerColor(kRed)
grpPt_K2.SetMarkerColor(kBlue)
grpPt_K.Draw("p")
grpPt_K1.Draw("p")
grpPt_K2.Draw("p")
legKCharge.AddEntry(grpPt_K1, "K^{+}"  , "lp")
legKCharge.AddEntry(grpPt_K , "K^{#pm}", "lp")
legKCharge.AddEntry(grpPt_K2, "K^{-}"  , "lp")
legKCharge.Draw()
cPt_K.Print("cPt_K.png")

cEta_K = TCanvas("cEta_K", "cEta_K", 500, 500)
hFrameEta.Draw()
grpEta_K  = f_K.Get("%s/AbsEta/ratio" % muonId)
grpEta_K1 = f_K1.Get("%s/AbsEta/ratio" % muonId)
grpEta_K2 = f_K2.Get("%s/AbsEta/ratio" % muonId)
grpEta_K.SetLineColor(kBlack)
grpEta_K1.SetLineColor(kRed)
grpEta_K2.SetLineColor(kBlue)
grpEta_K.SetMarkerColor(kBlack)
grpEta_K1.SetMarkerColor(kRed)
grpEta_K2.SetMarkerColor(kBlue)
grpEta_K.Draw("p")
grpEta_K1.Draw("p")
grpEta_K2.Draw("p")
legKCharge.Draw()
cEta_K.Print("cEta_K.png")

