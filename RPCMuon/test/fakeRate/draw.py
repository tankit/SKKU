#!/usr/bin/env python

from ROOT import *
from urllib import urlretrieve
if not os.path.exists('rootlogon.C'):
    urlretrieve('https://raw.github.com/cms-top-kr/tools/master/rootlogon.C', 'rootlogon.C')
gROOT.ProcessLine(".x rootlogon.C")
gROOT.ForceStyle()

#f_pi = TFile("pt_04to20/fit_Kshort.root")
#f_pi = TFile("pt_04to20/fit_Phi.root")
f_pi = TFile("fit_Kshort.root")
muonIds = ["RPCMuLoose", "RPCMuTight", "TMOneStationLoose", "TMOneStationTight", "TMTwoStationTest", "globalMuonTight"]
#muonIds = ["RPCMuLoose", "RPCMuTight", "TMOneStationLoose", "TMOneStationTight", "globalMuonTight"]# "TMTwoStationTest"]
colors = [kBlue, kAzure+1, kRed, kMagenta, kRed-1, kGreen+1, kBlack]
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
    table_pt_header = "| p_{T} bin |"
    table_pt_values = "| Fake rate (%) |"
    table_eta_header = "| abseta bin |"
    table_eta_values = "| Fake rate (%) |"

    grpPt  = f_pi.Get("%s/Pt/ratio" % muonId)
    grpEta = f_pi.Get("%s/AbsEta/ratio" % muonId)

    grpPt.SetLineColor(colors[muonIds.index(muonId)])
    grpPt.SetMarkerColor(colors[muonIds.index(muonId)])
    grpEta.SetLineColor(colors[muonIds.index(muonId)])
    grpEta.SetMarkerColor(colors[muonIds.index(muonId)])

    legMuonId.AddEntry(grpPt, muonId, "lp")

    cPt.cd()
    grpPt.Draw("p")

    for i in range(grpPt.GetN()):
        x = grpPt.GetX()[i]
        ex = grpPt.GetEX()[i]
        y = grpPt.GetY()[i]
        ey = grpPt.GetEY()[i]
        table_pt_header += " %g-%g |" % (x-ex, x+ex)
        table_pt_values += " %.2f+-%.2f |" % (y, ey)

    cEta.cd()
    grpEta.Draw("p")

    for i in range(grpEta.GetN()):
        x = grpEta.GetX()[i]
        ex = grpEta.GetEX()[i]
        y = grpEta.GetY()[i]
        ey = grpEta.GetEY()[i]
        table_eta_header += " %g-%g |" % (x-ex, x+ex)
        table_eta_values += " %.2f+-%.2f |" % (y, ey)

    print "=== %s ===" % muonId
    print table_pt_header
    print table_pt_values
    print table_eta_header
    print table_eta_values

cPt.cd()
legMuonId.Draw()
cEta.cd()
legMuonId.Draw()

cPt.Update()
cEta.Update()

cPt.Print("cFakerate_Pt.png")
cEta.Print("cFakerate_Eta.png")
cPt.Print("cFakerate_Pt.pdf")
cEta.Print("cFakerate_Eta.pdf")

