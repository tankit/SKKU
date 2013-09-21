#!/usr/bin/env python

import sys, os
from array import array

if len(sys.argv) < 3:
    print "hist.py INPUTFILE.root OUTPUTFILE.root"
    sys.exit()

srcFileName = sys.argv[1]
outFileName = sys.argv[2]

from ROOT import *
from urllib import urlretrieve
if not os.path.exists('rootlogon.C'):
    urlretrieve('https://raw.github.com/cms-top-kr/tools/master/rootlogon.C', 'rootlogon.C')
gROOT.ProcessLine(".x rootlogon.C")
gROOT.ForceStyle()

baseCut = "abseta < 1.6 && pt > 20"
muonIdDefs = [
    ("tightMuons", "tightMuons"),
    ("looseRPCMuons", "looseRPCMuons"),
    ("mediumRPCMuons", "mediumRPCMuons"),
    ("tightRPCMuons", "tightRPCMuons"),
    ("looseRPCInclusive", "looseRPCMuons || tightMuons"),
    ("mediumRPCInclusive", "mediumRPCMuons || tightMuons"),
    ("tightRPCInclusive", "tightRPCMuons || tightMuons"),
]
histDefs = [
    ("AbsEta", "Pseudorapidity |#eta|", "abseta", [0.0, 0.9, 1.2, 1.4, 1.6]),
    ("Pt"    , "Transverse momentum p_{T} (GeV/c)", "pt", [20, 30, 40, 50, 60, 100, 250])
]

massMin, massMax, binWidth = 75, 105, 0.2
nbin = int((massMax-massMin)/binWidth)

srcFile = TFile(srcFileName)
tree = srcFile.Get("muonEffs/fitter_tree")
outFile = TFile(outFileName, "RECREATE")
for muonIdName, muonIdCut in muonIdDefs:
    print "Projecting", muonIdName
    muonIdDir = outFile.mkdir(muonIdName)
    for varName, varTitle, varExp, bins in histDefs:
        varDir = muonIdDir.mkdir(varName)
        varDir.cd()
        binsArr = array('d', bins)
        hFrame = TH1F("hFrame", "%s %s;%s" % (muonIdName, varName, varTitle), len(bins)-1, binsArr)
        hFrame.Write()
        for i in range(len(bins)-1):
            binCut = "%s >= %f && %s < %f" % (varExp, bins[i], varExp, bins[i+1])
            binDir = varDir.mkdir("bin_%d" % i)
            binDir.cd()

            cutStrPass = TNamed("cut_pass", "(%s) &&  (%s) && (%s)" % (baseCut, muonIdCut, binCut))
            cutStrFail = TNamed("cut_fail", "(%s) && !(%s) && (%s)" % (baseCut, muonIdCut, binCut))
            hM_pass = TH1F("hM_pass", "Passing candidates;Mass (GeV/c^{2});Entries per %fGeV/c^{2}" % binWidth, nbin, massMin, massMax)
            hM_fail = TH1F("hM_fail", "Failing candidates;Mass (GeV/c^{2});Entries per %fGeV/c^{2}" % binWidth, nbin, massMin, massMax)
            tree.Draw("mass>>hM_pass", cutStrPass.GetTitle(), "goff")
            tree.Draw("mass>>hM_fail", cutStrFail.GetTitle(), "goff")

            cutStrPass.Write()
            cutStrFail.Write()
            hM_pass.Write()
            hM_fail.Write()
