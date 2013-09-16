#!/usr/bin/env python

import sys, os
from array import array
if len(sys.argv) < 4:
    print "hist.py : project fakerate tree to histograms"
    print "  Usage : hist.py MODENAME INPUTFILE PROBE_LEG_ID"
    print "hist_MODENAME_LEGID.root will be created."
    sys.exit()

modName = sys.argv[1]
srcFileName = sys.argv[2]

probeId = int(sys.argv[3])
tagId = (probeId+1)%2
outFileName = 'hist_%s_%d.root' % (modName, probeId)

modes = {
    "Kshort":(0.43, 0.56, 0.002),
    "Lambda":(1.08, 1.20, 0.002),
    "Phi"   :(1.00, 1.04, 0.001),
    "Jpsi"  :(2.80, 3.40, 0.010),
}

if modName not in modes:
    print "Mode not in modes. Choose among", modes.keys()

massMin, massMax, binWidth = modes[modName]
nbin = int(round((massMax-massMin)/binWidth))

from ROOT import *
if os.path.exists("rootlogon.C"):
    gROOT.ProcessLine(".x rootlogon.C")

muonTypes = {
    "RPCMuLoose" :"muonId%d_RPCMuLoose  == 1" % probeId,
    "RPCMuMedium":"muonId%d_RPCMuMedium == 1" % probeId,
    "RPCMuTight" :"muonId%d_RPCMuTight  == 1" % probeId,

    "TMOneStationLoose": "muonId%d_TMOneStationLoose == 1" % probeId,
    "TMOneStationTight": "muonId%d_TMOneStationTight == 1" % probeId,
    "TMTwoStationTest" : "muonId%d_TMTwoStationTest  == 1" % probeId,

    "globalMuon"      :"muonId%d_globalMuon       == 1" % probeId,
    "globalMuonTight" :"muonId%d_globalMuonTight  == 1" % probeId,
    "globalMuonMedium":"muonId%d_globalMuonMedium == 1" % probeId,
}

baseCut = "abs(track%d.eta()) < 1.6 && track%d.pt() > 4" % (probeId, probeId)
histDefs = [
    ("AbsEta", "Pseudorapidity |#eta|", "abs(track%d.eta())" % probeId, [0.0, 0.8, 1.2, 1.6]),
    ("Pt"    , "Transverse momentum p_{T} (GeV/c)", "track%d.pt()" % probeId, [4,5,7,10,20,30,50,500]),
]

srcFile = TFile(srcFileName)
tree = srcFile.Get("fake%s/tree" % modName)
outFile = TFile(outFileName, "RECREATE")
for muonType in muonTypes:
    print "Projecting", muonType
    muonIdCut = muonTypes[muonType]
    muonIdDir = outFile.mkdir(muonType)
    for varName, varTitle, varExp, bins in histDefs:
        varDir = muonIdDir.mkdir(varName)

        varDir.cd()
        binsArr = array('d', bins)
        hFrame = TH1F("hFrame", "%s %s;%s" % (muonType, varName, varTitle), len(bins)-1, binsArr)
        hFrame.Write()
        for i in range(len(bins)-1):
            binCut = "%s >= %f && %s < %f" % (varExp, bins[i], varExp, bins[i+1])
            binDir = varDir.mkdir("bin_%d" % i)
            binDir.cd()

            cutStrPass = TNamed("cut_pass", "(%s) &&  (%s) && (%s)" % (baseCut, muonIdCut, binCut))
            cutStrFail = TNamed("cut_fail", "(%s) && !(%s) && (%s)" % (baseCut, muonIdCut, binCut))
            hM_pass = TH1F("hM_pass", "Passing candidates;Mass (GeV/c^{2});Entries per 2MeV/c^{2}", nbin, massMin, massMax)
            hM_fail = TH1F("hM_fail", "Failing candidates;Mass (GeV/c^{2});Entries per 2MeV/c^{2}", nbin, massMin, massMax)
            tree.Draw("mass>>hM_pass", cutStrPass.GetTitle(), "goff")
            tree.Draw("mass>>hM_fail", cutStrFail.GetTitle(), "goff")

            cutStrPass.Write()
            cutStrFail.Write()
            hM_pass.Write()
            hM_fail.Write()
