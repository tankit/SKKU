#!/usr/bin/env python

import sys, os
from array import array
if len(sys.argv) < 3:
    print "hist.py INPUTFILE.root OUTPUTFILE.root"
    sys.exit(2)
srcFileName = sys.argv[1]
outFileName = sys.argv[2]

from ROOT import *
if os.path.exists("rootlogon.C"):
    gROOT.ProcessLine(".x rootlogon.C")

probeId = 1
tagId = (probeId+1)%2
modName = "fakeKshort"
massMin, massMax = 0.43, 0.56
nbin = int(round((massMax-massMin)/0.002))

muonTypes = {
    "RPCMuLoose" :"muonId%d_RPCMuLoose  == 1" % probeId,
    "RPCMuMedium":"muonId%d_RPCMuMedium == 1" % probeId,
    "RPCMuTight" :"muonId%d_RPCMuTight  == 1" % probeId,
}

baseCut = "abs(track%d.eta()) < 1.6 && track%d.pt() > 4" % (probeId, probeId)
histDefs = [
    ("AbsEta", "Pseudorapidity |#eta|", "abs(track%d.eta())" % probeId, [0.0, 0.8, 1.2, 1.6]),
    ("Pt"    , "Transverse momentum p_{T} (GeV/c)", "track%d.pt()" % probeId, [4,5,7,10,20,100]),
]

srcFile = TFile(srcFileName)
tree = srcFile.Get("%s/tree" % modName)
outFile = TFile(outFileName, "RECREATE")
for muonType in muonTypes:
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
            tree.Draw("vertexMass>>hM_pass", cutStrPass.GetTitle(), "goff")
            tree.Draw("vertexMass>>hM_fail", cutStrFail.GetTitle(), "goff")

            cutStrPass.Write()
            cutStrFail.Write()
            hM_pass.Write()
            hM_fail.Write()
