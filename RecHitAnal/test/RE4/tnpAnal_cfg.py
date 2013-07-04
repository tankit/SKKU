import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.tnpWithRE4 = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring("tnpTree_PRE_PO61.root"),
    InputDirectoryName = cms.string("muonEffs"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("tnpAnal_PRE_PO61.root"),
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "70.0", "110.0", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "20", "500", "GeV/c"),
        eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        phi = cms.vstring("Probe #phi", "-3.141592", "3.141592", ""),
        abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
    ),

    Categories = cms.PSet(
        PassingLooseMuons = cms.vstring("PassingLooseMuons", "dummy[pass=1,fail=0]"),
        PassingTightMuons = cms.vstring("PassingTightMuons", "dummy[pass=1,fail=0]"),
        PassingRPCMuons = cms.vstring("PassingRPCMuons", "dummy[pass=1,fail=0]"),
    ),

    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
            "Gaussian::signal(mass, mean[91.2, 89.0, 93.0], sigma[2.5, 0.5, 10.0])",
            "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
            "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
        gaussPlusQuadratic = cms.vstring(
            "Gaussian::signal(mass, mean[90,80,100], sigma[3,1,20])",
            "Chebychev::backgroundPass(mass, {cPass1[0,-5,5], cPass2[0,-5,5]})",
            "Chebychev::backgroundFail(mass, {cFail1[0,-5,5], cFail2[0,-5,5]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]",
        ),
        twoVoigtians = cms.vstring(
            "Voigtian::signalPass(mass, meanP[90,85,95], width[2.495], sigmaP[3,1,20])", ## allow different means and sigmas for
            "Voigtian::signalFail(mass, meanF[90,85,95], width, sigmaF[3,1,20])",    ## passing and failing probes
            #"Chebychev::backgroundPass(mass, {cPass1[0,-5,5], cPass2[0,-5,5]})",
            #"Chebychev::backgroundFail(mass, {cFail1[0,-5,5], cFail2[0,-5,5]})",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    Efficiencies = cms.PSet(
        loose__pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingLooseMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 500)
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        loose__eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingLooseMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-2.5, -1.85, -1.15, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.15, 1.85, 2.5),
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        loose_re4__pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingLooseMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 500),
                abseta = cms.vdouble( 1.15,1.85 )
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        tight__pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingTightMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 500)
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        tight__eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingTightMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-2.5, -1.85, -1.15, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.15, 1.85, 2.5),
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        tight_re4__pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingTightMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 500),
                abseta = cms.vdouble( 1.15,1.85 )
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        rpcmu__pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingRPCMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 500)
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        rpcmu__eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingRPCMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                eta = cms.vdouble(-2.5, -1.85, -1.15, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.15, 1.85, 2.5),
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        rpcmu_re4__pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingRPCMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 500),
                abseta = cms.vdouble( 1.15,1.85 )
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
    ),
)

process.tnpWithoutRE4 = process.tnpWithRE4.clone(
    InputFileNames = cms.vstring("tnpTree_START61.root"),
    OutputFileName = cms.string("tnpAnal_START61.root"),
)

process.fit = cms.Path(process.tnpWithRE4 + process.tnpWithoutRE4)
