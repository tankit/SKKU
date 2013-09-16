import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(
        'tnpTree.root',
    ),
    InputDirectoryName = cms.string("muonEffs"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("Efficiency.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(False),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "20", "500", "GeV/c"),
        #eta = cms.vstring("Probe #eta", "-1.8", "1.8", ""),
        abseta = cms.vstring("Probe #eta", "0", "1.8", ""),
        #phi = cms.vstring("Probe #phi", "-3.14", "3.14", "radian"),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        #PassingMediumMuons = cms.vstring("", "dummy[pass=1,fail=0]")
        #PassingMediumTightMuons = cms.vstring("PassingMediumTightMuons", "dummy[pass=1,fail=0]")
        #PassingMediumTightMuonsRPCMu = cms.vstring("PassingMediumTightMuonsRPCMu", "dummy[pass=1,fail=0]")
        #PassingTightMuons = cms.vstring("PassingTightMuons", "dummy[pass=1,fail=0]")
        looseRPCMuons = cms.vstring("looseRPCMuons", "dummy[pass=1,fail=0]")
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        twoVoigtians = cms.vstring(
            "Voigtian::signalPass(mass, meanP[90,87,93], width[2.495], sigmaP[3,1,20])", ## allow different means and sigmas for
            "Voigtian::signalFail(mass, meanF[90,87,93], width[2.495], sigmaF[3,1,20])",    ## passing and failing probes
            "Chebychev::backgroundPass(mass, {cPass1[0,-5,5], cPass2[0,-5,5]})",
            "Chebychev::backgroundFail(mass, {cFail1[0,-5,5], cFail2[0,-5,5]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory
        pt = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("looseRPCMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 250)
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        abseta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("looseRPCMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                abseta = cms.vdouble(0.0, 0.2, 0.4, 0.6, 0.8,
                                     1.0, 1.2, 1.4, 1.6, 1.8)
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
    ),
    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(50),
)


process.fit = cms.Path(process.TagProbeFitTreeAnalyzer)
