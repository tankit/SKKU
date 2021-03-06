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
    "TagProbeFitTree_Run173406.root"
    ),
    InputDirectoryName = cms.string("muonEffs"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("Efficiency.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "75.0", "115.0", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "20", "500", "GeV/c"),
        #eta = cms.vstring("Probe #eta", "-2.5", "2.5", "")
        eta = cms.vstring("Probe #eta", "-1.8", "1.8", "")
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        #PassingMediumMuons = cms.vstring("", "dummy[pass=1,fail=0]")
        PassingMediumMuons = cms.vstring("PassingMediumMuons", "dummy[pass=1,fail=0]")
        #passingMediumMuonsNoRPC = cms.vstring("passingMediumMuonsNoRPC", "dummy[pass=1,fail=0]")
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
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
            "Voigtian::signalPass(mass, meanP[90,80,100], width[2.495], sigmaP[3,1,20])", ## allow different means and sigmas for
            "Voigtian::signalFail(mass, meanF[90,80,100], width[2.495], sigmaF[3,1,20])",    ## passing and failing probes
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
            ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory
        pt = cms.PSet(
            #specifies the efficiency of which category and state to measure 
            EfficiencyCategoryAndState = cms.vstring("PassingMediumMuons","pass"),
            #EfficiencyCategoryAndState = cms.vstring("passingMediumMuonsNoRPC","pass"),

            #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
            UnbinnedVariables = cms.vstring("mass"),
            #specifies the binning of parameters
            BinnedVariables = cms.PSet(
                pt = cms.vdouble(20, 30, 40, 50, 60, 100, 250)
            ),
            #first string is the default followed by binRegExp - PDFname pairs
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
        eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("PassingMediumMuons","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                #eta = cms.vdouble(-1.8, -1.5, -1.2, -0.9, -0.6, -0.3,
                #                  0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8)
                #eta = cms.vdouble(-1.8, -1.6, -1.4, -1.2, -1.0, -0.8,
                #                  -0.6, -0.4, -0.2, 
                #                  0.0, 0.2, 0.4, 0.6, 0.8,
                #                  1.0, 1.2, 1.4, 1.6, 1.8)
                eta = cms.vdouble(-1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8,
                                  -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                                  0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                                  0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8)
            ),
            BinToPDFmap = cms.vstring("twoVoigtians")
        ),
    )
)


process.fit = cms.Path(process.TagProbeFitTreeAnalyzer)
