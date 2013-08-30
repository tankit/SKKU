import FWCore.ParameterSet.Config as cms

def lumiList( json ):
    import FWCore.PythonUtilities.LumiList as LumiList
    myLumis = LumiList.LumiList(filename = json ).getCMSSWString().split(',')
    return myLumis

def applyJSON( process, json ):
    # import PhysicsTools.PythonAnalysis.LumiList as LumiList
    # import FWCore.ParameterSet.Types as CfgTypes
    # myLumis = LumiList.LumiList(filename = json ).getCMSSWString().split(',')

    myLumis = lumiList( json )

    import FWCore.ParameterSet.Types as CfgTypes
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)

    # print process.source.lumisToProcess

