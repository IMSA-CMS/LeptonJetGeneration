print("You just lost the game")
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
from mcCmdLineOptions_cfi import registerDefaultMCOptions
registerDefaultMCOptions(options)
options.register ('zPrimeModel',
                  "zPrimeSSM",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "which Z' model to use")
options.register ('interferenceMode',
                  3,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          
                  "Z/gamma/Z' interference setting")
options.parseArguments()

print("Please let process be defined")
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')
print("I need to study for MVC")
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1000),
    limit = cms.untracked.int32(10000000)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)
process.source = cms.Source("EmptySource",
			    firstLuminosityBlock = cms.untracked.uint32(options.seed),
			    numberEventsInLuminosityBlock = cms.untracked.uint32(100)
)

process.RandomNumberGeneratorService.generator.initialSeed = options.seed*10000

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6','')
isDY = "on"
isCImumu = "off"
isCIee = "off"

if options.ciGen:
	isDY = "off"
	if options.pdgId == 11:
		isCIee = "on"
		isCImumu = "off"
	else:
		isCIee = "off"
		isCImumu = "on"


process.generator = cms.EDFilter("Pythia8GeneratorFilter",
        comEnergy =  cms.double(options.comEnergy*1000),
        crossSection = cms.untracked.double(1e10),
        filterEfficiency = cms.untracked.double(1),
        maxEventsToPrint = cms.untracked.int32(0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        PythiaParameters = cms.PSet(
                processParameters = cms.vstring(
                        'Main:timesAllowErrors    = 10000',
                        # 'ParticleDecays:limitTau0 = on',
                        # 'ParticleDecays:tauMax = 10',
                        # 'Tune:pp 5',
                        # 'Tune:ee 3',
 			# #'WeakSingleBoson:ffbar2ffbar(s:gmZ)= %s'%(isDY),
                        # #'ContactInteractions:QCffbar2eebar = %s'%(isCIee),
                        # #'ContactInteractions:QCffbar2mumubar = %s'%(isCImumu),
                        # '23:onMode = off',
                        # '23:onIfAny = '+str(options.pdgId),
			# #'PartonLevel:MPI = '+str(options.ULE),
			# #'PartonLevel:ISR = '+str(options.ISR),
			# #'PartonLevel:FSR = '+str(options.FSR),
                        # 'PartonLevel:MPI = off',
                        # 'PartonLevel:ISR = off',
                        # 'PartonLevel:FSR = off',
                        # # turn on for final
			# 'PhaseSpace:pTHatMin = '+str(options.pTMin),
                        # 'PhaseSpace:mHatMin = '+str(options.minMass),
                        # 'PhaseSpace:mHatMax = '+str(options.maxMass),
			# 'PhaseSpace:pTHatMax = '+str(options.pTMax),
			# #'ContactInteractions:etaLL = '+str(options.helicityLL),
			# #'ContactInteractions:etaLR = '+str(options.helicityLR),
			# #'ContactInteractions:etaRL = '+str(options.helicityRL),
			# #'ContactInteractions:etaRR = '+str(options.helicityRR),
			# #'ContactInteractions:Lambda = '+str(options.Lambda),
                        # 'SUSY:all = on',
                        # #'HiddenValley:all = on',
                        # #'SUSY:qqbar2squarkantisquark = on',
                        # #'SUSY:idA = 1000006',
                        # 'SLHA:readFrom = 2',
		
                        # '4900022:mayDecay = true',
                        # '1000022:mayDecay = true',
                        # '4900002:mayDecay = off',
                        # '4900004:mayDecay = off',
                        # '4900022:m0 = 1',    
                        # '4900022:0:meMode = 0', 		
		
                        # #'SLHA:minDecayDeltaM = 0.001',
		
                        # 'HiddenValley:Ngauge = 1',
                        # 'HiddenValley:doKinMix = on',
                        # 'HiddenValley:FSR = on',
                        # 'HiddenValley:alphaFSR = 0.3'
                         'Tune:pp 5',
            'Tune:ee 3',
            'PartonLevel:MPI = on',
            'PartonLevel:ISR = on',
            'PartonLevel:FSR = on',
            'SUSY:all = on',
                        'SLHA:readFrom = 2',
                        'SLHA:file = /uscms/home/hibarra/nobackup/softsusy/softsusy-4.1.12/inOutFiles/nmssmSLHAnoZ3Outputm0_300_m12_350',
                        #'print("There is no error here")', #could the error be here?
                        '4900021:m0 = 1000000',
                        '4900023:m0 = 1000000',
                        '4900101:m0 = 1000000',
                        '4900111:m0 = 1000000',
                        '4900113:m0 = 1000000',
                        '4900211:m0 = 1000000',
                        '4900213:m0 = 1000000',
                        '4900991:m0 = 1000000',
###Change ends here
                        '4900022:mayDecay = true',
                        '1000022:mayDecay = true',
                        '4900002:mayDecay = off',
                        '4900004:mayDecay = off',
                        '4900022:m0 = 0.3',
                        '4900022:0:meMode = 0',

                        'HiddenValley:Ngauge = 1',
                        'HiddenValley:doKinMix = on',
                        'HiddenValley:FSR = on',
                        'HiddenValley:alphaFSR = 30'),
                        #'print("there is no error here either")'), ##changing this from 0 to 30 didn't work
                #print ("alphaFSR error not here")
                parameterSets = cms.vstring('processParameters')
        )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)


process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
                                        SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('file:output.root'),
    outputCommands = process.AODSIMEventContent.outputCommands
)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)#*process.crossSecTreeMaker*process.pdfTreeMaker)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.AODSIMoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq


###-- Dump config -----------------------------------------------
file = open('GenTest.txt','w')
file.write(str(process.dumpPython()))
file.close()