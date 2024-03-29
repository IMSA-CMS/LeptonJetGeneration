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
            
                        # 'SUSY:all = on',
                        # 'SLHA:readFrom = 2',
                        # 'SLHA:file = /uscms/home/hibarra/nobackup/softsusy/softsusy-4.1.12/inOutFiles/nmssmSLHAnoZ3Outputm0_300_m12_350',

                        #'print("There is no error here")', #could the error be here?
                        #'SLHA:file = leptonJet.spc',

                        # Turn off SUSY based production mechanism [12/16/2021]
                        # 'SUSY:all = on',
                        # 'SLHA:readFrom = 2',
                        # 'SLHA:file = sps1aWithDecays.spc',

                        # Sketchy Higgs production mechanism
                        # Two quarks (or gluons) make a Higgs that decays into dark fermions
                        # Dark fermions radiate and decay into dark photons
                        # Dark photons turn into lepton jets
                        # Do this all in PYTHIA if we can
                        # Shut off other Higgs decay modes and make this decay 100% for now

                        'HiggsSM:all = on',
                        '25:onMode = off',
                        '25:m0 = 1000.0',


                        # 5) Set Higgs mass, width and branching ratios by hand.
                        # Values for 125 GeV SM Higgs taken from
                        # Handbook of LHC Higgs Cross Sections: 1. Inclusive Observables,
                        # S. Dittmaier et al., CERN-2011-002 [arXiv:1101.0593].
                        # Note: sum is not quite unity, and so BR's are rescaled accordingly.

                        # see whether this will work -- screw with things
                        '25:addChannel = 1 0.0 0 4900001 -4900001',
                        # '25:0:bRatio  = 0.0',
                        # '25:1:bRatio  = 0.0',
                        # '25:2:bRatio  = 0.0',
                        # '25:3:bRatio  = 0.0',
                        # '25:4:bRatio  = 0.0',
                        # '25:5:bRatio  = 0.0',
                        # '25:6:bRatio  = 0.0', # H -> e+ e-
                        # '25:7:bRatio  = 0.0', # H -> mu+ mu-
                        # '25:8:bRatio  = 1.0',
                        # '25:9:bRatio  = 0.0',
                        # '25:10:bRatio = 0.0',
                        # '25:11:bRatio = 0.0',
                        # '25:12:bRatio = 0.0',
                        # '25:13:bRatio = 0.0',
                        '25:76:bRatio = 1.0', # dark photon
                        # '25:0:meMode  = 100', # set meMode = 100 so that
                        # '25:1:meMode  = 100', #branching ratios are not
                        # '25:2:meMode  = 100', # overwritten at initialization
                        # '25:3:meMode  = 100',
                        # '25:4:meMode  = 100',
                        # '25:5:meMode  = 100',
                        # '25:6:meMode  = 100',
                        # '25:7:meMode  = 100',
                        # '25:8:meMode  = 100',
                        # '25:9:meMode  = 100',
                        # '25:10:meMode = 100',
                        # '25:11:meMode = 100',
                        # '25:12:meMode = 100',
                        # '25:13:meMode = 100',
                        '25:76:meMode = 100',
                        '25:onIfAny = 4900001',
                        '25:mMin = 1000.',
                        '25:mWidth = 10.',
                        # '25:tau0 = 1.97327e-13',

                        # 'HiggsSM:all = on',
                        # 'HiggsSM:gg2H = on',

                        # Is Higgs 25? 35? idk it shows multiple H_0
                        # '25:onMode = off',
                        # Decay to dark fermion only
                        # '25:oneChannel = 1 0.5 0 15 -15',
                        # '25:addChannel = 1 0.5 0 4900022 4900022',
                        # '25:onIfAny = 15',

                        
                        # meMode? 0
                        '4900001:oneChannel = 1 1.0 100 4900022 4900004',
                        '4900003:oneChannel = 1 1.0 100 4900022 4900022',

                        # '4900022:oneChannel = 1 1.0 100 11 -11',
                        # '4900022:addChannel = 1 0.5 100 13 -13',

                        # Is this how you force?
                        
                        # Do this for dark fermion
                        # '4900001:onMode = off',
                        # '4900001:onIfAll = 4900022 4900002',

                        # '4900003:onMode = off',
                        # '4900003:onIfAll = 4900022 4900022',

                        # '4900022:onMode = off',
                        # '4900022:onIfAll = 11 -11',
                        # '4900022:onIfAll = 13 -13',


                        # Dark Fermion is 4900001
                        '4900001:m0 = 50.0',

                        '4900002:m0 = 0.1',
                        '4900003:m0 = 4.5',
                        '4900004:m0 = 2.0',

                        '4900001:spinType = 2', # assume is fermion, changed from 1
                        '4900002:spinType = 1',
                        '4900003:spinType = 1',
                        '4900004:spinType = 2', # assume is HLSP

                        '4900001:chargeType = 0', # change from 0?
                        '4900002:chargeType = 0',
                        '4900003:chargeType = 0',
                        '4900004:chargeType = 0', # wouldn't dark be charged, change from 0

                        '4900001:colType = 0',
                        '4900002:colType = 0',
                        '4900003:colType = 0',
                        '4900004:colType = 0',

                        '4900001:name = darkPseudoScalar',
                        '4900002:name = darkLightHiggs',
                        '4900003:name = darkHeavyHiggs',
                        '4900004:name = darkFermion',

                        '4900001:antiname = darkPseudoScalarBar',
                        '4900002:antiname = darkLightHiggsBar',
                        '4900003:antiname = darkHeavyHiggsBar',
                        '4900004:antiname = darkFermionBar',

                        '4900005:m0 = 1000000',
                        '4900006:m0 = 1000000',
                        '4900011:m0 = 1000000',
                        '4900012:m0 = 1000000',
                        '4900013:m0 = 1000000',
                        '4900014:m0 = 1000000',
                        '4900015:m0 = 1000000',
                        '4900016:m0 = 1000000',
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
                        '4900001:mayDecay = true', # turn off for now
                        '4900002:mayDecay = off',
                        '4900004:mayDecay = off',
                        '4900022:m0 = 0.4',
                        # Just test because idk how it interprets this
                        '4900022:mWidth = 1.0',
                        '4900022:mMin = 5.0',
                        '4900022:mMax = 0.0',
                        '4900022:isResonance = true',
                        '4900022:tau0 = 1.97327e-13',
                        # Don't know what this means? 1.000e+00  5.000e+01  0.000e+00  1.97327e-13
                        # '4900022:0:meMode = 0',
                        '4900004:isResonance = false', # remove resonance for dark fermion
                        '4900004:mWidth = 0.0',
                        '4900004:mMin = 0.0',
                        '4900004:mMax = 0.0',
                        '4900004:tau0 = 1.97327e-3',


                        # I'm not sure what this means
                        'HiddenValley:Ngauge = 1',
                        'HiddenValley:doKinMix = on',
                        'HiddenValley:FSR = on',
                        'HiddenValley:alphaFSR = 0',

			'4900022:oneChannel = 1 0.526471 100 11 -11',
			'4900022:addChannel = 1 0.466296 100 13 -13',
		
                ),
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