# First try to reconstruct the channel B0 -> Ks Mu Mu
# work still in progress 

import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.BParkingNano.trgbits_cff import *

##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import *
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import muonTriggerMatchedTable, muonTriggerMatchedTables,  muonsBParkMCMatchForTable, selectedMuonsMCMatchEmbedded, muonBParkMCTable
from PhysicsTools.BParkingNano.tracksBPark_cff import tracksBParkMCMatchForTable, tracksBParkMCMatchEmbedded, tracksBParkMCTable

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
#from PhysicsTools.BParkingNano.tracksBPark_cff import *
from PhysicsTools.BParkingNano.common_cff import *

