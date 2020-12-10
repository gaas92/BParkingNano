#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "CommonTools/CandUtils/interface/Booster.h"
#include "Math/GenVector/Boost.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include <vector>
#include <memory> 
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"
#include <Math/VectorUtil.h>

#include "TVector3.h"
#include "TMatrixD.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "DataFormats/Math/interface/deltaR.h"
 
//GAAS
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include <DataFormats/Common/interface/View.h>
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"



class BToKsMuMuBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToKsMuMuBuilder(const edm::ParameterSet &cfg):
    k_selection_{cfg.getParameter<std::string>("kaonSelection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    //HCL
    dzCut_{cfg.getParameter<double>("dzCut")},
    trgMuonToken_{consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("trgMuon"))},
    vertexSrc_{consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>( "vertexCollection" ) ) },
    //HCL
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kaons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kaons") )},
    kaons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kaonsTransientTracks") )},
    isotracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
    isolostTracksToken_(consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))),
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},

    //GAAS
    v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(cfg.getParameter<edm::InputTag>("secundaryVerticesPtr"))),	       
    tracksCollection_label(consumes<edm::View<pat::PackedCandidate>>(cfg.getParameter<edm::InputTag>("tracks"))),
    muons_Label(consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
    l1_selection_{cfg.getParameter<std::string>("lep1Selection")},
    l2_selection_{cfg.getParameter<std::string>("lep2Selection")},
    DLB_pre_vtx_selection_{cfg.getParameter<std::string>("DLB_preVtxSelection")},
    DLB_post_vtx_selection_{cfg.getParameter<std::string>("DLB_postVtxSelection")},

    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )} {
      produces<pat::CompositeCandidateCollection>("Bcollection");
      //HCL
      produces<nanoaod::FlatTable>("VertexTable");
      //HCL
    }

  ~BToKsMuMuBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::CompositeCandidate> k_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit

  //HCL
  const double dzCut_;
  const edm::EDGetTokenT<pat::MuonCollection> trgMuonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_; 
  //HCL

  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kaons_;
  const edm::EDGetTokenT<TransientTrackCollection> kaons_ttracks_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 

  //GAAS
  const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> v0PtrCollection_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> tracksCollection_label;
  const edm::EDGetTokenT<pat::MuonCollection> muons_Label;
  const StringCutObjectSelector<Lepton> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<Lepton> l2_selection_; // cut on sub-leading lepton
  const StringCutObjectSelector<pat::CompositeCandidate> DLB_pre_vtx_selection_;  // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> DLB_post_vtx_selection_; // cut on the di-lepton after the SV fit

  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
};

//void BToKsMuMuBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & ) const {
void BToKsMuMuBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const & iSetup) const {

  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);
  
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kaons;
  evt.getByToken(kaons_, kaons);
  
  edm::Handle<TransientTrackCollection> kaons_ttracks;
  evt.getByToken(kaons_ttracks_, kaons_ttracks);  

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  





  //HCL PRIMERO OBTENGAMOS LOS VERTICES QUE NOS INTERESAN
  edm::Handle<reco::VertexCollection> vertexHandle;
  evt.getByToken(vertexSrc_, vertexHandle);

  edm::Handle<pat::MuonCollection> trgMuons;
  evt.getByToken(trgMuonToken_, trgMuons);
  
  //GAAS
  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
  evt.getByToken(v0PtrCollection_,  theV0PtrHandle);
  //GAAS Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
  edm::Handle<edm::View<pat::PackedCandidate> > thePATTrackHandle;
  evt.getByToken(tracksCollection_label,thePATTrackHandle);
  edm::Handle<pat::MuonCollection> thePATMuonHandle;
  evt.getByToken(muons_Label ,thePATMuonHandle);

  std::vector<float> vx,vy,vz, dzTrgMu;// ,dlenSig,pAngle;

  const reco::Vertex & vertex = vertexHandle->front(); 
  const pat::Muon & trgmu = trgMuons->front();
  dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
  vx.push_back(vertex.position().x());
  vy.push_back(vertex.position().y());
  vz.push_back(vertex.position().z());  
  
  int nv = 0;
  for (const reco::Vertex & vertex : *vertexHandle){
    for (const pat::Muon & trgmu : *trgMuons){
      if(nv==3) break;
      if(fabs(vertex.position().z()-trgmu.vz())<dzCut_) {
        dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
        vx.push_back(vertex.position().x());
        vy.push_back(vertex.position().y());
        vz.push_back(vertex.position().z());
        nv++;
      }
    }
  }
  // if (dzTrgMu.size()==0){
  //      const reco::Vertex & vertex = vertexHandle->front(); 
  //      const pat::Muon & trgmu = trgMuons->front();
  //       dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
  //       vx.push_back(vertex.position().x());
  //       vy.push_back(vertex.position().y());
  //       vz.push_back(vertex.position().z());	
  // }
  // output
  auto pvTable = std::make_unique<nanoaod::FlatTable>(dzTrgMu.size(),"PV",false);

  pvTable->addColumn<float>("dzTrgMu",dzTrgMu,"abs(vz-trgMuon.vz)",nanoaod::FlatTable::FloatColumn,10);
  pvTable->addColumn<float>("vx",vx,"vx in cm",nanoaod::FlatTable::FloatColumn,10);
  pvTable->addColumn<float>("vy",vy,"vy in cm",nanoaod::FlatTable::FloatColumn,10);
  pvTable->addColumn<float>("vz",vz,"vz in cm",nanoaod::FlatTable::FloatColumn,10);
  //std::cout << "Pasamos el almacenamiento delos PVs\n\n" ;
  //HCL

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();

  std::vector<int> used_lep1_id, used_lep2_id, used_trk_id;

  //std::cout << "Kaons Size:  "<< kaons->size() << std::endl;
  std::cout << "Dileptons Size:  "<< dileptons->size() << std::endl;
  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  
  //GAAS my reco, taken from Jhovanny
  int passmu = 0;
  for(pat::MuonCollection::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {
    for(pat::MuonCollection::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) {
        if(iMuon1==iMuon2) continue;
	    //opposite charge 
	    if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------
  
	    reco::TrackRef glbTrackP;	  
	    reco::TrackRef glbTrackM;	  
	    
	    if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	    if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
	    
	    if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	    if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
	    
	    if( glbTrackP.isNull() || glbTrackM.isNull() ){
	        //std::cout << "continue due to no track ref" << endl;
	        continue;
	    }
        pat::CompositeCandidate new_cand;
        //no pT cuts in Dileptonbuilder neither MuonTriggerSelector
        //if(iMuon1->track()->pt()<4.0) continue; 
	    //if(iMuon2->track()->pt()<4.0) continue;
        //pT cuts are in config_BtoKmumu_cff
        //if (!(iMuon1->track()->pt()<1.4 || iMuon1->track()->pt()<1.4)) continue;  //at least one muon must have pT > 1.5 GeV
        float dRm1m2 = reco::deltaR(*iMuon1, *iMuon2);
        //if (dRm1m2 < 0.02) continue; // DiMuonBuilder cuts at 0.03

        //no highPurity Cuts in DiLeptonBuilder neither in MuonTriggerSelector
	    //if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	    //if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
	    //Let's check the vertex and mass
	    reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	    reco::TransientTrack muon2TT((*theB).build(glbTrackM));
        // *****  Trajectory states to calculate DCA for the 2 muons *********************
	    FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	    FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	    if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;
        // Measure distance between tracks at their closest approach
	    ClosestApproachInRPhi cApp;
	    cApp.calculate(mu1State, mu2State);
	    if( !cApp.status() ) continue;
	    float dca = fabs( cApp.distance() );
        //if(dca < 1.5) continue; // abs(dz1 -dz2) < 1. in DiMuonBuilder preVtxSelection

        
        // *****  end DCA for the 2 muons *********************

	    //The mass of a muon and the insignificant mass sigma 
	    //to avoid singularities in the covariance matrix.
	    ParticleMass muon_mass = 0.10565837; //pdg mass
	    ParticleMass psi_mass = 3.096916;
	    float muon_sigma = muon_mass*1.e-6;
	    //float psi_sigma = psi_mass*1.e-6;

	    //Creating a KinematicParticleFactory
	    KinematicParticleFactoryFromTransientTrack pFactory;
    
	    //initial chi2 and ndf before kinematic fits.
	    float chi = 0.;
	    float ndf = 0.;
	    std::vector<RefCountedKinematicParticle> muonParticles;
	    try {
	      muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	      muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	    }
	    catch(...) { 
	      std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	      continue;
	    }

        KinematicParticleVertexFitter fitter;   

	    RefCountedKinematicTree psiVertexFitTree;
	    try {
	        psiVertexFitTree = fitter.fit(muonParticles); 
	    }
	    catch (...) { 
	        std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	        continue;
	    }
  
	    if (!psiVertexFitTree->isValid()){
	        //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	        continue; 
	    }
  
	    psiVertexFitTree->movePointerToTheTop();
	    
	    RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	    RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	    
	    if( psi_vFit_vertex_noMC->chiSquared() < 0 ){
	        //std::cout << "negative chisq from psi fit" << endl;
	        //continue;
	    }

        double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	    //if(J_Prob_tmp<0.01) //Jhovanny
	    if(J_Prob_tmp<1.e-5) //DiMuonBuilder postVtxSelection
	    {
	       continue;
	    }
	  
	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared()>999) continue; // DiMuonBuilder cuts at 998 
	   //if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue; //no cut at all

	   //  ***************  
	   if ( theV0PtrHandle->size()>0 && thePATMuonHandle->size()>=2 ){
	        for ( std::vector<reco::VertexCompositePtrCandidate>::const_iterator iVee = theV0PtrHandle->begin();   iVee != theV0PtrHandle->end(); ++iVee ){
                //get Lam tracks from V0 candidate
		        std::vector<pat::PackedCandidate> v0daughters;
		        std::vector<reco::Track> theDaughterTracks;
		        v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(0))) );
		        v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(1))) );
		     		     
		        for(unsigned int j = 0; j < v0daughters.size(); ++j){
			        theDaughterTracks.push_back(v0daughters[j].pseudoTrack());
		        }
                new_cand.setP4(iMuon1->p4() + iMuon2->p4() + iVee->daughter(0)->p4() + iVee->daughter(1)->p4());
                new_cand.setCharge(iMuon1->charge() + iMuon2->charge() + iVee->daughter(0)->charge() + iVee->daughter(0)->charge());

                //Now let's see if these two tracks make a vertex
		        reco::TransientTrack pion1TT((*theB).build(theDaughterTracks[0]));
		        reco::TransientTrack pion2TT((*theB).build(theDaughterTracks[1]));		     
		     
		        ParticleMass pion_mass = 0.13957018;
		        ParticleMass Ks0_mass = 0.497614;
		        float pion_sigma = pion_mass*1.e-6;
		        float Ks0_sigma = Ks0_mass*1.e-6;
		     
		        //initial chi2 and ndf before kinematic fits.
		        float chi = 0.;
		        float ndf = 0.;
                std::vector<RefCountedKinematicParticle> pionParticles;
		        // vector<RefCountedKinematicParticle> muonParticles;
		        try {
		            pionParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
		            pionParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
		        }
		        catch(...) {
		            std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		            continue;
		        }
		        
		        RefCountedKinematicTree Ks0VertexFitTree;
		        try{
		            Ks0VertexFitTree = fitter.fit(pionParticles); 
		        }
		        catch(...) {
		            std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		            continue;
		        }
		        if (!Ks0VertexFitTree->isValid()) {
			        //std::cout << "invalid vertex from the Ks0 vertex fit" << std::endl;
			        continue; 
		        }
		        Ks0VertexFitTree->movePointerToTheTop(); 
                RefCountedKinematicParticle Ks0_vFit_noMC = Ks0VertexFitTree->currentParticle();
		        RefCountedKinematicVertex Ks0_vFit_vertex_noMC = Ks0VertexFitTree->currentDecayVertex();
    
		        if( Ks0_vFit_vertex_noMC->chiSquared() < 0 ){ 
			        //std::cout << "negative chisq from ks fit" << endl;
			        continue;
		        }
    
		        //some loose cuts go here
		        if(Ks0_vFit_vertex_noMC->chiSquared()>50) continue;
		        if(Ks0_vFit_noMC->currentState().mass()<0.45 || Ks0_vFit_noMC->currentState().mass()>0.55) continue;
    
		        Ks0VertexFitTree->movePointerToTheFirstChild();
		        RefCountedKinematicParticle T1CandMC = Ks0VertexFitTree->currentParticle();
    
		        Ks0VertexFitTree->movePointerToTheNextChild();
		        RefCountedKinematicParticle T2CandMC = Ks0VertexFitTree->currentParticle();
    
		        //  Ks0  mass constrain
		        // do mass constrained vertex fit
		        // creating the constraint with a small sigma to put in the resulting covariance 
		        // matrix in order to avoid singularities
		        // JPsi mass constraint is applied in the final B fit
    
		        KinematicParticleFitter csFitterKs;
		        KinematicConstraint * ks_c = new MassKinematicConstraint(Ks0_mass,Ks0_sigma); //remember to kill 
		        // add mass constraint to the ks0 fit to do a constrained fit:  
    
		        Ks0VertexFitTree = csFitterKs.fit(ks_c,Ks0VertexFitTree);
		        if (!Ks0VertexFitTree->isValid()){
		            //std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
		            continue; 
		        }
                Ks0VertexFitTree->movePointerToTheTop();
		        RefCountedKinematicParticle ks0_vFit_withMC = Ks0VertexFitTree->currentParticle();
		     
		        //Now we are ready to combine!
		     
		        std::vector<RefCountedKinematicParticle> vFitMCParticles;
		        vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		        vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		        vFitMCParticles.push_back(ks0_vFit_withMC);
                //no mass constrain 
			    KinematicParticleVertexFitter kcvFitter;
			    RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles);
		        if (!vertexFitTree->isValid()) {
		            //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		            continue;
		        }
                vertexFitTree->movePointerToTheTop();		     
		     
		        RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		        RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		        if (!bDecayVertexMC->vertexIsValid()){
		            //std::cout << "B MC fit vertex is not valid" << endl;
		            continue;
		        }
    
		        if(bCandMC->currentState().mass()<4.8 || bCandMC->currentState().mass()>6.0) continue;
    
		        //if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ){
			    //    //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
			    //    continue;
		        //}
		        //std::cout << "pass 461 continues ... "<< std::endl;
		        double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		        //if(B_Prob_tmp<0.01) //Jhovanny
		        if(B_Prob_tmp<0.001){  //Horacio hardcoded
			        continue;
		        }
                // get children from final B fit
		        vertexFitTree->movePointerToTheFirstChild();
		        RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
		        vertexFitTree->movePointerToTheNextChild();
		        RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
		        
		        vertexFitTree->movePointerToTheNextChild();
		        RefCountedKinematicParticle Ks0CandMC = vertexFitTree->currentParticle();
		        
		        KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		        KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
		        KinematicParameters psiMupKP;
		        KinematicParameters psiMumKP;
	            
		        if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
		        if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
		        if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
		        if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;

                GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
				                    mu1CandMC->currentState().globalMomentum().y(),
 				                    mu1CandMC->currentState().globalMomentum().z());

 		        GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
				                    mu2CandMC->currentState().globalMomentum().y(),
 				                    mu2CandMC->currentState().globalMomentum().z());

 		        GlobalVector Ks0p1vec(T1CandMC->currentState().globalMomentum().x(),
				                      T1CandMC->currentState().globalMomentum().y(),
 				                      T1CandMC->currentState().globalMomentum().z());

 		        GlobalVector Ks0p2vec(T2CandMC->currentState().globalMomentum().x(),
				                      T2CandMC->currentState().globalMomentum().y(),
					                  T2CandMC->currentState().globalMomentum().z());

                KinematicParameters Ks0Pi1KP = T1CandMC->currentState().kinematicParameters();
		        KinematicParameters Ks0Pi2KP = T2CandMC->currentState().kinematicParameters();
		        KinematicParameters Ks0PipKP;
		        KinematicParameters Ks0PimKP;

                if ( T1CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi1KP;
		        if ( T1CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi1KP;
		        if ( T2CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi2KP;
		        if ( T2CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi2KP;	 

		   // fill candidate variables now                      
            }// en V0 Tracks
        }// end if dimuon&& V0Tracks   
        passmu++;

    }    
  }// end muon loop   
  std::cout << "pass mu: "<< passmu <<std::endl;

  for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx) {
    edm::Ptr<pat::CompositeCandidate> k_ptr(kaons, k_idx);
    
   // std::cout << "K DCASig : " <<k_ptr->userFloat("DCASig") << std::endl;
   // std::cout << "K nValidHits : " <<k_ptr->userInt("nValidHits") << std::endl;
    
   // try {
    if( !k_selection_(*k_ptr) ) continue;
    
    math::PtEtaPhiMLorentzVector k_p4(
      k_ptr->pt(), 
      k_ptr->eta(),
      k_ptr->phi(),
      K_MASS
      );


    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      
     // std::cout << "Dimu id : " << ll_idx << std::endl;

      edm::Ptr<pat::CompositeCandidate> ll_prt(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_prt->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_prt->userCand("l2");
      int l1_idx = ll_prt->userInt("l1_idx");
      int l2_idx = ll_prt->userInt("l2_idx");
    
      //const pat::Muon muon1 = dynamic_cast<const pat::Muon>(ll_prt->userCand("l1"));
      //const pat::Muon muon2 = dynamic_cast<const pat::Muon>(ll_prt->userCand("l2"));

      pat::CompositeCandidate cand;
      cand.setP4(ll_prt->p4() + k_p4);
      cand.setCharge(ll_prt->charge() + k_ptr->charge());
      // Use UserCands as they should not use memory but keep the Ptr itself
      // Put the lepton passing the corresponding selection
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("K", k_ptr);
      cand.addUserCand("dilepton", ll_prt);

      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("k_idx", k_idx);
     
      auto dr_info = min_max_dr({l1_ptr, l2_ptr, k_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);

      float dr1 = reco::deltaR(*l1_ptr, *k_ptr);
      float dr2 = reco::deltaR(*l2_ptr, *k_ptr);
      cand.addUserFloat("dr_l1", dr1);
      cand.addUserFloat("dr_l2", dr2);
      cand.addUserFloat("dz_l1", fabs( l1_ptr->vz()-k_ptr->vz() ) );
      cand.addUserFloat("dz_l2", fabs( l2_ptr->vz()-k_ptr->vz() ) );

      // TODO add meaningful variables
      // Variables Pre fitting
      // muon1
      // cand.addUserFloat("l1_vx", l1_ptr->vx());
      // cand.addUserFloat("l1_vy", l1_ptr->vy());
      // cand.addUserFloat("l1_vz", l1_ptr->vz());
      // // muon2
      // cand.addUserFloat("l2_vx", l2_ptr->vx());
      // cand.addUserFloat("l2_vy", l2_ptr->vy());
      // cand.addUserFloat("l2_vz", l2_ptr->vz());
      // //kaon 
      // cand.addUserFloat("k_vx", k_ptr->vx());
      // cand.addUserFloat("k_vy", k_ptr->vy());
      // cand.addUserFloat("k_vz", k_ptr->vz());
      // //dimuon
      // cand.addUserFloat("dmuon_vx", ll_prt->vx());
      // cand.addUserFloat("dmuon_vy", ll_prt->vy());
      // cand.addUserFloat("dmuon_vz", ll_prt->vz());
      

      if( !pre_vtx_selection_(cand) ) continue;
    
      //HCL
      //std::cout << "BBB PreVertex Fitting    ";
    
      KinVtxFitter fitter(
        {leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx), kaons_ttracks->at(k_idx)},
        {l1_ptr->mass(), l2_ptr->mass(), K_MASS},
        {LEP_SIGMA, LEP_SIGMA, K_SIGMA} //some small sigma for the lepton mass
        );
      //HCL Se puede extraer el Vertex???
      //std::cout << "BBB PostVertex Fitting    \n";
      // RefCountedKinematicVertex myVertex = fitter.fitted_vtx_();
      // std::cout << myVertex->vertexIsValid << "\n\n";


      if(!fitter.success()) continue; // hardcoded, but do we need otherwise?
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );
      
     // std::cout << "----->  SUCCESS!!!\n";
      used_lep1_id.emplace_back(l1_idx);
      used_lep2_id.emplace_back(l2_idx);
      used_trk_id.emplace_back(k_idx);
      cand.addUserInt("sv_OK" , fitter.success()); 
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof()); // float??
      cand.addUserFloat("sv_prob", fitter.prob());
 
      cand.addUserFloat("fitted_mll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass());
      cand.addUserFloat("fitted_pt_ll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).pt());
      cand.addUserFloat("fitted_eta_ll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).eta());
      cand.addUserFloat("fitted_phi_ll" , (fitter.daughter_p4(0) + fitter.daughter_p4(1)).phi());
      
      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      // cand.addUserFloat("fitted_px"  , fitter.fitted_candidate().globalMomentum().x()); 
      // cand.addUserFloat("fitted_py"  , fitter.fitted_candidate().globalMomentum().y()); 
      // cand.addUserFloat("fitted_pz"  , fitter.fitted_candidate().globalMomentum().z()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fitter.fitted_candidate().mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)));      
      cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, cand.p4())
        );
      cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );
      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());

      // //Almacenemos la informacion del beamspot 
      // reco::BeamSpot bs = *beamspot;
      // cand.addUserFloat("beamSpot_x", bs.x(cand.vz()));
      // cand.addUserFloat("beamSpot_y", bs.y(cand.vz()));

      cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());
      cand.addUserFloat("vtx_ex" , sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      cand.addUserFloat("vtx_ey" , sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      cand.addUserFloat("vtx_ez" , sqrt(fitter.fitted_vtx_uncertainty().czz()));
      try{
        cand.addUserFloat("vtx_eyx", fitter.fitted_vtx_uncertainty().cyx());
        cand.addUserFloat("vtx_ezx", fitter.fitted_vtx_uncertainty().czx());
        cand.addUserFloat("vtx_ezy", fitter.fitted_vtx_uncertainty().czy());
      }
      catch(...){
        cand.addUserFloat("vtx_eyx", -1);
        cand.addUserFloat("vtx_ezx", -1);
        cand.addUserFloat("vtx_ezy", -1);
      }
      cand.addUserFloat("fitted_l1_pt" , fitter.daughter_p4(0).pt()); 
      cand.addUserFloat("fitted_l1_eta", fitter.daughter_p4(0).eta());
      cand.addUserFloat("fitted_l1_phi", fitter.daughter_p4(0).phi());
      cand.addUserFloat("l1_charge", l1_ptr->charge());

      cand.addUserFloat("fitted_l2_pt" , fitter.daughter_p4(1).pt()); 
      cand.addUserFloat("fitted_l2_eta", fitter.daughter_p4(1).eta());
      cand.addUserFloat("fitted_l2_phi", fitter.daughter_p4(1).phi());
      cand.addUserFloat("l2_charge", l2_ptr->charge());

      cand.addUserFloat("fitted_k_pt"  , fitter.daughter_p4(2).pt()); 
      cand.addUserFloat("fitted_k_eta" , fitter.daughter_p4(2).eta());
      cand.addUserFloat("fitted_k_phi" , fitter.daughter_p4(2).phi());
      cand.addUserFloat("k_charge", k_ptr->charge());
      

      if( !post_vtx_selection_(cand) ) continue;        

      //compute isolation
      float l1_iso03 = 0;
      float l1_iso04 = 0;
      float l2_iso03 = 0;
      float l2_iso04 = 0;
      float k_iso03  = 0;
      float k_iso04  = 0;
      float b_iso03  = 0;
      float b_iso04  = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon
        if (k_ptr->userCand("cand") ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id(), iTrk) || 
            track_to_lepton_match(l2_ptr, iso_tracks.id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk.eta(), trk.phi());
        float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk.eta(), trk.phi());
        float dr_to_k  = deltaR(cand.userFloat("fitted_k_eta") , cand.userFloat("fitted_k_phi") , trk.eta(), trk.phi());
        float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_k < 0.4){
          k_iso04 += trk.pt();
          if (dr_to_k < 0.3) k_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      } 


      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
      cand.addUserFloat("k_iso03" , k_iso03 );
      cand.addUserFloat("k_iso04" , k_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );






      // Aqui creemos el boost al CM del dilepton
      math::XYZTLorentzVector dilep = ll_prt->p4();
      ROOT::Math::Boost cmboost(dilep.BoostToCM());

      math::XYZTLorentzVector kaonCM(  cmboost( k_ptr->p4() )  );
      math::XYZTLorentzVector muonCM1, muonCM2;

      //where the thetal is the angle between the
      //direction of the m-(m+) lepton and the K+(K-)
      if (l1_ptr->charge()==k_ptr->charge()){
        muonCM1 = cmboost(fitter.daughter_p4(1)) ;
        muonCM2 = cmboost(fitter.daughter_p4(0)) ;
      }
      else {
        muonCM1 = cmboost(fitter.daughter_p4(0)) ;
        muonCM2 = cmboost(fitter.daughter_p4(1)) ;
              } 

      float costhetaL = ( muonCM1.x()*muonCM2.x() 
                         + muonCM1.y()*muonCM2.y() 
                         + muonCM1.z()*muonCM2.z() ) / (muonCM1.P()*muonCM2.P() );

      float costhetaKL = ( muonCM1.x()*kaonCM.x()
                         + muonCM1.y()*kaonCM.y()
                         + muonCM1.z()*kaonCM.z() ) / (muonCM1.P()*kaonCM.P() );

      cand.addUserFloat("cosTheta_mm", costhetaL);
      cand.addUserFloat("cosTheta_km", costhetaKL);


      std::vector<float> cosAlpha(4, -2);
      std::vector<float> lxy_pv(4,-1);
      std::vector<float> errP(4,-1);

      // CHEQUEMOS QUE EL CANDIDATO A "B" SATIFAGA LOS CORTES EN COS(ALPHA) Y SIGNIFICANCIA
      for( unsigned int iPV=0; iPV<dzTrgMu.size(); ++iPV ){ 
        //Debemos calcular el vector que une los vertices primario y el ajustado
        math::XYZVector vPS(cand.vx()-vx.at(iPV), cand.vy()-vy.at(iPV), cand.vz()-vz.at(iPV));
        //Momento espacial del candidato
        math::XYZVector Bp(fitter.fitted_candidate().globalMomentum().x(), fitter.fitted_candidate().globalMomentum().y(), fitter.fitted_candidate().globalMomentum().z());
        //CosALPHA
        cosAlpha[iPV] = vPS.Dot(Bp)/(vPS.R()*Bp.R());
        
        //Para significancia:  
        GlobalError err = fitter.fitted_vtx_uncertainty();
        GlobalPoint delta(cand.vx()-vx.at(iPV), cand.vy()-vy.at(iPV), 0.);  

        lxy_pv[iPV] = delta.perp();
        errP[iPV] = sqrt(err.rerr(delta));

      }

      cand.addUserFloat("cosAlpha0", cosAlpha[0]);
      cand.addUserFloat("cosAlpha1", cosAlpha[1]);
      cand.addUserFloat("cosAlpha2", cosAlpha[2]);
      cand.addUserFloat("cosAlpha3", cosAlpha[3]);

      cand.addUserFloat("lxy_pv0", lxy_pv[0]);
      cand.addUserFloat("lxy_pv1", lxy_pv[1]);
      cand.addUserFloat("lxy_pv2", lxy_pv[2]);
      cand.addUserFloat("lxy_pv3", lxy_pv[3]);

      cand.addUserFloat("significance0", lxy_pv[0]/errP[0]);
      cand.addUserFloat("significance1", lxy_pv[1]/errP[1]);
      cand.addUserFloat("significance2", lxy_pv[2]/errP[2]);
      cand.addUserFloat("significance3", lxy_pv[3]/errP[3]);
      

      TVector3 pv, sv, pT;
      pT.SetXYZ(fit_p4.px(),fit_p4.py(),0.0);
      pv.SetXYZ(vx.at(0),vy.at(0),vz.at(0));
      sv.SetXYZ(cand.vx(),cand.vy(),cand.vz());


      TMatrix ESV(3,3);
      TMatrix EPV(3,3);

      ESV(0,0) = fitter.fitted_vtx_uncertainty().cxx();
      ESV(1,1) = fitter.fitted_vtx_uncertainty().cyy();
      ESV(2,2) = fitter.fitted_vtx_uncertainty().czz();
      ESV(0,1) = fitter.fitted_vtx_uncertainty().cyx();
      ESV(0,2) = fitter.fitted_vtx_uncertainty().czx();
      ESV(1,2) = fitter.fitted_vtx_uncertainty().czy();
       
      EPV(0,0) = vertexHandle->front().covariance(0,0);
      EPV(1,1) = vertexHandle->front().covariance(1,1);
      EPV(2,2) = vertexHandle->front().covariance(2,2);
      EPV(0,1) = vertexHandle->front().covariance(0,1);
      EPV(0,2) = vertexHandle->front().covariance(0,2);
      EPV(1,2) = vertexHandle->front().covariance(1,2);



      double ct, ect;
      V0_Lifetime(pv,sv,EPV,ESV, 5.27932, pT, ct, ect);

      
      cand.addUserFloat("PDL", ct);
      cand.addUserFloat("ePDL", ect);



      // VARIABLES DE TABLA DE KAONES 
      cand.addUserFloat("k_DCASig", k_ptr->userFloat("DCASig") );
      cand.addUserFloat("k_nValidHits", k_ptr->userInt("nValidHits") );
      cand.addUserInt("k_isMatchedToMuon", k_ptr->userInt("isMatchedToMuon") );
      cand.addUserInt("k_isMatchedToLooseMuon", k_ptr->userInt("isMatchedToLooseMuon") );
      cand.addUserInt("k_isMatchedToSoftMuon", k_ptr->userInt("isMatchedToSoftMuon") );
      cand.addUserInt("k_isMatchedToMediumMuon", k_ptr->userInt("isMatchedToMediumMuon") );
      cand.addUserInt("k_isMatchedToEle", k_ptr->userInt("isMatchedToEle") ); 

      cand.addUserInt("k_HighPurity", k_ptr->userInt("trackHighPurity") ); 
      cand.addUserInt("k_numberOfHits", k_ptr->userInt("numberOfHits") ); 
      cand.addUserInt("k_numberOfPixelHits", k_ptr->userInt("numberOfPixelHits") ); 
      cand.addUserInt("k_lostInnerHits", k_ptr->userInt("lostInnerHits") ); 

      // std::cout<<"\nBToKsMuMuBuilder\n-------------\ntrackHighPurity  " << k_ptr->userInt("trackHighPurity") << "\n";
      // std::cout<<"numberOfHits  " << k_ptr->userInt("numberOfHits") << "\n";
      // std::cout<<"numberOfPixelHits  " << k_ptr->userInt("numberOfPixelHits") << "\n";
      // std::cout<<"lostInnerHits  " << k_ptr->userInt("lostInnerHits") << "\n";
      // std::cout << "------------\n\n";

 
      // VARIABLES DE TABLA DE MUONES
      // cand.addUserInt("l1_isPFMuon", l1_ptr->userInt("isPFMuon"));
      // cand.addUserInt("l1_isGlobalMuon", l1_ptr->userInt("isGlobalMuon"));
      // cand.addUserInt("l1_isTrackerMuon", l1_ptr->userInt("isTrackerMuon"));
      // cand.addUserInt("l1_isTriggering", l1_ptr->userInt("isTriggering"));
      // cand.addUserInt("l1_isSoft", l1_ptr->userInt("isSoft"));

      // cand.addUserInt("l2_isPFMuon", l2_ptr->userInt("isPFMuon"));
      // cand.addUserInt("l2_isGlobalMuon", l2_ptr->userInt("isGlobalMuon"));
      // cand.addUserInt("l2_isTrackerMuon", l2_ptr->userInt("isTrackerMuon"));
      // cand.addUserInt("l2_isTriggering", l2_ptr->userInt("isTriggering"));
      // cand.addUserInt("l2_isSoft", l2_ptr->userInt("isSoft"));



      ret_val->push_back(cand);




    } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
  
 //   }
   // catch(...){
     // std::cout << "Error!!!\n  "; 
    //}  

  } // for(size_t k_idx = 0; k_idx < kaons->size(); ++k_idx)



  for (auto & cand: *ret_val){
    cand.addUserInt("n_k_used", std::count(used_trk_id.begin(),used_trk_id.end(),cand.userInt("k_idx")));
    cand.addUserInt("n_l1_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l1_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l1_idx")));
    cand.addUserInt("n_l2_used", std::count(used_lep1_id.begin(),used_lep1_id.end(),cand.userInt("l2_idx"))+std::count(used_lep2_id.begin(),used_lep2_id.end(),cand.userInt("l2_idx")));
  }


  evt.put(std::move(pvTable), "VertexTable");
  evt.put(std::move(ret_val), "Bcollection");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToKsMuMuBuilder);