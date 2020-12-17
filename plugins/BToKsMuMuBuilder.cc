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
  const StringCutObjectSelector<pat::Muon> l1_selection_; // cut on leading lepton
  const StringCutObjectSelector<pat::Muon> l2_selection_; // cut on sub-leading lepton
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
  reco::Vertex bestVtx;
  
  for (const reco::Vertex & vertex : *vertexHandle){
    for (const pat::Muon & trgmu : *trgMuons){
      if(nv==3) break;
      if(fabs(vertex.position().z()-trgmu.vz())<dzCut_) {
        dzTrgMu.push_back(fabs(vertex.position().z()-trgmu.vz()));
        bestVtx = vertex;
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
	    if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; 
        
        //taken from DiLeptonBuilder
        if(!l1_selection_((iMuon1->pt() > iMuon2->pt() ? *iMuon1 : *iMuon2))) continue;
        if(!l2_selection_((iMuon1->pt() < iMuon2->pt() ? *iMuon1 : *iMuon2))) continue;
        
        pat::CompositeCandidate lepton_pair;
        lepton_pair.setP4(iMuon1->p4() + iMuon2->p4());
        lepton_pair.setCharge(iMuon1->charge() + iMuon2->charge());
        lepton_pair.addUserFloat("lep_deltaR", reco::deltaR(*iMuon1, *iMuon2));
        
        // Use UserCands as they should not use memory but keep the Ptr itself
        float dmz = (iMuon1->vz() - iMuon2->vz());
        lepton_pair.addUserFloat("lept_DZ", dmz);
        if( !DLB_pre_vtx_selection_(lepton_pair) ) continue;



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
        pat::CompositeCandidate b_cand;
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
	    //ParticleMass psi_mass = 3.096916;
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
	        //continue; 
	    }
      try{
	      psiVertexFitTree->movePointerToTheTop();
      }
      catch (...){
        std::cout<< "PsiVertexEmpty" << std::endl;
        std::cout<< "muonParticles size: " << muonParticles.size() << std::endl;
        continue;
      }
	    RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	    RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	    
	    if( psi_vFit_vertex_noMC->chiSquared() < 0 ){
	        //std::cout << "negative chisq from psi fit" << endl;
	        //continue;
	    }

        double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	    //if(J_Prob_tmp<0.01) //Jhovanny
	    //if(J_Prob_tmp<1.e-5) //DiMuonBuilder postVtxSelection
	    //{
	    //   continue;
	    //}
	  
	   //some loose cuts go here
       lepton_pair.addUserFloat("sv_chi2", psi_vFit_vertex_noMC->chiSquared());
       lepton_pair.addUserFloat("sv_ndof", psi_vFit_vertex_noMC->degreesOfFreedom()); // float??
       lepton_pair.addUserFloat("sv_prob", J_Prob_tmp);
       lepton_pair.addUserFloat("fitted_mass", psiVertexFitTree->isValid() ? psi_vFit_noMC->currentState().mass() : -1);
       lepton_pair.addUserFloat("fitted_massErr", psiVertexFitTree->isValid()  ? sqrt(psi_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6)) : -1);
       if( !DLB_post_vtx_selection_(lepton_pair) ) continue;
	   //if(psi_vFit_vertex_noMC->chiSquared()>999) continue; // DiMuonBuilder cuts at 998 
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
                b_cand.setP4(iMuon1->p4() + iMuon2->p4() + iVee->daughter(0)->p4() + iVee->daughter(1)->p4());
                b_cand.setCharge(iMuon1->charge() + iMuon2->charge() + iVee->daughter(0)->charge() + iVee->daughter(0)->charge());
                if( !pre_vtx_selection_(b_cand) ) continue;

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
            try{
		          Ks0VertexFitTree->movePointerToTheTop();
            }
            catch (...){
              std::cout<< "K0 vertex Failed!"<< std::endl;
              std::cout<< "pionParticles size: "<<pionParticles.size() << std::endl; 
              continue;
            }   
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
            try{
              Ks0VertexFitTree->movePointerToTheTop();
            }
            catch(...){
              std::cout<< "Ks0Vertex Failed! (for MC)"<< std::endl;
              continue;
            }
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
            try{
              vertexFitTree->movePointerToTheTop();		     
            }
            catch(...){
              std::cout<< "VertexFitTree Failed!" << std::endl;
              std::cout<< "vFitMCParticles size: "<<vFitMCParticles.size() << std::endl;
              continue;
            }
		        RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		        RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		        //if (!bDecayVertexMC->vertexIsValid()){
		        //    //std::cout << "B MC fit vertex is not valid" << endl;
		        //    continue;
		        //}
    
		        //if(bCandMC->currentState().mass()<4.8 || bCandMC->currentState().mass()>6.0) continue;
    
		        //if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ){
			    //    //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
			    //    continue;
		        //}
		        //std::cout << "pass 461 continues ... "<< std::endl;
		        double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		        //if(B_Prob_tmp<0.01) //Jhovanny
		        //if(B_Prob_tmp<0.001){  //Horacio hardcoded
			    //    continue;
		        //}
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

                b_cand.addUserInt("Bvtx_OK" , bDecayVertexMC->vertexIsValid()); //sv_OK
                b_cand.addUserFloat("Bvtx_chi2", bDecayVertexMC->chiSquared()); //sv_chi2
                b_cand.addUserFloat("Bvtx_ndof", bDecayVertexMC->degreesOfFreedom()); // sv_ndof
                b_cand.addUserFloat("Bvtx_prob", B_Prob_tmp); // sv_prob
 
                b_cand.addUserFloat("fitted_mll" , psi_vFit_noMC->currentState().mass());
                b_cand.addUserFloat("fitted_pt_ll" , psi_vFit_noMC->currentState().globalMomentum().perp());
                b_cand.addUserFloat("fitted_eta_ll" , psi_vFit_noMC->currentState().globalMomentum().eta());
                b_cand.addUserFloat("fitted_phi_ll" , psi_vFit_noMC->currentState().globalMomentum().phi());
      
                b_cand.addUserFloat("Bfitted_pt"  , bCandMC->currentState().globalMomentum().perp()); 
                // cand.addUserFloat("fitted_px"  , fitter.fitted_candidate().globalMomentum().x()); 
                // cand.addUserFloat("fitted_py"  , fitter.fitted_candidate().globalMomentum().y()); 
                // cand.addUserFloat("fitted_pz"  , fitter.fitted_candidate().globalMomentum().z()); 
                b_cand.addUserFloat("Bfitted_eta" , bCandMC->currentState().globalMomentum().eta());
                b_cand.addUserFloat("Bfitted_phi" , bCandMC->currentState().globalMomentum().phi());
                b_cand.addUserFloat("Bfitted_mass", bCandMC->currentState().mass());      
                b_cand.addUserFloat("Bfitted_massErr", sqrt(bCandMC->currentState().kinematicParametersError().matrix()(6,6)));
                //std::cout << "cos2D: "<< cos_theta_2D(fitter, *beamspot, cand.p4()) << std::endl;
                //TLorentzVector testVect;
                auto B_vect = math::PtEtaPhiMLorentzVector(psi_vFit_noMC->currentState().globalMomentum().perp(),
                                                           psi_vFit_noMC->currentState().globalMomentum().eta(),
                                                           psi_vFit_noMC->currentState().globalMomentum().phi(),
                                                           psi_vFit_noMC->currentState().mass());
                //testVect.SetPtEtaPhiM(cand.pt(), cand.eta(), cand.phi(), cand.mass());
                GlobalPoint b_gp = bDecayVertexMC->position();
                b_cand.addUserFloat("Bcos_theta_2D", (bDecayVertexMC->vertexIsValid()) ?  my_cos_theta_2D(b_gp, *beamspot, b_cand.p4()) : -2 );
                b_cand.addUserFloat("Bfitted_cos_theta_2D", (bDecayVertexMC->vertexIsValid()) ?  my_cos_theta_2D(b_gp, *beamspot, B_vect) : -2 );
                //std::cout << "my cos2D: "<< my_cos_theta_2D(gp, *beamspot, testVect) << std::endl;
                //FIX
                if( !post_vtx_selection_(b_cand) ) continue;
                // fill candidate variables now                      

                GlobalError b_gp_err = bDecayVertexMC->error();
                auto lxy = my_l_xy(b_gp, b_gp_err, *beamspot);
                b_cand.addUserFloat("B_l_xy", lxy.value());
                b_cand.addUserFloat("B_l_xy_unc", lxy.error());

                b_cand.addUserFloat("Bvtx_x", b_cand.vx());
                b_cand.addUserFloat("Bvtx_y", b_cand.vy());
                b_cand.addUserFloat("Bvtx_z", b_cand.vz());
                b_cand.addUserFloat("Bvtx_ex" , sqrt(b_gp_err.cxx()));
                b_cand.addUserFloat("Bvtx_ey" , sqrt(b_gp_err.cyy()));
                b_cand.addUserFloat("Bvtx_ez" , sqrt(b_gp_err.czz()));
                try{
                  b_cand.addUserFloat("Bvtx_eyx", b_gp_err.cyx());
                  b_cand.addUserFloat("Bvtx_ezx", b_gp_err.czx());
                  b_cand.addUserFloat("Bvtx_ezy", b_gp_err.czy());
                }
                catch(...){
                  b_cand.addUserFloat("Bvtx_eyx", -1);
                  b_cand.addUserFloat("Bvtx_ezx", -1);
                  b_cand.addUserFloat("Bvtx_ezy", -1);
                }
                b_cand.addUserFloat("Bfitted_l1_pt" , mu1CandMC->currentState().globalMomentum().perp()); 
                b_cand.addUserFloat("Bfitted_l1_eta", mu1CandMC->currentState().globalMomentum().eta());
                b_cand.addUserFloat("Bfitted_l1_phi", mu1CandMC->currentState().globalMomentum().phi());
                b_cand.addUserFloat("Bl1_charge", mu1CandMC->currentState().particleCharge());

                b_cand.addUserFloat("Bfitted_l2_pt" , mu2CandMC->currentState().globalMomentum().perp()); 
                b_cand.addUserFloat("Bfitted_l2_eta", mu2CandMC->currentState().globalMomentum().eta());
                b_cand.addUserFloat("Bfitted_l2_phi", mu2CandMC->currentState().globalMomentum().phi());
                b_cand.addUserFloat("Bl2_charge", mu2CandMC->currentState().particleCharge());

                b_cand.addUserFloat("Bfitted_ks_pt"  , Ks0CandMC->currentState().globalMomentum().perp());
                b_cand.addUserFloat("Bfitted_ks_eta" , Ks0CandMC->currentState().globalMomentum().eta());
                b_cand.addUserFloat("Bfitted_ks_phi" , Ks0CandMC->currentState().globalMomentum().phi());
                b_cand.addUserFloat("Bfitted_ks_m" , Ks0_vFit_vertex_noMC->currentState().mass());
                b_cand.addUserFloat("Bks_charge", Ks0CandMC->currentState().particleCharge());            
                //compute isolation
                float l1_iso03 = iMuon1->trackIso();
                float l2_iso03 = iMuon2->trackIso();
                float l1_PFiso03 = getMuPFIso03(*iMuon1);
                float l1_PFiso04 = getMuPFIso04(*iMuon1);
                float l2_PFiso03 = getMuPFIso03(*iMuon2);
                float l2_PFiso04 = getMuPFIso04(*iMuon2);
                b_cand.addUserFloat("l1_iso03", l1_iso03);
                b_cand.addUserFloat("l1_PFiso03", l1_PFiso03);
                b_cand.addUserFloat("l1_PFiso04", l1_PFiso04);
                b_cand.addUserFloat("l2_iso03", l2_iso03);
                b_cand.addUserFloat("l2_PFiso03", l2_PFiso03);
                b_cand.addUserFloat("l2_PFiso04", l2_PFiso04);

                // Aqui creemos el boost al CM del dilepton
                auto dilep = math::XYZTLorentzVector(psi_vFit_noMC->currentState().globalMomentum().x(),
                                                     psi_vFit_noMC->currentState().globalMomentum().y(),
                                                     psi_vFit_noMC->currentState().globalMomentum().z(),
                                                     psi_vFit_noMC->currentState().mass());
                
                auto k0vec = math::XYZTLorentzVector(Ks0CandMC->currentState().globalMomentum().x(),
                                                     Ks0CandMC->currentState().globalMomentum().y(),
                                                     Ks0CandMC->currentState().globalMomentum().z(),
                                                     Ks0CandMC->currentState().mass());

                ROOT::Math::Boost cmboost(dilep.BoostToCM());
                
                math::XYZTLorentzVector kaonCM(  cmboost( k0vec )  );
                math::XYZTLorentzVector dimuCM(cmboost(getTLV(psi_vFit_noMC)));
                math::XYZTLorentzVector muonCMp, muonCMn; 
                //where the thetal is the angle between the
                //direction of the m-(m+) lepton and the K+(K-)
                //in this case we can calculate three anlges M-/Ks, M+/Ks, Dim/Ks
                if( mu1CandMC->currentState().particleCharge() > 0){ 
                    muonCMp = cmboost( getTLV(mu1CandMC) );
                    muonCMn = cmboost( getTLV(mu2CandMC) );
                }
                else {
                    muonCMp = cmboost( getTLV(mu2CandMC) );
                    muonCMn = cmboost( getTLV(mu1CandMC) );
                }
                float costhetaL = ( muonCMp.x()*muonCMn.x() 
                                  + muonCMp.y()*muonCMn.y() 
                                  + muonCMp.z()*muonCMn.z() ) / (muonCMp.P()*muonCMn.P() );

                float costhetaLpKs = ( muonCMp.x()*kaonCM.x() 
                                     + muonCMp.y()*kaonCM.y() 
                                     + muonCMp.z()*kaonCM.z() ) / (muonCMp.P()*kaonCM.P() );

                float costhetaLnKs = ( muonCMn.x()*kaonCM.x() 
                                     + muonCMn.y()*kaonCM.y() 
                                     + muonCMn.z()*kaonCM.z() ) / (muonCMn.P()*kaonCM.P() );

                float costhetaKsDM = ( dimuCM.x()*kaonCM.x()
                                     + dimuCM.y()*kaonCM.y()
                                     + dimuCM.z()*kaonCM.z() ) / (dimuCM.P()*kaonCM.P() );
                
                std::vector<float> cosAlpha(4, -2);
                std::vector<float> lxy_pv(4,-1);
                std::vector<float> errP(4,-1);
                
                b_cand.addUserFloat("costhetaLep", costhetaL);
                b_cand.addUserFloat("costhetaLpKs", costhetaLpKs);
                b_cand.addUserFloat("costhetaLnKs", costhetaLnKs);
                b_cand.addUserFloat("costhetaKsDM", costhetaKsDM);
                b_cand.addUserFloat("dRm1m2", dRm1m2);
                // CHEQUEMOS QUE EL CANDIDATO A "B" SATIFAGA LOS CORTES EN COS(ALPHA) Y SIGNIFICANCIA
                for( unsigned int iPV=0; iPV<dzTrgMu.size(); ++iPV ){ 
                    //Debemos calcular el vector que une los vertices primario y el ajustado
                    math::XYZVector vPS(b_gp.x()-vx.at(iPV), b_gp.y()-vy.at(iPV), b_gp.z()-vz.at(iPV));
                    //Momento espacial del candidato?????????
                    math::XYZVector Bp(b_gp.x(), b_gp.y(), b_gp.z());
                    //CosALPHA
                    cosAlpha[iPV] = vPS.Dot(Bp)/(vPS.R()*Bp.R());
        
                    //Para significancia:  
                    GlobalError err = b_gp_err;
                    GlobalPoint delta(b_gp.x()-vx.at(iPV), b_gp.y()-vy.at(iPV), 0.);  

                    lxy_pv[iPV] = delta.perp();
                    errP[iPV] = sqrt(err.rerr(delta));
                }

                b_cand.addUserFloat("cosAlpha0", cosAlpha[0]);
                b_cand.addUserFloat("cosAlpha1", cosAlpha[1]);
                b_cand.addUserFloat("cosAlpha2", cosAlpha[2]);
                b_cand.addUserFloat("cosAlpha3", cosAlpha[3]);

                b_cand.addUserFloat("lxy_pv0", lxy_pv[0]);
                b_cand.addUserFloat("lxy_pv1", lxy_pv[1]);
                b_cand.addUserFloat("lxy_pv2", lxy_pv[2]);
                b_cand.addUserFloat("lxy_pv3", lxy_pv[3]);

                b_cand.addUserFloat("significance0", lxy_pv[0]/errP[0]);
                b_cand.addUserFloat("significance1", lxy_pv[1]/errP[1]);
                b_cand.addUserFloat("significance2", lxy_pv[2]/errP[2]);
                b_cand.addUserFloat("significance3", lxy_pv[3]/errP[3]);
                
                //life time
                TVector3 pv, Bvtx, BpT;

                BpT.SetXYZ(bCandMC->currentState().globalMomentum().x(), bCandMC->currentState().globalMomentum().y(),0.0);
                pv.SetXYZ(vx.at(0),vy.at(0),vz.at(0));
                Bvtx.SetXYZ(b_gp.x(), b_gp.y(), b_gp.z());
                TMatrix ESV(3,3);
                TMatrix EPV(3,3);
                ESV(0,0) = bDecayVertexMC->error().cxx();
                ESV(1,1) = bDecayVertexMC->error().cyy();
                ESV(2,2) = bDecayVertexMC->error().czz();
                ESV(0,1) = bDecayVertexMC->error().cyx();
                ESV(0,2) = bDecayVertexMC->error().czx();
                ESV(1,2) = bDecayVertexMC->error().czy();
                EPV(0,0) = vertexHandle->front().covariance(0,0);
                EPV(1,1) = vertexHandle->front().covariance(1,1);
                EPV(2,2) = vertexHandle->front().covariance(2,2);
                EPV(0,1) = vertexHandle->front().covariance(0,1);
                EPV(0,2) = vertexHandle->front().covariance(0,2);
                EPV(1,2) = vertexHandle->front().covariance(1,2);
                //Vertice primario y error GAAS
                b_cand.addUserFloat("PV_x", vertexHandle->front().x());
                b_cand.addUserFloat("PV_y", vertexHandle->front().y());
                b_cand.addUserFloat("PV_z", vertexHandle->front().z());

                b_cand.addUserFloat("PV_ex",  vertexHandle->front().covariance(0,0));
                b_cand.addUserFloat("PV_ey",  vertexHandle->front().covariance(1,1));
                b_cand.addUserFloat("PV_ez",  vertexHandle->front().covariance(2,2));
                b_cand.addUserFloat("PV_eyx", vertexHandle->front().covariance(0,1));
                b_cand.addUserFloat("PV_ezx", vertexHandle->front().covariance(0,2));
                b_cand.addUserFloat("PV_ezy", vertexHandle->front().covariance(1,2));
                //info Ks0 y Piones
                b_cand.addUserFloat("Ks0_pt1",  Ks0p1vec.perp());
                b_cand.addUserFloat("Ks0_px1",  Ks0Pi1KP.momentum().x());
                b_cand.addUserFloat("Ks0_py1",  Ks0Pi1KP.momentum().y());
                b_cand.addUserFloat("Ks0_pz1",  Ks0Pi1KP.momentum().z());
                b_cand.addUserFloat("Ks0_px1_track", v0daughters[0].px());
                b_cand.addUserFloat("Ks0_py1_track", v0daughters[0].py());
                b_cand.addUserFloat("Ks0_pz1_track", v0daughters[0].pz());
                b_cand.addUserFloat("Ks0_p1_ch", T1CandMC->currentState().particleCharge());
                
                b_cand.addUserFloat("Ks0_pt2",  Ks0p2vec.perp());
                b_cand.addUserFloat("Ks0_px2",  Ks0Pi2KP.momentum().x());
                b_cand.addUserFloat("Ks0_py2",  Ks0Pi2KP.momentum().y());
                b_cand.addUserFloat("Ks0_pz2",  Ks0Pi2KP.momentum().z());
                b_cand.addUserFloat("Ks0_px2_track", v0daughters[1].px());
                b_cand.addUserFloat("Ks0_py2_track", v0daughters[1].py());
                b_cand.addUserFloat("Ks0_pz2_track", v0daughters[1].pz());
                b_cand.addUserFloat("Ks0_p2_ch", T2CandMC->currentState().particleCharge());


                double Bct, Bect;
                V0_Lifetime(pv, Bvtx, EPV, ESV, 5.27961, BpT, Bct, Bect);
                b_cand.addUserFloat("B_PDL", Bct);
                b_cand.addUserFloat("eB_PDL", Bect);
                
                TVector3 Kvtx, KpT;
                KpT.SetXYZ(Ks0CandMC->currentState().globalMomentum().x(), Ks0CandMC->currentState().globalMomentum().y(),0.0);
                Kvtx.SetXYZ(Ks0_vFit_vertex_noMC->position().x(), Ks0_vFit_vertex_noMC->position().y(), Ks0_vFit_vertex_noMC->position().z());
                TMatrix EKsV(3,3);
                EKsV(0,0) = Ks0_vFit_vertex_noMC->error().cxx();
                EKsV(1,1) = Ks0_vFit_vertex_noMC->error().cyy();
                EKsV(2,2) = Ks0_vFit_vertex_noMC->error().czz();
                EKsV(0,1) = Ks0_vFit_vertex_noMC->error().cyx();
                EKsV(0,2) = Ks0_vFit_vertex_noMC->error().czx();
                EKsV(1,2) = Ks0_vFit_vertex_noMC->error().czy();
                double Kct, Kect;
                V0_Lifetime(Bvtx, Kvtx, ESV, EKsV, 0.49761, KpT, Kct, Kect);
                b_cand.addUserFloat("K_PDL", Kct);
                b_cand.addUserFloat("eK_PDL", Kect);

                // Quality Vars 
                b_cand.addUserFloat("pi1_nValidPixelHits", theDaughterTracks[0].hitPattern().numberOfValidPixelHits());
                b_cand.addUserFloat("pi1_nPixelLWM",       theDaughterTracks[0].hitPattern().pixelLayersWithMeasurement());
                b_cand.addUserFloat("pi1_nTrackerLWM",     theDaughterTracks[0].hitPattern().trackerLayersWithMeasurement());
                b_cand.addUserFloat("pi2_nValidPixelHits", theDaughterTracks[1].hitPattern().numberOfValidPixelHits());
                b_cand.addUserFloat("pi2_nPixelLWM",       theDaughterTracks[1].hitPattern().pixelLayersWithMeasurement());
                b_cand.addUserFloat("pi2_nTrackerLWM",     theDaughterTracks[1].hitPattern().trackerLayersWithMeasurement());

                b_cand.addUserFloat("p1_HighPurity",     theDaughterTracks[0].quality(reco::TrackBase::highPurity));
                b_cand.addUserFloat("p2_HighPurity",     theDaughterTracks[1].quality(reco::TrackBase::highPurity));

                b_cand.addUserFloat("mu1_soft",  iMuon1->isSoftMuon(bestVtx)); 
                b_cand.addUserFloat("mu1_tight", iMuon1->isTightMuon(bestVtx));
                b_cand.addUserFloat("mu1_PF", iMuon1->isPFMuon());
                b_cand.addUserFloat("mu1_loose", muon::isLooseMuon(*iMuon1));
                b_cand.addUserFloat("mu1_medium", muon::isMediumMuon(*iMuon1));
                b_cand.addUserFloat("mu1_global", iMuon1->isGlobalMuon() );

                b_cand.addUserFloat("mu2_soft",  iMuon2->isSoftMuon(bestVtx)); 
                b_cand.addUserFloat("mu2_tight", iMuon2->isTightMuon(bestVtx));
                b_cand.addUserFloat("mu2_PF", iMuon2->isPFMuon());
                b_cand.addUserFloat("mu2_loose", muon::isLooseMuon(*iMuon2));
                b_cand.addUserFloat("mu2_medium", muon::isMediumMuon(*iMuon2));
                b_cand.addUserFloat("mu2_global", iMuon2->isGlobalMuon() );
                 

		        b_cand.addUserFloat("mum_C2", glbTrackM->normalizedChi2() );
		        b_cand.addUserFloat("mum_nValidHits", glbTrackM->numberOfValidHits() );
		        b_cand.addUserFloat("mum_nValidPixelHits", glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		        b_cand.addUserFloat("mup_C2", glbTrackP->normalizedChi2() );
		        b_cand.addUserFloat("mup_nValidHits", glbTrackP->numberOfValidHits() );
		        b_cand.addUserFloat("mup_nValidPixelHits", glbTrackP->hitPattern().numberOfValidPixelHits() );
                b_cand.addUserFloat("mum_dxy", glbTrackM->dxy(bestVtx.position()) );
		        b_cand.addUserFloat("mup_dxy", glbTrackP->dxy(bestVtx.position()) );
		        b_cand.addUserFloat("mum_dz" , glbTrackM->dz(bestVtx.position()) );
		        b_cand.addUserFloat("mup_dz" , glbTrackP->dz(bestVtx.position()) );
		        b_cand.addUserFloat("dimuon_dca" , dca);
            

                ret_val->push_back(b_cand);
            }// end V0 Tracks
        }// end if dimuon&& V0Tracks   
        passmu++;

    }    
  }// end muon loop   

  evt.put(std::move(pvTable), "VertexTable");
  evt.put(std::move(ret_val), "Bcollection");
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToKsMuMuBuilder);