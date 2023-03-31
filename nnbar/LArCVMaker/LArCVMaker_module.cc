// framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// data product includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

#include "nusimdata/SimulationBase/MCTruth.h"
//#include "lardataobj/Simulation/SupernovaTruth.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"


// root includes
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

// c++ includes
#include <vector>
#include <iterator>
#include <typeinfo>
#include <memory>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>

// local includes
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/EventROI.h"
#include "larcv/core/DataFormat/IOManager.h"

#include <iostream>
#include <fstream>

//#define DEBUG

#ifndef FOR
#define FOR(i, size) for (unsigned int i = 0; i < size; ++i)
#endif

namespace nnbar {

class LArCVMaker : public art::EDAnalyzer {

public:
	enum {
		kKaon = 321,
		kK0 = 311,
		kK0S = 310,
		kK0L = 130,
		kProton = 2212,
		kNeutron = 2112,
		kPi = 211,
		kPi0 = 111,
		kMu = 13,
		kEl = 11,
		kGamma = 22
	};

	std::set<unsigned int> fVisibleParticles;
	std::map<unsigned int, const char* > fParticleLabel;

public:

  explicit LArCVMaker(fhicl::ParameterSet const & pset);
  void analyze(art::Event const& evt);
  void beginJob();
  void endJob();


private:

  void ClearData();
  void ResetROI();
  void SetROISize();
  int FindAPAWithNeutrino(std::vector<int> apas, art::Event const & evt);
  int FindTPCWithNeutrino(std::vector<int> apas, art::Event const & evt);
  int FindROI(int apa, int plane);
	int Preselect(art::Event const& evt);
  larcv::IOManager fMgr;

  std::string fWireModuleLabel;
  std::string fMCTruthModuleLabel;
  std::string fG4ModuleLabel;

  std::string fHitModuleLabel;
  std::string fClusterModuleLabel;
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;

  std::string fDataFileName;
  std::string fSpectrumFileName;

  int fMaxTick;
  int fADCCut;
  int fEventType;

  int fFirstWire;
  int fLastWire;
  int fFirstTick;
  int fLastTick;

  const int fNumberChannels[3] = { 800, 800, 960 };
  const int fFirstChannel[3] = { 0, 800, 1600 };
  const int fLastChannel[3] = { 799, 1599, 2559 };

  int fEvent;
  int fAPA;
  int fNumberWires;
  int fNumberTicks;

	int fDoPreselection;
	int fNHitsMax;
	int fNHitsMin;
	int fNClustersMax;
	int fNClustersMin;
	int fNTracksMax;
	int fNTracksMin;
	int fNShowersMax;
	float fTrackLenMax;

  std::map<int, std::vector<float> > fWireMap;
  //std::ofstream pdg;
  TH1D* hADCSpectrum;
  TFile* SpectrumFile;
  std::ofstream pdg;
}; // class LArCVMaker

LArCVMaker::LArCVMaker(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMgr(larcv::IOManager::kWRITE),
    fWireModuleLabel(pset.get<std::string>("WireModuleLabel")),
    fMCTruthModuleLabel(pset.get<std::string>("MCTruthModuleLabel")),
    fG4ModuleLabel(pset.get<std::string>("G4ModuleLabel")),
    fHitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel")),
    fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    fShowerModuleLabel(pset.get<std::string>("ShowerModuleLabel")),
    fDataFileName(pset.get<std::string>("DataFileName")),
    fSpectrumFileName(pset.get<std::string>("SpectrumFileName")),
    fMaxTick(pset.get<int>("MaxTick")),
    fADCCut(pset.get<int>("ADCCut")),
    fEventType(pset.get<int>("EventType")),
    fDoPreselection(pset.get<int>("DoPreselection")),
    fNHitsMax(pset.get<int>("NHitsMax")),
    fNHitsMin(pset.get<int>("NHitsMin")),
    fNClustersMax(pset.get<int>("NClustersMax")),
    fNClustersMin(pset.get<int>("NClustersMin")),
    fNTracksMax(pset.get<int>("NTracksMax")),
    fNTracksMin(pset.get<int>("NTracksMin")),
    fNShowersMax(pset.get<int>("NShowersMax")),
    fTrackLenMax(pset.get<float>("TrackLenMax"))
{

	unsigned int tmp_list[] = {
		kKaon,
		kK0,
		kK0S,
		kK0L,
		kProton,
		kNeutron,
		kPi,
		kPi0,
		kMu,
		kEl,
		kGamma
	};

	fVisibleParticles = std::set<unsigned int>(tmp_list, tmp_list+11);

	fParticleLabel = {
		{kKaon, "K"},
		{kK0, "K0"},
		{kK0S, "K0S"},
		{kK0L, "K0L"},
		{kProton, "p"},
		{kNeutron, "n"},
		{kPi, "pi"},
		{kPi0, "pi0"},
		{kMu, "mu"},
		{kEl, "el"},
		{kGamma, "gamma"}
	};

	for (auto pair: fParticleLabel) {
		std::cout<<pair.first<<" "<<pair.second<<std::endl;
	}


	std::cout<<"LArCVMaker: will do preselection = "<<fDoPreselection<<std::endl;
} // function LArCVMaker::LArCVMaker



void LArCVMaker::beginJob() {
  std::ofstream pdg;

  std::string filename;
  if (std::getenv("PROCESS") != nullptr) filename = "larcv_" + std::string(std::getenv("PROCESS")) + ".root";
  else filename = fDataFileName;
  fMgr.set_out_file(filename);
  fMgr.initialize();
  SpectrumFile = new TFile(fSpectrumFileName.c_str(),"RECREATE");
  hADCSpectrum = new TH1D("hADCSpectrum","ADC Spectrum Collection; ADC; Entries",4096, 0., 4096.);
} // function LArCVMaker::beginJob




void LArCVMaker::endJob() {

  SpectrumFile->cd();
  hADCSpectrum->Write();
  SpectrumFile->Close();
  fMgr.finalize();
} // function LArCVMaker::endJob




void LArCVMaker::ClearData() {

  ResetROI();
  fAPA = -1;
  fWireMap.clear();
} // function LArCVMaker::ClearData

void LArCVMaker::ResetROI() {

  fFirstWire = -1;
  fLastWire = -1;
  fFirstTick = -1;
  fLastTick = -1;
  fNumberWires = -1;
  fNumberTicks = -1;
} // function LArCVMaker::ResetROI

void LArCVMaker::SetROISize() {

  fNumberWires = fLastWire - fFirstWire + 1;
  fNumberTicks = fLastTick - fFirstTick + 1;
}




int LArCVMaker::FindAPAWithNeutrino(std::vector<int> apas, art::Event const & evt) {

  int fVertexAPA =-1;
  int best_apa = -1;

  art::Handle<std::vector<simb::MCTruth>> TruthListHandle;
  std::vector<art::Ptr<simb::MCTruth>> TruthList;
  if (evt.getByLabel("marley",TruthListHandle))//

  art::fill_ptr_vector(TruthList,TruthListHandle);
  art::Ptr<simb::MCTruth> mct = TruthList.at(0);

  for ( auto i = 0; i < mct->NParticles(); i++ ) {
      std::cout <<"MC Truth particle "<<i<<":"<< std::endl;
      std::cout <<"  mother: "<< mct->GetParticle(i).Mother() << std::endl;
      std::cout <<"  pdg: " <<mct->GetParticle(i).PdgCode() << std::endl;
  }
  TVector3  vertex_position = mct->GetParticle(0).Position(0).Vect();
  //std::cout<<"mother code picked: "<<mct->GetParticle(0).Mother();
  //std::cout<<"position: "<<vertex_position.x() <<"," <<vertex_position.y() <<"," <<vertex_position.z() <<std::endl;


  art::ServiceHandle<geo::Geometry> geo;
  unsigned int itpc = 0;
  for (auto tpc_it = geo->begin<geo::TPCGeo>(); tpc_it != geo->end<geo::TPCGeo>(); ++tpc_it) {
		const geo::TPCGeo & tpc = *tpc_it;

		if (tpc.ContainsPosition(vertex_position)) {
			fVertexAPA = std::floor((float)itpc/2); // This assumes 1x2x6 workspace geometry, where 1 APA is shared by 2 TPCs.
			break;
		}
		++itpc;
  }


  best_apa = fVertexAPA;
  //std::cout<< best_apa <<std::endl;

  return best_apa;

} // function LArCVMaker::FindBestAPA



int LArCVMaker::FindTPCWithNeutrino(std::vector<int> tpcs, art::Event const & evt) {

  int fVertexTPC =-1;
  int best_tpc = -1;

  art::Handle<std::vector<simb::MCTruth>> TruthListHandle;
  std::vector<art::Ptr<simb::MCTruth>> TruthList;
  if (evt.getByLabel(fMCTruthModuleLabel, TruthListHandle))

  art::fill_ptr_vector(TruthList,TruthListHandle);
  art::Ptr<simb::MCTruth> mct = TruthList.at(0);

  // int npi = 0;
  // int npi0 = 0;
  // int nprot = 0;
  // int nneutron = 0;
  // double EKaon = 0.;
  // double EhadTot = 0.;

	std::vector<unsigned int> primaries;
	std::vector<float> primary_energies;
	//	unsigned int ikaon = 999;

  for ( auto i = 0; i < mct->NParticles(); i++ ) {
      if ( mct->GetParticle(i).StatusCode() != 1 ) continue;
      auto mom = mct->GetParticle(i).Momentum();

      int pdg = mct->GetParticle(i).PdgCode();
			if (fVisibleParticles.find(abs(pdg)) == fVisibleParticles.end()) continue;

			if (pdg == kEl && mom.E() - mom.M() < 0.6e-3) continue;

			primaries.push_back(abs(pdg));
			primary_energies.push_back(mom.E() - mom.M()); // in GeV

			// if (pdg == kKaon)
			// 	ikaon = i;

#ifdef DEBUG
      std::cout <<"MC Truth particle "<<i<<":"<< std::endl;
      std::cout <<"  pdg: " <<mct->GetParticle(i).PdgCode() << std::endl;
      std::cout <<"  status: "<< mct->GetParticle(i).StatusCode() << std::endl;
      std::cout <<"  KinE: " << (mom.E() - mom.M()) << std::endl;
      std::cout <<"  track id: "<< mct->GetParticle(i).TrackId() << std::endl;
      std::cout <<"  mother: "<< mct->GetParticle(i).Mother() << std::endl;
#endif
  }

  // print out truth summary
	std::cout<<"Primaries: "<<primaries.size()<<std::endl;
  std::cout<<"p -> ";
	std::cout.precision(3);
	for(size_t i = 0; i < primaries.size(); ++i) {
		//	std::cout<<i<<" "<<primaries[i]<<" "<<primary_energies[i]<<std::endl;
		if (i > 0)
			std::cout<<" + ";
		std::cout<<fParticleLabel[primaries[i]]<<"("<<(primary_energies[i]*1e3)<<" MeV)";
	}
	std::cout<<std::endl;


  TVector3  vertex_position = mct->GetParticle(0).Position(0).Vect();
  //std::cout<<"mother code picked: "<<mct->GetParticle(0).Mother();
  //std::cout<<"position: "<<vertex_position.x() <<"," <<vertex_position.y() <<"," <<vertex_position.z() <<std::endl;


  art::ServiceHandle<geo::Geometry> geo;
  unsigned int itpc = 0;
  for (auto tpc_it = geo->begin<geo::TPCGeo>(); tpc_it != geo->end<geo::TPCGeo>(); ++tpc_it) {
		const geo::TPCGeo & tpc = *tpc_it;

		if (tpc.ContainsPosition(vertex_position)) {
			fVertexTPC = itpc;
			break;
		}
		++itpc;
  }

  best_tpc = fVertexTPC;
  //std::cout<< best_tpc <<std::endl;


  // Extract info about output K destiny
  art::Handle<std::vector<simb::MCParticle>> G4ListHandle;
  std::vector<art::Ptr<simb::MCParticle>> PartList;
  if (evt.getByLabel(fG4ModuleLabel, G4ListHandle))
      art::fill_ptr_vector(PartList, G4ListHandle);

	std::vector<unsigned int> kdaughters;
	std::vector<float> daughter_energies;

  // npi = 0;
  // npi0 = 0;
  for (art::Ptr<simb::MCParticle> part: PartList) {
		if ( abs(part->PdgCode()) != kKaon ) continue; // interested in the kaon
		if ( strcmp(part->EndProcess().c_str(), "Decay") ) continue; // this is not the decaying kaon... maybe inelastic scattering?

		auto trkid = part->TrackId();

#ifdef DEBUG
		auto mom = part->Momentum();
		std::cout << "Kaon end process: " << part->EndProcess() << std::endl;
		std::cout <<"    pdg: " <<part->PdgCode() << std::endl;
		std::cout <<"    status: "<< part->StatusCode() << std::endl;
		std::cout <<"    KinE: " << (mom.E() - mom.M()) << std::endl;
		std::cout <<"    track id: "<< part->TrackId() << std::endl;
		std::cout <<"    mother: "<< part->Mother() << std::endl;

		// print out daughters
		std::cout <<"Kaon daughters (" << part->NumberDaughters() << "): " << std::endl;
		std::cout <<"(particle list size: " << PartList.size() << ")" << std::endl;
		for (int idaugh = 0; idaugh < part->NumberDaughters(); ++idaugh) {
			auto d = PartList.at(part->Daughter(idaugh)-1);
			auto mom = d->Momentum();
			std::cout <<"  Daughter " << idaugh << "(" << part->Daughter(idaugh) << "):" << std::endl;
			std::cout <<"    pdg: " <<d->PdgCode() << std::endl;
			std::cout <<"    status: "<< d->StatusCode() << std::endl;
			std::cout <<"    KinE: " << (mom.E() - mom.M()) << std::endl;
			std::cout <<"    track id: "<< d->TrackId() << std::endl;
			std::cout <<"    mother: "<< d->Mother() << std::endl;
		}
#endif

		for (int idaugh = 0; idaugh < part->NumberDaughters(); ++idaugh) {
			auto d = PartList.at(part->Daughter(idaugh)-1);
			if (d->Mother() != trkid) continue; // this is not kaon's direct descendant

			auto pdg = abs(d->PdgCode());
			auto mom = d->Momentum();
			auto kine = mom.E() - mom.M();
			if (pdg == kEl && kine < 0.6e-3) continue;
			if (fVisibleParticles.find(abs(pdg)) == fVisibleParticles.end()) continue;

			kdaughters.push_back(pdg);
			daughter_energies.push_back(kine);
		}

		break;
  }

	if (kdaughters.size() > 0) {
		// print out kaon decay
		std::cout<<"K -> ";
		for(long unsigned i = 0; i < kdaughters.size(); ++i) {
			if (i > 0)
				std::cout<<" + ";
			std::cout<<fParticleLabel[kdaughters[i]]<<"("<<(daughter_energies[i]*1e3)<<" MeV)";
		}
		std::cout<<std::endl;
	}

  return best_tpc;

} // function LArCVMaker::FindBestAPA

int LArCVMaker::FindROI(int best_apa, int plane) {

  ResetROI();

  //std::cout << "setting ROI as full APA" <<std::endl;

  fFirstTick=0;
  fLastTick=4487;

  if (plane==0) {
    fFirstWire = (2560*best_apa);
    fLastWire = (2560*best_apa)+799;
  }

  if (plane==1) {
    fFirstWire = (2560*best_apa)+800;
    fLastWire = (2560*best_apa)+1599;
  }

  if (plane==2) {
    fFirstWire = (2560*best_apa)+1600;
    fLastWire = (2560*best_apa)+2559;
  }
  //std::cout <<"ROI set!"<<std::endl;


  SetROISize();

  int downsample = 1;
  return downsample;
} // function LArCVMaker::FindROI

void LArCVMaker::analyze(art::Event const & evt) {

  ClearData();


  fEvent = evt.event();

  // Do pre-selection??
	if (fDoPreselection && !Preselect(evt)) {
		//#ifdef DEBUG
		std::cout<<"LArCVMaker: Event does not meet selection criteria. Skipping."<<std::endl;
		//#endif
		return;
	}


  fMgr.set_id(evt.id().run(),evt.id().subRun(),evt.id().event());

  //int fVertexAPA =-1;

  // get wire objects
  art::Handle<std::vector<raw::RawDigit>> wireh;
  evt.getByLabel(fWireModuleLabel,wireh);

  // initialize ROI finding variables
  std::vector<int> apas;

  // initialize map of pedestals
  std::map<int,float> pedestalMap;

  // fill wire map
  for (std::vector<raw::RawDigit>::const_iterator it = wireh->begin();
      it != wireh->end(); ++it) {
    const raw::RawDigit & rawr = *it;


    raw::RawDigit::ADCvector_t adc(rawr.Samples());
    raw::Uncompress(rawr.ADCs(), adc, rawr.Compression());

    // get pedestal for the channel
    float pedestal = rawr.GetPedestal();
    pedestalMap.insert(std::pair<int,float>(rawr.Channel(),pedestal));

    // record adc for each channel
    fWireMap.insert(std::pair<int,std::vector<float>>(rawr.Channel(),std::vector<float>(adc.begin(),adc.end())));
    int apa = std::floor(rawr.Channel()/2560);
    if (std::find(apas.begin(),apas.end(),apa) == apas.end())
      apas.push_back(apa);
  }

  // find best APA
  if (apas.size() == 0) {
    std::cout << "Skipping event. No activity inside the TPC!" << std::endl;
    return;
  }
  //int best_apa = FindAPAWithNeutrino(apas,evt);

  int best_tpc = FindTPCWithNeutrino(apas,evt);
  //int best_tpc = rand() % 24;//for rad.

  std::cout << "best_tpc %2 " << best_tpc%2 << std::endl;

  int best_apa = std::floor((float)best_tpc/2);

  if (best_apa != -1) fAPA = best_apa;
  else {
    std::cout << "Skipping event. Could not find good APA!" << std::endl;
    return;
  }
  std::cout << fAPA<<std::endl;

  // check for problems
  for (int it_plane = 0; it_plane < 3; ++it_plane) {
    if (FindROI(best_apa,it_plane) == -1) {
      std::cout << "Skipping event. Could not find good ROI in APA!" << std::endl;
      return;
    }
  }

  // produce image
  const int nticks_skip = 80; // number of ticks to skip to cover for gap between collection planes on each side of the APA

  auto images = (larcv::EventImage2D*)(fMgr.get_data(larcv::kProductImage2D, "tpc"));
  std::cout << std::endl;
  for (int it_plane = 0; it_plane < 3; ++it_plane) { // run over all 3 planes
    int downsample = FindROI(best_apa,it_plane);
    std::cout << "downsampling? " << downsample << std::endl;
    std::cout << "PLANE " << it_plane << " IMAGE" << std::endl;
    std::cout << "Original image resolution " << fNumberWires << "x" << fNumberTicks << std::endl;

		// std::cout<<"DBG: fNumberWires = "<<fNumberWires<<", it_plane = "<<it_plane<<std::endl;
		int rows = (it_plane == 2) ? fNumberWires/2               : fNumberWires;
		int cols = (it_plane == 2) ? 2*fNumberTicks + nticks_skip : fNumberTicks;
		// std::cout<<"DBG: rows = "<<rows<<", cols = "<<cols<<std::endl;
    larcv::Image2D image(rows, cols);// read TPCs on both sides of APA

		if (it_plane == 2) {
			for (int it_tpc = 0; it_tpc < 2; ++it_tpc) {
				int which_tpc = (best_tpc/2)*2 + it_tpc;
				for (int it_channel = 0; it_channel < fNumberWires/2; ++it_channel) {
					int channel = it_channel + fFirstWire;
					if (which_tpc%2 == 1){
						channel = it_channel + fFirstWire + 480; // other side of the APA
					}
					for (int it_tick = 0; it_tick < fNumberTicks; ++it_tick) {
						int tick = it_tick + fFirstTick;
						if (fWireMap.find(channel) != fWireMap.end()) {
							float adc = fWireMap[channel][tick];
							image.set_pixel(it_channel,
															fNumberTicks + nticks_skip + it_tick*(1 - 2*(which_tpc%2)) - (nticks_skip+1)*(which_tpc%2),
															(adc > 0.)?adc - pedestalMap[channel]:0.);
							//  if (it_plane ==2)
							if (adc!=0.) {
								hADCSpectrum->Fill(adc);
							}
							// if (fWireMap[channel][tick])
							//     std::cout << "it_channel : "<< it_channel << " , it_tick : " << it_tick << " , value : " << fWireMap[channel][tick] <<std::endl;
						}
					}
				}
			}
		} else {
			for (int it_channel = 0; it_channel < fNumberWires; ++it_channel) {
				int channel = it_channel + fFirstWire;
				for (int it_tick = 0; it_tick < fNumberTicks; ++it_tick) {
					int tick = it_tick + fFirstTick;
					if (fWireMap.find(channel) != fWireMap.end()) {
						float adc = fWireMap[channel][tick];
						image.set_pixel(it_channel,
														it_tick,
														(adc > 0)?adc - pedestalMap[channel]:0);
						//  if (it_plane ==2)
						if (adc!=0) {
							hADCSpectrum->Fill(adc);
						}
						// if (fWireMap[channel][tick])
						//     std::cout << "it_channel : "<< it_channel << " , it_tick : " << it_tick << " , value : " << fWireMap[channel][tick] <<std::endl;
					}
				}
			}
		}

    //yj commented this out june 20th 2019
    // image.compress(fNumberWires/downsample,fNumberTicks/(4*downsample));
    //std::cout << " => downsampling to " << fNumberWires/downsample << "x" << fNumberTicks/(4*downsample) << "." << std::endl << std::endl;
    //image.resize(600,600,0);
    //std::cout << "resized to 600 x 600" << std::endl;
    images->Emplace(std::move(image));
    std::cout << "emplace for plane "<<it_plane<<" done" << std::endl;
  }

  auto roi_v = (larcv::EventROI*)(fMgr.get_data(larcv::kProductROI, "tpc"));

  larcv::ROI roi((larcv::ROIType_t)fEventType);


  //  art::Handle<std::vector<simb::MCTruth>> TruthListHandle;
  //  std::vector<art::Ptr<simb::MCTruth>> TruthList;

  //  if (evt.getByLabel("nnbar",TruthListHandle))
  //    art::fill_ptr_vector(TruthList,TruthListHandle);

  //  art::Ptr<simb::MCTruth> mct = TruthList.at(0);
  //  std::cout<< "mct nparticles " << mct->NParticles() << std::endl;


  std::cout << "autoroi done" << std::endl;
  roi_v->Emplace(std::move(roi));
  std::cout << "second emplace done" << std::endl;
  fMgr.save_entry();
  std::cout << "save entry done" << std::endl;
} // function LArCVMaker::analyze

int LArCVMaker::Preselect(art::Event const& evt)
{
	int nhits = 0;
	int nclusters = 0;
	int ntracks = 0;
	int nshowers = 0;
	float maxtracklen = 0.;

#ifdef DEBUG

		auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	if (hitListHandle)
		nhits = hitListHandle->size();
	auto clusterListHandle = evt.getHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);
	if (clusterListHandle)
		nclusters = clusterListHandle->size();


	auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
	if (trackListHandle)
		ntracks = trackListHandle->size();

	auto showerListHandle = evt.getHandle< std::vector<recob::Shower> >(fShowerModuleLabel);
	if (showerListHandle)
		nshowers = showerListHandle->size();

	// iterate over reco tracks and get max track length
	std::vector<art::Ptr<recob::Track> > tracklist;
	if (trackListHandle)
		art::fill_ptr_vector(tracklist, trackListHandle);

	for (auto trk: tracklist) {
		if (trk->Length() > maxtracklen)
			maxtracklen = trk->Length();
	}

	std::cout<<"LArCVMaker: Basic variables:"<<std::endl;
	std::cout<<"            nhits = "<<nhits<<std::endl
					 <<"            nclusters = "<<nclusters<<std::endl
					 <<"            ntracks = "<<ntracks<<std::endl
					 <<"            nshowers = "<<nshowers<<std::endl
					 <<"            maxtracklen = "<<maxtracklen<<std::endl;

	if (nhits > fNHitsMax || nhits < fNHitsMin) return 0;
	if (nclusters > fNClustersMax || nclusters < fNClustersMin) return 0;
	if (ntracks > fNTracksMax || ntracks < fNTracksMin) return 0;
	if (nshowers > fNShowersMax) return 0;
	if (maxtracklen > fTrackLenMax) return 0;

#else

	// get product handles for hits, tracks, clusters and showers
	//
	// std::vector<art::Ptr<recob::Hit> > hitlist;
	auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	if (hitListHandle)
		nhits = hitListHandle->size();

	if (nhits > fNHitsMax || nhits < fNHitsMin) return 0;

	auto clusterListHandle = evt.getHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);
	if (clusterListHandle)
		nclusters = clusterListHandle->size();

	if (nclusters > fNClustersMax || nclusters < fNClustersMin) return 0;

	auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
	if (trackListHandle)
		ntracks = trackListHandle->size();

	if (ntracks > fNTracksMax || ntracks < fNTracksMin) return 0;


	auto showerListHandle = evt.getHandle< std::vector<recob::Shower> >(fShowerModuleLabel);
	if (showerListHandle)
		nshowers = showerListHandle->size();

	if (nshowers > fNShowersMax) return 0;

	// iterate over reco tracks and get max track length
	std::vector<art::Ptr<recob::Track> > tracklist;
	if (trackListHandle)
		art::fill_ptr_vector(tracklist, trackListHandle);

	for (auto trk: tracklist) {
		if (trk->Length() > maxtracklen)
			maxtracklen = trk->Length();
	}

	if (maxtracklen > fTrackLenMax) return 0;

#endif
	return 1;
}


DEFINE_ART_MODULE(LArCVMaker)

} // namespace nnbar

// Local Variables:
// tab-width: 2
// c-basic-offset: 2
// End:
