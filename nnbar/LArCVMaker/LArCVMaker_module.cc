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
#include "lardataobj/Simulation/SupernovaTruth.h"

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

namespace nnbar {

class LArCVMaker : public art::EDAnalyzer {

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

  larcv::IOManager fMgr;

  std::string fWireModuleLabel;
  std::string fMCTruthModuleLabel;
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

  std::map<int,std::vector<float>> fWireMap;
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
    fDataFileName(pset.get<std::string>("DataFileName")),
    fSpectrumFileName(pset.get<std::string>("SpectrumFileName")),
    fMaxTick(pset.get<int>("MaxTick")),
    fADCCut(pset.get<int>("ADCCut")),
    fEventType(pset.get<int>("EventType"))
{} // function LArCVMaker::LArCVMaker

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
    //std::cout <<"mother: "<< mct->GetParticle(i).Mother() << std::endl;
    //std::cout << "pdg: " <<mct->GetParticle(i).PdgCode() <<std::endl;
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

  for ( auto i = 0; i < mct->NParticles(); i++ ) {
    //std::cout <<"mother: "<< mct->GetParticle(i).Mother() << std::endl;
    //std::cout << "pdg: " <<mct->GetParticle(i).PdgCode() <<std::endl;
  }
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

  fMgr.set_id(evt.id().run(),evt.id().subRun(),evt.id().event());

  //int fVertexAPA =-1;

  // get wire objects
  art::Handle<std::vector<raw::RawDigit>> wireh;
  evt.getByLabel(fWireModuleLabel,wireh);

  // initialize ROI finding variables
  std::vector<int> apas;

  // fill wire map
  for (std::vector<raw::RawDigit>::const_iterator it = wireh->begin();
      it != wireh->end(); ++it) {
    const raw::RawDigit & rawr = *it;


    raw::RawDigit::ADCvector_t adc(rawr.Samples());
    raw::Uncompress(rawr.ADCs(), adc, rawr.Compression());

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
  auto images = (larcv::EventImage2D*)(fMgr.get_data(larcv::kProductImage2D, "tpc"));
  std::cout << std::endl;
  for (int it_plane = 2; it_plane < 3; ++it_plane) {
    int downsample = FindROI(best_apa,it_plane);
    std::cout << "downsampleing? " << downsample << std::endl;
    std::cout << "PLANE " << it_plane << " IMAGE" << std::endl;
    std::cout << "Original image resolution " << fNumberWires << "x" << fNumberTicks;
    larcv::Image2D image(fNumberWires/2,fNumberTicks);//fNumberWires ->
    for (int it_channel = 0; it_channel < fNumberWires/2; ++it_channel) {
      int channel = it_channel + fFirstWire;
      if (best_tpc%2 == 1){
	channel = it_channel + fFirstWire + 480;
      }
      for (int it_tick = 0; it_tick < fNumberTicks; ++it_tick) {
        int tick = it_tick + fFirstTick;
        if (fWireMap.find(channel) != fWireMap.end()) {
	  image.set_pixel(it_channel,it_tick,fWireMap[channel][tick]);

	  //  if (it_plane ==2)
	  if (fWireMap[channel][tick]!=0) {
	    hADCSpectrum->Fill(fWireMap[channel][tick]);
	  }
	  //std::cout << "it_channel : "<< it_channel << " , it_tick : " << it_tick << " , value : " << fWireMap[channel][tick] <<std::endl;
	}
      }
    }
    //yj commented this out june 20th 2019
    // image.compress(fNumberWires/downsample,fNumberTicks/(4*downsample));
    //std::cout << " => downsampling to " << fNumberWires/downsample << "x" << fNumberTicks/(4*downsample) << "." << std::endl << std::endl;
    //image.resize(600,600,0);
    //std::cout << "resized to 600 x 600" << std::endl;
    images->Emplace(std::move(image));
    std::cout << "emplace done" << std::endl;
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

DEFINE_ART_MODULE(LArCVMaker)

} // namespace nnbar
