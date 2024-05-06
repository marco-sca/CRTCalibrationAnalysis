/**
 * @file   CRTCalibrationAnalysis_module.cc
 * @brief  Access CRT data and reco products and compare to MCTruth info 
 * @author Chris Hilgenberg (Chris.Hilgenberg@colostate.edu)
 * 
 * The last revision of this code was done in October 2018 with LArSoft v07_06_01.
 */

// LArSoft includes
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"

// C++ includes
#include <map>
#include <vector>
#include <string>
#include <set>
#include <cmath>
#include <iostream>
#include <utility>
#include <array>
#include <algorithm>

// CRT data products
#include "sbnobj/ICARUS/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTCommonUtils.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;
using std::to_string;

//Is it good procedure to define in this way the struct? How can I define it in the scope of the function that use this struct type and keep it private?
/*
struct ChannelInfo {
   int feb;
   int channel;
	 int adc;
   bool broken;
   ChannelInfo(int feb, int channel, std::string name)
     	: feb(feb), channel(channel), name(name) {}
};
*/
namespace icarus {
namespace crt {

  class CRTCalibrationAnalysis : public art::EDAnalyzer
  {
  public:

    struct Config {
      
      // Save some typing:
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> CRTDAQLabel {
        Name("CRTDAQLabel"),
        Comment("tag of the input data product with calibrated CRT data")
        };
   

    }; // Config

    using Parameters = art::EDAnalyzer::Table<Config>;
    
  
    /// Constructor: configures the module (see the Config structure above)
    explicit CRTCalibrationAnalysis(Parameters const& config);

    virtual void beginJob() override;
    virtual void beginRun(const art::Run& run) override;
    virtual void analyze (const art::Event& event) override;
    void endJob() override;
    void endRun(const art::Run&) override;

  private:

    art::ServiceHandle<art::TFileService> tfs;

		//std::vector<ChannelInfo> brokenChannels;
		map<uint8_t, TH1F*> pedSigHistos;
		map<uint8_t, TH1F*> pedResHistos;
    map<uint8_t,vector<TH1F*>*> macToHistos;
		map<uint8_t,vector<TH1F*>*> macToHistos_nonTriggeringPed;
    map<uint8_t,vector<TH1F*>*> macToHistos_ped;
    map<uint8_t,vector<TH1F*>*> macToHistos_sig;
    TRandom* rnd;
    // The parameters we'll read from the .fcl file.
    art::InputTag fCRTDAQProducerLabel;

    // Other variables that will be shared between different methods.
    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    int                      fTriggerOffset;     ///< (units of ticks) time of expected neutrino event
    CRTCommonUtils* fCrtutils;  
  }; // class CRTCalibrationAnalysis


  //-----------------------------------------------------------------------
   
  CRTCalibrationAnalysis::CRTCalibrationAnalysis(Parameters const& config)
    : EDAnalyzer(config)
    , fCRTDAQProducerLabel(config().CRTDAQLabel())
    , fCrtutils(new CRTCommonUtils())
  {
    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();
    // The same for detector TDC clock services.
    // Access to detector properties.
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    fTriggerOffset = trigger_offset(clockData);

    for(int feb=1; feb<232; feb++){
      
 			string hname = "hadc_"+to_string(feb);
			string htitle = "raw charge: mac5 "+to_string(feb);
			string signame = hname + "_signal_sum";
			string resname = hname + "_reset_sum";
			pedSigHistos[feb] = tfs->make<TH1F>(signame.c_str(), htitle.c_str(), 131040, 0, 131040);
			pedResHistos[feb] = tfs->make<TH1F>(resname.c_str(), htitle.c_str(), 131040, 0, 131040);
      macToHistos[feb] = new vector<TH1F*>();
			macToHistos_nonTriggeringPed[feb] = new vector<TH1F*>();
      macToHistos_ped[feb] = new vector<TH1F*>();
      macToHistos_sig[feb] = new vector<TH1F*>();

      for(int ch=0; ch<32; ch++){

				string chname = hname + "_"+to_string(ch);
  			string chtitle = htitle + ", ch. "+to_string(ch);
				string ch_signame = chname + "_signal";
				string ch_pedname = chname + "_pedestal";
				string ch_pedname_nonTriggering = ch_pedname + "_non_triggering";
				macToHistos[feb]->push_back(tfs->make<TH1F>(chname.c_str(),chtitle.c_str(),4100,0,4100));
				macToHistos_nonTriggeringPed[feb]->push_back(tfs->make<TH1F>(ch_pedname_nonTriggering.c_str(),ch_pedname_nonTriggering.c_str(),4100,0,4100));
				macToHistos_ped[feb]->push_back(tfs->make<TH1F>(ch_pedname.c_str(),ch_pedname.c_str(),4100,0,4100));
				macToHistos_sig[feb]->push_back(tfs->make<TH1F>(ch_signame.c_str(),ch_signame.c_str(),4100,0,4100));
//				if(feb > 100) brokenChannels.push_back(ChannelInfo(feb, ch, chname));
      }
    }
    rnd = new TRandom();
  }
  
  //-----------------------------------------------------------------------
  void CRTCalibrationAnalysis::beginJob()
  {
  }
   
  void CRTCalibrationAnalysis::beginRun(const art::Run& /*run*/)
  {
  }
  
  void CRTCalibrationAnalysis::endRun(const art::Run& /*run*/)
  {
  }

   void CRTCalibrationAnalysis::endJob()
  {
  }

  //-----------------------------------------------------------------------

	//Function to have an updated list of broken channels
/*	void listOfBrokenChannels(art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle, std::vector<ChannelInfo>& channelsList){
	for(auto const& ch : channelsList) std::cout << "FEB: " << ch.feb << " CHANNEL: " << ch.channel << " IS IN THE LIST" << std::endl;
	int found_channels = 0;
    channelsList.erase(
        std::remove_if(channelsList.begin(), channelsList.end(), [&crtDAQHandle, &found_channels](const ChannelInfo& ch) {
	                for (const auto& febdat : *crtDAQHandle) {
										if(febdat.fMac5 <100) continue;
	                    if (febdat.fAdc[ch.channel] > 500) { //Does this logic make sense? NO
	found_channels++;
	                        return true;  
	                    }
	                 } 
									 return false;  
									}	), channelsList.end());
	std::cout << "I'M IN THE FUNCTION" << std::endl;
    for(auto const& ch : channelsList) std::cout << "FEB: " << ch.feb << " CHANNEL: " << ch.channel << " IS BROKEN!" << std::endl;
	}
*/

  // Function we need and we don't know where to put
  ULong64_t getTimeOfReset(art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle, uint flag_value) {
	  ULong64_t time = 0;
    //int current_feb = 0;
    for ( auto const& febdat : (*crtDAQHandle) ) {
      if (febdat.fFlags == flag_value) {
       // if (current_feb == febdat.fMac5) {
       //   continue;
       // }
        std::cout << "Found a reset hit on feb " << +febdat.fMac5  
                  << " at the time " << febdat.fTs0  << std::endl;
	      time = febdat.fTs0;

	// current_feb = febdat.fMac5;
        break;
      }
    }
	return time;
  }

vector<icarus::crt::CRTData> fixFlags(vector<icarus::crt::CRTData> crtDAQHandle, uint flag_value, ULong64_t time) {
	vector<icarus::crt::CRTData> new_data;
  for ( auto& febdat : crtDAQHandle ) {
    if (febdat.fMac5 > 100) {
      new_data.push_back(febdat);
      if (new_data.back().fFlags == 3) {
        Long64_t diff = abs((Long64_t)new_data.back().fTs0 - (Long64_t)time);
        if (diff < 100) {
          new_data.back().fFlags = flag_value;
        }
      }
    }
  }
  return new_data;
}

vector<icarus::crt::CRTData> fixFlags(art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle, uint flag_value, ULong64_t time) {
	vector<icarus::crt::CRTData> new_data;
  for ( auto& febdat : (*crtDAQHandle) ) {
    new_data.push_back(febdat);
  }
  return fixFlags(new_data, flag_value, time);
}




  //-----------------------------------------------------------------------
  void CRTCalibrationAnalysis::analyze(const art::Event& event) 
  {
    MF_LOG_DEBUG("CRTCalibrationAnalysis") << "beginning analyis" << '\n';


    art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    bool isCRTDAQ = event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle);

		//I check the state of the channels to make a list of broken ones		
		//listOfBrokenChannels(crtDAQHandle, brokenChannels);

		//For all the febs with mac5 lower than 100 I don't make any difference
    if (isCRTDAQ)  {
      MF_LOG_DEBUG("CRTCalibrationAnalysis") << "about to loop over CRTDAQ entries" << '\n';
    	int adc_sum;
			for ( auto const& febdat : (*crtDAQHandle) ) {
      	if (febdat.fMac5 < 100) {
					adc_sum=0;
					for(int ch=0; ch<32; ch++) {
						adc_sum+=febdat.fAdc[ch];
  	    	  macToHistos[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
  	    	  macToHistos_nonTriggeringPed[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
  	    	  macToHistos_ped[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
  	    	  macToHistos_sig[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
        	}
					pedSigHistos[febdat.fMac5]->Fill(adc_sum);
					pedResHistos[febdat.fMac5]->Fill(adc_sum);
      	}	
    	}

/*
     // Search of first resert hit and storing of its time
     ULong64_t time9 = getTimeOfReset(crtDAQHandle, 9);
     ULong64_t time7 = getTimeOfReset(crtDAQHandle, 7);

     // Fixing of unflagged reset hits
     vector<icarus::crt::CRTData> new_data = fixFlags(crtDAQHandle, 9, time9);
                                  new_data = fixFlags(new_data,     7, time7);
     
*/
     // Actual stuff with hits of the TopCRT (you'll later have to uncomment the fix on the hit)
     int c = 0; 
     int current_feb = 0;
//		 int pedChannels = 10;//Lower values per layer in the non-triggering channel logic for the pedestal
     for ( auto const& febdat : (*crtDAQHandle)/*new_data*/ ) {
        if (febdat.fMac5 > 100) {
					adc_sum=0;
	
	      	if (febdat.fFlags == 3) {
          	std::map<int, int> top_layer;

          	for (int ch = 0; ch < 16; ch++) {
  						top_layer[febdat.fAdc[ch]] = ch; //adc, ch with ch=[0,15]
							adc_sum+=febdat.fAdc[ch];
          	}

	        	auto max = std::max_element(top_layer.begin(), top_layer.end());
	    			macToHistos[febdat.fMac5]->at(max->second)->Fill( max->first );
	    			if(max->first > 275) macToHistos_sig[febdat.fMac5]->at(max->second)->Fill( max->first );	
	    			top_layer.erase(max);

		    		max = std::max_element(top_layer.begin(), top_layer.end());
		    		macToHistos[febdat.fMac5]->at(max->second)->Fill( max->first );
		    		if(max->first > 275) macToHistos_sig[febdat.fMac5]->at(max->second)->Fill( max->first );	
	          top_layer.erase(max);
					
          	std::map<int, int> bot_layer;

	    	  	for (int ch = 16; ch < 32; ch++) {
  						bot_layer[febdat.fAdc[ch]] = ch;
							adc_sum+=febdat.fAdc[ch];
          	}

	    			max = std::max_element(bot_layer.begin(), bot_layer.end());
	    			macToHistos[febdat.fMac5]->at(max->second)->Fill( max->first );
	    			if(max->first > 275) macToHistos_sig[febdat.fMac5]->at(max->second)->Fill( max->first );	
	    			bot_layer.erase(max);

	    			max = std::max_element(bot_layer.begin(), bot_layer.end());
	    			macToHistos[febdat.fMac5]->at(max->second)->Fill( max->first );
	    			if(max->first > 275) macToHistos_sig[febdat.fMac5]->at(max->second)->Fill( max->first );	
						bot_layer.erase(max);


					//In total I skip the higher 6 values per layer 4+2 erased
						max = std::max_element(top_layer.begin(), top_layer.end());
	    			top_layer.erase(max);
						max = std::max_element(top_layer.begin(), top_layer.end());
	    			top_layer.erase(max);
						max = std::max_element(top_layer.begin(), top_layer.end());
	    			top_layer.erase(max);
						max = std::max_element(top_layer.begin(), top_layer.end());
	    			top_layer.erase(max);

	    			max = std::max_element(bot_layer.begin(), bot_layer.end());	
						bot_layer.erase(max);
	    			max = std::max_element(bot_layer.begin(), bot_layer.end());
						bot_layer.erase(max);
	    			max = std::max_element(bot_layer.begin(), bot_layer.end());
						bot_layer.erase(max);
	    			max = std::max_element(bot_layer.begin(), bot_layer.end());	
						bot_layer.erase(max);

  					for (auto it = top_layer.begin(); it != top_layer.end(); ++it) {
    					macToHistos_nonTriggeringPed[febdat.fMac5]->at(it->second)->Fill(it->first);
  					}
  					for (auto it = bot_layer.begin(); it != bot_layer.end(); ++it) {
    					macToHistos_nonTriggeringPed[febdat.fMac5]->at(it->second)->Fill(it->first);
  					}


					//Instead of erasing here you tryed with a counter on the number of channels for the pedestal = last 10
/*					
  					int counter = 0;
  					for (auto it = top_layer.begin(); it != top_layer.end() && counter < pedChannels; ++it, ++counter) {
    					macToHistos_nonTriggeringPed[febdat.fMac5]->at(it->second)->Fill(it->first);
  					}
  
  					counter = 0; // Resetting counter for the bottom layer
  					for (auto it = bot_layer.begin(); it != bot_layer.end() && counter < pedChannels; ++it, ++counter) {
    				macToHistos_nonTriggeringPed[febdat.fMac5]->at(it->second)->Fill(it->first);
  }
*/
				
						pedSigHistos[febdat.fMac5]->Fill(adc_sum);

        }//SIGNAL HIT


  	    if (febdat.fFlags == 9 || febdat.fFlags == 7) {
		      c++;
	  //The first time it enters in the if() this isn't true, but if for the next hit the condition is true, then you have two or more hits flagged as reset in the same event. 
	  // Do you want to keep this check???
          if (current_feb == febdat.fMac5) {
						continue;
	      	}

          for(int ch=0; ch<32; ch++) {
  	        macToHistos[febdat.fMac5]->at(ch)->Fill(febdat.fAdc[ch]);
	    			macToHistos_ped[febdat.fMac5]->at(ch)->Fill(febdat.fAdc[ch]);	
						adc_sum+=febdat.fAdc[ch];
           }
					
					pedResHistos[febdat.fMac5]->Fill(adc_sum);
          current_feb = febdat.fMac5;
        }//RESET HIT

      }
    }

		
    std::cout << "HERE IS THE NUMBER OF RESET HITS: " << c << std::endl;
  }//if crtdetsim products present

    else 
      mf::LogError("CRTCalibrationAnalysis") << "CRTDAQ products not found!" << std::endl;
    
  } // CRTCalibrationAnalysis::analyze()
  
  DEFINE_ART_MODULE(CRTCalibrationAnalysis)
} // namespace crt
} // namespace icarus

