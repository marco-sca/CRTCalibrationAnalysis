/**
 * CRT Calibration Analysis Module
 * 
 * This file contains the CRT Calibration Analysis module to access CRT data and reco products.
 * 
 * Part of the initialization code, including includes and the definition of the class 
 * and methods, was implemented by Chris Hilgenberg (Chris.Hilgenberg@colostate.edu).
 * It was last revised in October 2018 with LArSoft v07_06_01.
 * 
 * The implementation details and the bulk of the logic, were developed by me.
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

//-----------------------------------------------------------------------
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
    
  
    		// Constructor: configures the module (see the Config structure above)
    		explicit CRTCalibrationAnalysis(Parameters const& config);

    		virtual void beginJob() override;
    		virtual void beginRun(const art::Run& run) override;
    		virtual void analyze (const art::Event& event) override;
    		void endJob() override;
    		void endRun(const art::Run&) override;

  	private:

    		art::ServiceHandle<art::TFileService> tfs;

		map<uint8_t, TH1F*> allChannels_adcSum_histograms;
		map<uint8_t, TH1F*> allChannels_resetHits_adcSum_histograms;
    		map<uint8_t,vector<TH1F*>*> channelSpectrum_histograms;
		map<uint8_t,vector<TH1F*>*> channelSpectrum_pedestal_noTrig_histograms;
    		map<uint8_t,vector<TH1F*>*> channelSpectrum_pedestal_resetHits_histograms;
    		map<uint8_t,vector<TH1F*>*> channelSpectrum_onlySignal_histograms;
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
			allChannels_adcSum_histograms[feb] = tfs->make<TH1F>(signame.c_str(), htitle.c_str(), 131040, 0, 131040);
			allChannels_resetHits_adcSum_histograms[feb] = tfs->make<TH1F>(resname.c_str(), htitle.c_str(), 131040, 0, 131040);
      			channelSpectrum_histograms[feb] = new vector<TH1F*>();
			channelSpectrum_pedestal_noTrig_histograms[feb] = new vector<TH1F*>();
      			channelSpectrum_pedestal_resetHits_histograms[feb] = new vector<TH1F*>();
      			channelSpectrum_onlySignal_histograms[feb] = new vector<TH1F*>();

      			for(int ch=0; ch<32; ch++){

				string chname = hname + "_"+to_string(ch);
  				string chtitle = htitle + ", ch. "+to_string(ch);
				string ch_signame = chname + "_signal";
				string ch_pedname = chname + "_pedestal";
				string ch_pedname_nonTriggering = ch_pedname + "_non_triggering";
				channelSpectrum_histograms[feb]->push_back(tfs->make<TH1F>(chname.c_str(),chtitle.c_str(),4100,0,4100));
				channelSpectrum_pedestal_noTrig_histograms[feb]->push_back(tfs->make<TH1F>(ch_pedname_nonTriggering.c_str(),ch_pedname_nonTriggering.c_str(),4100,0,4100));
				channelSpectrum_pedestal_resetHits_histograms[feb]->push_back(tfs->make<TH1F>(ch_pedname.c_str(),ch_pedname.c_str(),4100,0,4100));
				channelSpectrum_onlySignal_histograms[feb]->push_back(tfs->make<TH1F>(ch_signame.c_str(),ch_signame.c_str(),4100,0,4100));
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
//------------------------------------------------------------------------------------------------------
  	void CRTCalibrationAnalysis::analyze(const art::Event& event) 
  	{
    		MF_LOG_DEBUG("CRTCalibrationAnalysis") << "beginning analyis" << '\n';

		art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    		bool isCRTDAQ = event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle);

		//For all the febs with mac5 lower than 100 I don't make any difference
    		if (isCRTDAQ)  {
      			MF_LOG_DEBUG("CRTCalibrationAnalysis") << "about to loop over CRTDAQ entries" << '\n';
    			int adc_sum;
			
			for ( auto const& febdat : (*crtDAQHandle) ) {
      				if (febdat.fMac5 < 100) {
					adc_sum=0;
					for(int ch=0; ch<32; ch++) {
						adc_sum+=febdat.fAdc[ch];
  	    	  				channelSpectrum_histograms[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
  	    	  				channelSpectrum_pedestal_noTrig_histograms[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
  	    	  				channelSpectrum_pedestal_resetHits_histograms[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
  	    	  				channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(ch)->Fill( febdat.fAdc[ch] );
        				}
					allChannels_adcSum_histograms[febdat.fMac5]->Fill(adc_sum);
					allChannels_resetHits_adcSum_histograms[febdat.fMac5]->Fill(adc_sum);
      				}	
    			}
			
     			int c = 0; 
     			int current_feb = 0;
			
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
	    					channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
	    					if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
	    					top_layer.erase(max);

		    				max = std::max_element(top_layer.begin(), top_layer.end());
		    				channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
		    				if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
	          				top_layer.erase(max);
					
          					std::map<int, int> bot_layer;

	    	  				for (int ch = 16; ch < 32; ch++) {
  							bot_layer[febdat.fAdc[ch]] = ch;
							adc_sum+=febdat.fAdc[ch];
          					}

	    					max = std::max_element(bot_layer.begin(), bot_layer.end());
	    					channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
	    					if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
	    					bot_layer.erase(max);

	    					max = std::max_element(bot_layer.begin(), bot_layer.end());
	    					channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
	    					if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
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
    							channelSpectrum_pedestal_noTrig_histograms[febdat.fMac5]->at(it->second)->Fill(it->first);
  						}
						
  						for (auto it = bot_layer.begin(); it != bot_layer.end(); ++it) {
    							channelSpectrum_pedestal_noTrig_histograms[febdat.fMac5]->at(it->second)->Fill(it->first);
  						}


					//Instead of erasing here you tryed with a counter on the number of channels for the pedestal = last 10

				
						allChannels_adcSum_histograms[febdat.fMac5]->Fill(adc_sum);

        }//SIGNAL HIT


  	    if (febdat.fFlags == 9 || febdat.fFlags == 7) {
		      c++;
	  //The first time it enters in the if() this isn't true, but if for the next hit the condition is true, then you have two or more hits flagged as reset in the same event. 
          if (current_feb == febdat.fMac5) {
						continue;
	      	}

          for(int ch=0; ch<32; ch++) {
  	        channelSpectrum_histograms[febdat.fMac5]->at(ch)->Fill(febdat.fAdc[ch]);
	    			channelSpectrum_pedestal_resetHits_histograms[febdat.fMac5]->at(ch)->Fill(febdat.fAdc[ch]);	
						adc_sum+=febdat.fAdc[ch];
           }
					
					allChannels_resetHits_adcSum_histograms[febdat.fMac5]->Fill(adc_sum);
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

