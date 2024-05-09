/**
 * CRT Calibration Analysis Module
 * 
 * This file contains the CRT Calibration Analysis module to access CRT raw data and reco products.
 * The module reads CRT data and processed outputs to perform calibration analysis and diagnostics.
 * 
 * It was originally authored by Chris Hilgenberg and has been updated and maintained to include more features.
 * Specifically part of the initialization code, including includes and the definition of the class 
 * and methods, was implemented by Chris Hilgenberg (Chris.Hilgenberg@colostate.edu).
 * It was last revised in October 2018 with LArSoft v07_06_01. 
 * 
 * The implementation details and the bulk of the logic, were developed by me.
 */

// LArSoft includes for handling CRT data and geometry within the LArSoft framework.
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ART framework core classes, for dealing with event data and services
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries including message logging and parameter set handling
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "cetlib/pow.h" // cet::sum_of_squares()

// ROOT libraries used for histogramming and other data processing functions
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"

// Standard C++ libraries
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
	// Declaration of the CRTCalibrationAnalysis class, inheriting from art::EDAnalyzer
	class CRTCalibrationAnalysis : public art::EDAnalyzer
  	{
  	public:
		// Config struct, to keep parameters from FHiCL files
    		struct Config {
      
      			using Name = fhicl::Name;
      			using Comment = fhicl::Comment;
      
      			fhicl::Atom<art::InputTag> CRTDAQLabel {
        			Name("CRTDAQLabel"),
        			Comment("tag of the input data product with calibrated CRT data")
        		};
   

    		}; // Config

    		using Parameters = art::EDAnalyzer::Table<Config>;
    
  
    		// Constructor: Initializes the CRTCalibrationAnalysis module using configuration parameters
    		explicit CRTCalibrationAnalysis(Parameters const& config);

		// Called once at the start of the job; use to initialize resources
    		virtual void beginJob() override;

		// Called at the beginning of each run, use for run-specific initialization
    		virtual void beginRun(const art::Run& run) override;

		// Processes each event, main analysis logic is implemented here
    		virtual void analyze (const art::Event& event) override;

		// Called once at the end of the job, use for cleanup and finalizing results
    		void endJob() override;

		// Called at the end of each run, use for run-specific cleanup
    		void endRun(const art::Run&) override;

  	private:
		// Handle to ROOT file service
    		art::ServiceHandle<art::TFileService> tfs;

		map<uint8_t, TH1F*> allChannels_adcSum_histograms; // Histograms containing ADC sums over all channels, the map key is the FEB number
		map<uint8_t, TH1F*> allChannels_resetHits_adcSum_histograms; // Histograms for reset hits ADC sums, per FEB number
    		map<uint8_t,vector<TH1F*>*> channelSpectrum_histograms; // Channels' spectra histograms, the map key is the FEB number
		map<uint8_t,vector<TH1F*>*> channelSpectrum_pedestal_noTrig_histograms; // Channels' pedestal histograms, filled with the ADC counts of the channels that did not trigger the event
    		map<uint8_t,vector<TH1F*>*> channelSpectrum_pedestal_resetHits_histograms; // Channels' pedestal histograms, filled with the ADC counts of that channels when a 'reset hit' is found
    		map<uint8_t,vector<TH1F*>*> channelSpectrum_onlySignal_histograms; // Signal-only histograms
    		TRandom* rnd;

    		// // Label for CRT data product, with the parameters we'll read from the .fcl file.
    		art::InputTag fCRTDAQProducerLabel;

    		// Other variables that will be shared between different methods.
    		geo::GeometryCore const* fGeometryService;   // Pointer to Geometry provider
    		int                      fTriggerOffset;     // Trigger offset in units of ticks (time of expected neutrino event)
    		CRTCommonUtils* fCrtutils;  // CRT utilities
  	}; // class CRTCalibrationAnalysis

//-----------------------------------------------------------------------
   	// The actual constructor implementation
	CRTCalibrationAnalysis::CRTCalibrationAnalysis(Parameters const& config)
    	: EDAnalyzer(config)
    	, fCRTDAQProducerLabel(config().CRTDAQLabel())
    	, fCrtutils(new CRTCommonUtils())
 	{
    		// Obtain a pointer to the geometry service, which provides access to detector geometry
    		fGeometryService = lar::providerFrom<geo::Geometry>();
    		
		// Obtain detector-specific time information from the DetectorClocksService
    		auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    		fTriggerOffset = trigger_offset(clockData);

    		for(int feb=1; feb<232; feb++){
      			// For each FEB generate the histograms names
 			string hname = "hadc_"+to_string(feb);
			string htitle = "raw charge: mac5 "+to_string(feb);
			string signame = hname + "_signal_sum";
			string resname = hname + "_reset_sum";

			// Create the histograms 
			allChannels_adcSum_histograms[feb] = tfs->make<TH1F>(signame.c_str(), htitle.c_str(), 131040, 0, 131040);
			allChannels_resetHits_adcSum_histograms[feb] = tfs->make<TH1F>(resname.c_str(), htitle.c_str(), 131040, 0, 131040);
      			channelSpectrum_histograms[feb] = new vector<TH1F*>(); // The vector will contain the 32 histograms, one per channel
			channelSpectrum_pedestal_noTrig_histograms[feb] = new vector<TH1F*>();
      			channelSpectrum_pedestal_resetHits_histograms[feb] = new vector<TH1F*>();
      			channelSpectrum_onlySignal_histograms[feb] = new vector<TH1F*>();

      			for(int ch=0; ch<32; ch++){
				// Append the channel number in the histograms names
				string chname = hname + "_"+to_string(ch);
  				string chtitle = htitle + ", ch. "+to_string(ch);
				// Append the histogram type: signal or pedestal (noise)
				string ch_signame = chname + "_signal";
				string ch_pedname = chname + "_pedestal"; // The default 'pedestal' histograms will contain the spectrum obtained by the reset hits
				string ch_pedname_nonTriggering = ch_pedname + "_non_triggering"; // The 'pedestal' histograms obtained with a non-triggering channels logic have a specific name

				// Allocate histograms for each Front-End Board (FEB) and each channel
				channelSpectrum_histograms[feb]->push_back(tfs->make<TH1F>(chname.c_str(),chtitle.c_str(),4100,0,4100));
				channelSpectrum_pedestal_noTrig_histograms[feb]->push_back(tfs->make<TH1F>(ch_pedname_nonTriggering.c_str(),ch_pedname_nonTriggering.c_str(),4100,0,4100));
				channelSpectrum_pedestal_resetHits_histograms[feb]->push_back(tfs->make<TH1F>(ch_pedname.c_str(),ch_pedname.c_str(),4100,0,4100));
				channelSpectrum_onlySignal_histograms[feb]->push_back(tfs->make<TH1F>(ch_signame.c_str(),ch_signame.c_str(),4100,0,4100));
      			}
    		}
		
		rnd = new TRandom();
  	}
  
//-----------------------------------------------------------------------
    	// Begin job lifecycle methods to setup analysis
  	void CRTCalibrationAnalysis::beginJob()
  	{
  	}
   
  	void CRTCalibrationAnalysis::beginRun(const art::Run& /*run*/)
  	{
  	}

  	// Handle end of run and job cleanup
  	void CRTCalibrationAnalysis::endRun(const art::Run& /*run*/)
  	{
  	}

  	void CRTCalibrationAnalysis::endJob()
  	{
  	}
//------------------------------------------------------------------------------------------------------
	// 'analyze' method processes each event, filling histograms and performing analysis
  	void CRTCalibrationAnalysis::analyze(const art::Event& event) 
  	{
    		MF_LOG_DEBUG("CRTCalibrationAnalysis") << "beginning analyis" << '\n';

		art::Handle<vector<icarus::crt::CRTData>> crtDAQHandle;
    		bool isCRTDAQ = event.getByLabel(fCRTDAQProducerLabel, crtDAQHandle);

		// Proceed if CRT data is found in the event
    		if (isCRTDAQ)  {
      			MF_LOG_DEBUG("CRTCalibrationAnalysis") << "about to loop over CRTDAQ entries" << '\n';
    			int adc_sum;

			// Loop through each CRT data entry for FEBs with mac5 < 108 (side CRT)
			for ( auto const& febdat : (*crtDAQHandle) ) {
				//For all the febs with mac5 lower than 108 (Side CRT) I don't make any difference and only fill them with the ADC counts measured per channel
      				if (febdat.fMac5 < 108) {
					adc_sum=0; // Reset ADC sum for the new FEB data
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

			// Counter for the number of reset hits (signal) found in an event, as sometimes we have more than one
     			int resetHits_counter = 0; 
     			int current_feb = 0;

			// Processing events at the Top CRT (with a FEB MAC5 ID greater than 107)
     			for ( auto const& febdat : (*crtDAQHandle)) {
        			if (febdat.fMac5 > 107) {
					adc_sum=0; // Reset ADC sum for the new FEB data
					// Check if the data corresponds to a signal hit based on flag value
	      				if (febdat.fFlags == 3) {
          					std::map<int, int> top_layer; // Map to hold the module's top layer channel data: ADC values and their corresponding channels
						
						// Populate the map for the top layer channels (0-15)
						for (int ch = 0; ch < 16; ch++) {
  							top_layer[febdat.fAdc[ch]] = ch; // Map ADC value to channel number
							adc_sum+=febdat.fAdc[ch]; // Map ADC value to channel number
          					}

						// Find and fill histograms with the first maximum ADC value from the top layer
	        				auto max = std::max_element(top_layer.begin(), top_layer.end());
	    					channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
	    					if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first ); // Fill the channel's 'signal only' histogram for ADC counts values > 275  
	    					top_layer.erase(max);

						// Repeat for the second maximum ADC value
		    				max = std::max_element(top_layer.begin(), top_layer.end());
		    				channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
		    				if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
	          				top_layer.erase(max);

						// Similar processing for the bottom layer channels (16-31)
          					std::map<int, int> bot_layer;

						// Populate the map for the bottom layer channels
	    	  				for (int ch = 16; ch < 32; ch++) {
  							bot_layer[febdat.fAdc[ch]] = ch;
							adc_sum+=febdat.fAdc[ch];
          					}

						// Process the maximum values in the bottom layer similarly
	    					max = std::max_element(bot_layer.begin(), bot_layer.end());
	    					channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
	    					if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
	    					bot_layer.erase(max);

	    					max = std::max_element(bot_layer.begin(), bot_layer.end());
	    					channelSpectrum_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );
	    					if(max->first > 275) channelSpectrum_onlySignal_histograms[febdat.fMac5]->at(max->second)->Fill( max->first );	
						bot_layer.erase(max);


						// In total I skip the higher 6 values per layer
						// Top layer
						max = std::max_element(top_layer.begin(), top_layer.end());
	    					top_layer.erase(max);
						max = std::max_element(top_layer.begin(), top_layer.end());
	    					top_layer.erase(max);
						max = std::max_element(top_layer.begin(), top_layer.end());
	    					top_layer.erase(max);
						max = std::max_element(top_layer.begin(), top_layer.end());
	    					top_layer.erase(max);

						// Bottom layer
	    					max = std::max_element(bot_layer.begin(), bot_layer.end());	
						bot_layer.erase(max);
	    					max = std::max_element(bot_layer.begin(), bot_layer.end());
						bot_layer.erase(max);
	    					max = std::max_element(bot_layer.begin(), bot_layer.end());
						bot_layer.erase(max);
	    					max = std::max_element(bot_layer.begin(), bot_layer.end());	
						bot_layer.erase(max);

						// The lower 10 ADC values populate the 'pedestal' histograms, in a 'non-triggering channels' logic 
						// In the top layer
  						for (auto it = top_layer.begin(); it != top_layer.end(); ++it) {
    							channelSpectrum_pedestal_noTrig_histograms[febdat.fMac5]->at(it->second)->Fill(it->first);
  						}

						// Same for bottom layer
  						for (auto it = bot_layer.begin(); it != bot_layer.end(); ++it) {
    							channelSpectrum_pedestal_noTrig_histograms[febdat.fMac5]->at(it->second)->Fill(it->first);
  						}

						// Fill the histogram containing the sum of ADC count values of all the channels in this FEB
						allChannels_adcSum_histograms[febdat.fMac5]->Fill(adc_sum);

        				}//SIGNAL HIT

					// Handle reset hits identified by specific flags, to populate the pedestal
  	    				if (febdat.fFlags == 9 || febdat.fFlags == 7) {
		      				resetHits_counter++; // Increment counter for reset hits
						
	  					// Avoid duplicate processing of the same FEB within a single event, useful if there are more than one reset hit in this event 
          					if (current_feb == febdat.fMac5) {
							continue;
	      					}

						// Process each channel for reset hits
          					for(int ch=0; ch<32; ch++) {
  	        					channelSpectrum_histograms[febdat.fMac5]->at(ch)->Fill(febdat.fAdc[ch]);
	    						channelSpectrum_pedestal_resetHits_histograms[febdat.fMac5]->at(ch)->Fill(febdat.fAdc[ch]); // Fill the other 'pedestal' histogram, to compare the two logics	
							adc_sum+=febdat.fAdc[ch];
          			 		}
					
						allChannels_resetHits_adcSum_histograms[febdat.fMac5]->Fill(adc_sum);
          					current_feb = febdat.fMac5; // Update the current FEB to avoid re-processing
                    			}//RESET HIT
        			}

      			}

			// Print the number of reset hits to check if and when more that one are found, in a specific event
    			std::cout << "HERE IS THE NUMBER OF RESET HITS: " << resetHits_counter << std::endl;
    		
		}else{ 
      			mf::LogError("CRTCalibrationAnalysis") << "CRTDAQ products not found!" << std::endl;
		}//if crtdetsim products present
		
  	} // CRTCalibrationAnalysis::analyze()
  
  	DEFINE_ART_MODULE(CRTCalibrationAnalysis)
} // namespace crt
} // namespace icarus

