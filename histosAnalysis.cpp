// ROOT includes
#include "TH1F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TObjArray.h"
#include "TObjString.h"

// C++ includes
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//Define here the hard values for the fit
bool save_png = true;
const int chi2_red_pedestal = 10;
const int chi2_red_signal = 2;
const int pedestal_rebin = 5;
const int pedestal_smooth = 0;
const int signal_rebin = 10;
const int signal_smooth = 2;
const int adc_last_pe = 1000;
const int range_pedestal_histo = 500;
const int chi2_fit_increase = 10;
const int chi2_red_limit = 500;
const int max_peaks = 15;
const int peaks_of_interest = 5;
const int max_nr_of_fits = 10;
const int rearrangement_tries = 2;
const double sigmaXPeakSearch = 1;
const double peaks_threshold = 0.10;
std::string dataset_specific_name = "afterCut";

//Contains gain information of a channel, including pedestal values, slope, peaks, and associated errors.
struct gain {
	double pedestal;
  	double pedestal_error;

  	double slope;
  	double slope_error;

  	std::vector<double> peaks;
  	std::vector<double> peaks_error;

  	std::vector<int> index;

  	std::string feb;
  	std::string channel;
};

//Represents the result of a linear fit, containing the slope (channel's gain), error, and reduced chi-squared value.
struct lFit {
	double slope;
	double error;
	double red_chi2;
};

//--------------------------------------------------------------------------------------------------------------------------
/**
 * This function creates 2 separate files for the pedestal and the signal histograms
 *
 * INPUT PARAMETERS
 * fIn: Input file containing all the histograms generated after the decoding stage;
 * fInSig: File that will contain only the histograms representing the integrated adc spectrum generated by the signal 
 * fInPed: File that will contain only the histograms representing the pedestal of each channel of a CRT module
 *
 * WORKFLOW
 * The ROOT file produced at the decoding stage contains a directory named "CRTCalibrationAnalysis";
 * inside we can find directly the histograms or some subdirectories containing the histograms, this depends on how
 * many histograms are created during the decoding of raw data;
 * I search for the names of the histograms, that have the following format "hadc_channelNumber_febNumber_type", where the
 * type can be "signal" or "pedestal";
 * once the histogram is found is saved in the respective file: fInSig if is of type signal, fInPed if is of type pedestal.
 *
 * */
void separateHistograms(const char* fIn, TFile *fInSig, TFile *fInPed) {
	TFile *inputFile = TFile::Open(fIn, "READ");

	TDirectory *mainDir = (TDirectory*)inputFile->Get("CRTCalibrationAnalysis");
	TList *list = mainDir->GetListOfKeys();
	bool hasSubDirs = false;
	
	// Here I search for subdirectories inside CRTCalibrationAnalysis
    	TIter iter(list);
    	TKey *key;
    	while ((key = (TKey*)iter())) {
    		if (strcmp(key->GetClassName(), "TDirectoryFile") == 0) {
        		hasSubDirs = true;
            		break;
        	}
    	}

	// I loop over and search all the histograms for the 32 channels in the 231 CRT modules/Front-End-Boards (FEB) 
    	for (int feb = 1; feb <= 231; ++feb) {
    		for (int ch = 0; ch <= 31; ++ch) {
        		TString signalHistName = Form("hadc_%d_%d_signal", feb, ch);
            		TString pedestalHistName = Form("hadc_%d_%d_pedestal", feb, ch);
            		TH1 *signalHist = nullptr;
            		TH1 *pedestalHist = nullptr;

			// If a subdirectory is found, I point at it and search for the histogram using its name
            		if (hasSubDirs) {
            			TIter dirIter(list);		
				while ((key = (TKey*)dirIter())) {
					// Double-check that I am pointing to a subdirectory
					if (strcmp(key->GetClassName(), "TDirectoryFile") != 0) {
           	        			continue;
                  			}
                  			TDirectory *subDir = (TDirectory*)mainDir->Get(key->GetName());
                  			signalHist = (TH1*)subDir->Get(signalHistName);
                  			pedestalHist = (TH1*)subDir->Get(pedestalHistName);
                  			if (signalHist || pedestalHist) {
						std::cout << "Failed search of: " << signalHistName << " or " << pedestalHistName << std::endl;
                    				break;
                    			}
				}	
			// Otherwise I search directly for the histogram using its name
            		} else {
            			signalHist = (TH1*)mainDir->Get(signalHistName);
            			pedestalHist = (TH1*)mainDir->Get(pedestalHistName);
            		}

			// I then move the histogram in the right file depending on the type: signal or pedestal
            		if(signalHist) {
				std::cout << "Moving " << signalHistName << " to " << fInSig << std::endl;
                		fInSig->cd();
                		signalHist->Write();
            		}

            		if (pedestalHist) {
				std::cout << "Moving " << pedestalHistName << " to " << fInPed << std::endl;
          		  	fInPed->cd();
                		pedestalHist->Write();
            		}
        	}
   	 }
	inputFile->Close();
}

/**
 * This function smooths and rebins the histogram, eliminating also the underflow and overflow
 *
 * INPUT PARAMETERS
 * histogram: The histogram of type TH1I (1D histogram with integer bin content)
 * rebin: Number of bins to merge together in the new histogram
 * smooth: Number of times the smoothing procedure is repeated
 *
 * WORKFLOW
 * The last bin of the histogram is set, using a "for" loop, to zero to solve the overflow problem.
 * Then the first bin is also set to 0, solving the underflow problem. 
 * Overflow/underflow: the first and last bin have very high values. 
 * The last bin (4095) contains all the collected signal with ADC counts > 4095.
 * Finally I use the Rebin() and Smooth() functions of ROOT, that help reducing the noise and statistical fluctuations
 * in the data. 
 *
 * */
void rebinAndSmooth(TH1I* histogram, int rebin, int smooth) {
		for ( int nbin = histogram->GetNbinsX(); nbin > 0; nbin--){
				if( histogram->GetBinContent(nbin) != 0) {
						histogram->SetBinContent(nbin, 0);
						break;
				}
		}

		histogram->SetBinContent(1,0);
		histogram->Rebin(rebin);
		histogram->Smooth(smooth);
}

/**
 * This function is used to extract the sigma and mean values for the peaks of a given histogram, using recursive gaussian fits
 * of the peaks
 *
 * INPUT PARAMETERS
 * histogram: histogram from which we will extract the peaks' mean and sigma values
 * xpeaks_sorted: vector containing the peaks' positions on the x-axis (in ADC counts)
 * requested_red_chi2: requested reduced chi squared value for the fit
 * histo_type: can assume 2 values, "ped" or "sig", depending on the histogram passed as argument
 *
 * WORKFLOW
 * Iterating through each peak x-axis position in the xpeaks_sorted vector, I calculate the left and right distance between peaks. 
 * For the first peak (pedestal) the left distance is equal to its distance from the origin of the x-axis. For the last peak and 
 * the pedestal, the right distance is computed using adc_last_pe.
 * These distances are used to determine the fit range, that depends on which distance (right or left) is smaller. 
 * The fit range is computed as the peaks distance divided by 2 for the signal (p.e) peaks and by 1 for the pedestal.
 * A fitting loop is used to search the best Gaussian fit in the specified range.
 * If the fit meets the following criteria:
 * - Reduced chi-squared is below the requested value
 * - Mean falls within the range of the peak
 * - Mean + sigma is below a certain threshold
 * add the mean and sigma to the "means" vector and exit the loop. Otherwise shrink the fit range and repeat the fit. 
 * If there are no more bins in the selected range or the fit was repeated more than "max_nr_of_fits" times the loop ends without 
 * finding a suitable fit. The current mean and sigma are then added to the "means" vector.
 *
 * RETURNS
 * The vector "means", containing the values of mean and sigma for the peaks found in the histogram
 *  
 * */
std::vector<std::pair<double, double>> getPeaksMean(TH1I* histogram,  std::vector<double> xpeaks_sorted, int requested_red_chi2, std::string histo_type) {
	std::vector<std::pair<double, double>> means;
	int low_binx, up_binx, den_fit_range, nr_of_fits;
	double left_distance, right_distance, half_mean_distance, low_x, up_x;
	double mean, sigma, red_chi2;
	bool fitted;

	// Loop over each peak x-axis position and compute the left and right distances between peaks	
	for(int npeak = 0; npeak < xpeaks_sorted.size(); npeak++){
    	TAxis* xaxis = histogram->GetXaxis();
    		
	if (npeak == 0) {
      		left_distance = xpeaks_sorted[npeak];
    	} else {
      		left_distance = xpeaks_sorted[npeak] - xpeaks_sorted[npeak - 1];
    	}

    	if (npeak == xpeaks_sorted.size() - 1) {
      		right_distance = adc_last_pe - xpeaks_sorted[npeak];
    	} else {
      		right_distance = xpeaks_sorted[npeak + 1] - xpeaks_sorted[npeak];
    	}

	// The fit range is computed as the smaller peaks' distance (right or left), divided by 2 for the signal (p.e) peaks and by 1 for the pedestal. For the signal the set distance is the adc counts distance between peaks (not true for the pedestal, as its spectrum has a single peak.
	if(histo_type == "ped") den_fit_range = 1;
	if(histo_type == "sig") den_fit_range = 2;
	if(left_distance < right_distance){
      		half_mean_distance = left_distance / den_fit_range;
    	} else {
      		half_mean_distance = right_distance / den_fit_range;
    	}

	low_binx = xaxis->FindBin(xpeaks_sorted[npeak] - half_mean_distance);
    	up_binx  = xaxis->FindBin(xpeaks_sorted[npeak] + half_mean_distance);

	fitted = false;
    	nr_of_fits = 0;
    	while (!fitted) {

		// Create a Gaussian fit function and fit it to the data within the specified range.
      		low_x = xaxis->GetBinCenter(low_binx);
      		up_x  = xaxis->GetBinCenter(up_binx);
      		TF1* gaussian_fit = new TF1("Gaussian Fit", "gaus", low_x, up_x);
      		gaussian_fit->SetLineColor(2);
      		histogram->Fit(gaussian_fit, "QR+");

		// Get the mean, sigma, and reduced chi-squared value from the Gaussian fit.
      		mean = gaussian_fit->GetParameter(1);
      		sigma = gaussian_fit->GetParameter(2);
      		red_chi2 = gaussian_fit->GetChisquare() / gaussian_fit->GetNDF();
	  		
		// Check if the fit meets the criteria.
		if (red_chi2 < requested_red_chi2 &&
          		mean > xpeaks_sorted[npeak] - half_mean_distance &&
          		mean < xpeaks_sorted[npeak] + half_mean_distance && 
          		mean + sigma < adc_last_pe) {
        			means.push_back(std::make_pair(mean, sigma));
        			fitted = true;
      		}

		// Increment/decrement bin indices to shrink the fit range, and check termination conditions.
		low_binx++;
		up_binx--;
		if (up_binx - low_binx < 2 || nr_of_fits >= max_nr_of_fits){ 
			fitted = true;
			means.push_back(std::make_pair(mean, sigma));
		}
			
		nr_of_fits++;
		delete gaussian_fit;
		} 
	}

	return means;
}

/**
 * The function is used to extract the ordered peak positions from a histogram
 *
 * INPUT PARAMETERS
 * histogram: pointer to the TH1I histogram from which peak positions will be extracted
 * npeaks: the desired number of peaks to detect
 *
 * WORKFLOW
 * The following variables initialized: nr_peaks_found to store the number of detected peaks, xpeaks for storing the 
 * x-positions of peaks, and xpeaks_sorted for the sorted peak positions.
 * A TSpectrum object is created to analyze the spectrum of the histogram with the specified number of peaks and set 
 * the x-axis range for peak search within the histogram.
 * TSpectrum is used to search for peaks, providing the parameters:
 *  - histogram: histogram to analyze;
 *  - sigmaXPeakSearch, width of the Gaussian kernel applied during the peak search;
 *  - "", empty string for optional parameters;
 *  - peaks_threshold, threshold value controlling which peaks are identified as significant (peaks with amplitudes 
 *     below this threshold are not reported).
 * The number of peaks found is then stored in nr_peaks_found and the x-positions of the detected peaks are extracted, 
 * copied in the xpeaks_sorted array and sorted.
 *
 * RETURNS
 * The vector xpeaks_sorted containing the x-axis values of the detected peaks in ascending order.
 *
 * */
std::vector<double> getPeaks(TH1I* histogram, int npeaks) {
	int nr_peaks_found;
	double* xpeaks;
	std::vector<double> xpeaks_sorted;

	// Create a TSpectrum analyzer with the desired number of peaks
	TSpectrum* spectrum = new TSpectrum(npeaks);

	// Search for peaks in the histogram, in a specified range, using the TSpectrum analyzer
  	histogram->GetXaxis()->SetRangeUser(0, adc_last_pe);
  	nr_peaks_found = spectrum->Search(histogram, sigmaXPeakSearch, "", peaks_threshold);

	// Get the x-positions of the detected peaks and copy them into a sorted vector
  	xpeaks = spectrum->GetPositionX();
  	for (int peak_nr = 0; peak_nr < nr_peaks_found; peak_nr++) {
    		xpeaks_sorted.push_back(xpeaks[peak_nr]);
  	}
  	std::sort(xpeaks_sorted.begin(), xpeaks_sorted.end());

  	delete spectrum;

  	return xpeaks_sorted;
}

/**
 * The function creates a graph with peak data (points) and associated sigma (error)
 *
 * INPUT PARAMETES
 * means: vector containing the mean and sigma values for peaks
 * indices: vector containing the indices associated to each peak (peak's number in the distribution)
 *
 * WORKFLOW
 * A TGraphErrors graph is created with the number of data points equal to the size of "means" vector.
 * A loop over each peak adds then data points to the graph, where the sigma of the fitted peaks are represented with error bars.
 *
 * RETURNS
 * The filled TGraphErrors "graph"
 * */
TGraphErrors fillGraph(std::vector<std::pair<double, double>> means, std::vector<int> indices) {
	TGraphErrors graph(means.size());

	// Loop through each peak mean and fill the graph with data points and sigmas (represented as error bars)
  	for (int peak_mean_nr = 0; peak_mean_nr < means.size(); peak_mean_nr++) {
    		graph.SetPoint(peak_mean_nr, indices[peak_mean_nr], means[peak_mean_nr].first);
    		
		// Set the error bars for the data points, where the x-axis error is set to 0
		graph.SetPointError(peak_mean_nr, 0, means[peak_mean_nr].second);
  	}

  	return graph;
}

/**
 * The function performs a linear fit on given data points
 *
 * INPUT PARAMETERS
 * means: vector of pairs containing the mean and sigma values of the gaussian fit on the peaks
 * indices: vector of integers containing the index related to each peak (peak number)
 *
 * WORKFLOW
 * A linear fit function of type "pol1" (first degree polynomial) is created, and a graph with error bars (sigma) is 
 * filled with data points (peak mean values).
 * The linear function is fitted to the data in the graph and a struct of type "lFit" is defined to store the fit 
 * results (slope, error, reduced chi-squared).
 *
 * RETURNS
 * An "lFit" struct containing the slope, error, and reduced chi-squared value for the linear fit.
 *
 * */
lFit linearFit(std::vector<std::pair<double, double>> means, std::vector<int> indices) {
  	// Create a linear (first degree polynomial) fit function named "lfit"
	TF1* linear_fit = new TF1("lfit", "pol1");

	// Create a TgraphErrors and fill it with peaks mean values (data points) and sigmas (error bars)
  	auto graph = fillGraph(means, indices);

	//Fit the linear function (blue color) to the data in the graph (Q: suppress graphic outputs)
  	linear_fit->SetLineColor(2);
  	graph.Fit(linear_fit, "Q");

	//Create an lFit struct and store the fit results inside it
  	lFit fit_info;
  	fit_info.slope = linear_fit->GetParameter(1);
  	fit_info.error = linear_fit->GetParError(1);
  	fit_info.red_chi2 = linear_fit->GetChisquare() / linear_fit->GetNDF();

  	delete linear_fit;
  	
	return fit_info;
}

/**
 * The function aims at adjusting the peak indices that might have been detected incorrectly,
 * in an attempt to improve the quality of the linear fit results. Peaks can be missed during detection or have incorrect indices.
 *
 * INPUT PARAMETERS
 * means: vector of pairs containing the mean and sigma values of the gaussian fit on the peaks
 * indices: vector of integers containing the index related to each peak (peak number)
 * fit: linear fit results including slope, error and reduced chi-squared
 *
 * WORKFLOW
 * A copy of the input indices "adjustedIndices" and a flag "adjusted" are created.
 * Then iterating through peaks, "rearrangement_tries" times:
 * - check if the left peak distance is significantly greater than the right peak distance;
 * - if true, increment indices to the right to adjust for the possible missed peak;
 * - check if the right peak distance is significantly greater that the left peak distance;
 * - if true, increment indices to the right to adjust for the possible missed peak;
 * - if any adjustment is made, recompute the linear fit with the new indices and if the new fit has lower reduced chi-squared update the fit results
 *
 * RETURNS
 * A vector of adjusted indices for the detected peaks.
 *
 * */
std::vector<int> rearrangeIndices(std::vector<std::pair<double, double>> means, std::vector<int> indices, lFit fit) {
  
	const double INCREMENT_FACTOR = 0.33;

	std::vector<int> adjustedIndices = indices;
	bool adjusted = false;

	for(int peak = 1; peak < means.size(); peak++) {
		double leftPeakDistance = means[peak].first - means[peak - 1].first;
                double rightPeakDistance = means[peak + 1].first - means[peak].first;
		for(int tries = 1; tries <= rearrangement_tries; tries++) {
			if (leftPeakDistance > (tries + INCREMENT_FACTOR) * rightPeakDistance){
				for (int i = peak; i < means.size(); i++){
					adjustedIndices[i]++;
				}
				adjusted = true;
			}

			if (rightPeakDistance > (tries + INCREMENT_FACTOR) * leftPeakDistance){
                                for (int i = peak + 1; i < means.size(); i++){
                                        adjustedIndices[i]++;
                                }
				adjusted = true;
			}

			if(adjusted == true) {
				lFit fitWithAdjustment = linearFit(means, adjustedIndices);
        			if (fitWithAdjustment.red_chi2 < fit.red_chi2) {
          				fit = fitWithAdjustment;
          				indices = adjustedIndices;
				}
				break;
			}
		}
	}
	
	return indices;
}

/**
 * Saves a canvas displaying the linear fit results of gain analysis.
 *
 * INPUT PARAMETERS
 * histo: Histogram object to be drawn on the canvas.
 * means: Vector containing mean and sigma values for peaks.
 * indices: Vector containing indices corresponding to each peak.
 * png_name: String specifying the name of the output PNG file.
 *
 * WORKFLOW
 * The function configures and draws histograms and fitted graphs on a canvas, organizing data presentation
 * in a manner that highlights the results of the linear fit. The canvas is then saved as a PNG file.
 */
void saveCanvas(TH1I* histo, std::vector<std::pair<double, double>> means, std::vector<int> indices, std::string png_name) {
	TCanvas* canvas = new TCanvas("canvas", "Channel's gain fit");
	canvas->Divide(1, 2);
	canvas->cd(1);
	histo->Draw("");

	// Prepare vectors for error band
    	canvas->cd(2);
	std::vector<double> x, y, ey;
	for (int i = 0; i < means.size(); ++i) {
    		x.push_back(indices[i]);
    		y.push_back(means[i].first);
    		ey.push_back(means[i].second);
	}

	// Create the error graph
	TGraphErrors grError(x.size(), &x[0], &y[0], nullptr, &ey[0]);

	// Customize the error graph to look like a band
	grError.SetFillColorAlpha(kBlue, 0.1);
	grError.SetFillStyle(3002); 
    	grError.SetLineColor(kRed);
    	grError.SetLineWidth(2);

    	// Draw the error graph first to make it a background
    	grError.Draw("A3"); 
    	canvas->Update();

    	// Then, overlay it with your line graphs
    	auto graph = fillGraph(means, indices);
    	graph.Draw("LP same");

    	canvas->Update();
  	TF1* lfit = new TF1("lfit", "pol1");
    	lfit->SetLineColor(2);
    	graph.Fit(lfit, "Q");
    	graph.Draw("LP same");
	canvas->Update();  
  	canvas->SaveAs(png_name.c_str());
  	delete canvas;
}

/**
 * Generates canvas displaying distributions of slopes and chi-squared values for quality assessment.
 *
 * INPUT PARAMETERS
 * slopes: Vector containing slope values derived from linear fits.
 * chi2s: Vector containing chi-squared values corresponding to each fit.
 *
 * WORKFLOW
 * The function creates histograms for both slope distributions and reduced chi-squared values.
 * These histograms are visualized on a split canvas to facilitate comparative analysis of fitting quality.
 */
void slopeChi2Canvas(std::vector<double> slopes, std::vector<double> chi2s, TString inputDir) {
	int nbins = 500;
	int max_range_slopes = 250;
	int max_range_chi2s = 1;
	std::string png_name;

	TH1D* slopes_histo = new TH1D("Slopes", "Slopes", nbins, 0, max_range_slopes);
  	TH1D* red_chi2s_histo = new TH1D("Reduced Chi2", "Reduced Chi2", nbins, 0, max_range_chi2s);

  	for (int i = 0; i < slopes.size(); i++) {
    		if (!isnan(slopes[i])) {
      			slopes_histo->Fill(slopes[i]);
      			red_chi2s_histo->Fill(chi2s[i]);
    		} else {
      			slopes_histo->Fill(-1);
     			red_chi2s_histo->Fill(-1);
    		}
  	}

  	TCanvas *canvas = new TCanvas("canvas", "Gain and reduced chi squared distributions");
  	canvas->Divide(1, 2);
  	canvas->cd(1);
  	red_chi2s_histo->Draw();
  	canvas->cd(2);
  	slopes_histo->Draw();

	png_name = std::string(inputDir.Data()) + "/" + dataset_specific_name + "_SlopeChi2" + ".png";
}

/**
 * The function retrieves the mean and sigma information from the pedestal's histograms
 *
 * INPUT PARAMETERS
 * fInPed: Input file, containing the pedestal distributions of all channels
 * inputDir: String specifying the directory containing the fInPed file
 *
 * WORKFLOW
 * Initialize vectors to store pedestal values (pedestals), peaks values (means) and sorted peak positions (xpeaks_sorted), 
 * as well as other variables and strings for histogram and FEB names.
 * An iterator is created to loop through keys in the input ROOT file. Looping through keys, the pedestal histograms are 
 * analysed. For each histogram, the FEB number specified in the name is extracted. 
 * If the FEB number is within a certain range, analyze and fit the pedestal histogram:
 * - perform a rebinning and smoothing of the histogram;
 * - determine the pedestal peak position and add it to xpeaks_sorted;
 * - call the "getPeaksMean" function to obtain mean and sigma values for the pedestal peak;
 * the histogram will now have the superposed gaussian fit function.
 * If no mean value is found, increase the chi-squared requirement and retry the fit.
 * If the fit still doesn't work after multiple attempts, consider it as failed and set the mean and sigma values to 0.
 * Finally save the pedestal peak information and histogram as an image (if save_png is true).
 *
 * RETURNS 
 * The vector "pedestals", containing the mean and sigma values for the pedestal peaks in the histograms
 *
 * */
std::vector<std::pair<double, double>> pedestal(TFile *fInPed, TString inputDir) {
  
	std::vector<std::pair<double, double>> pedestals;
	std::vector<std::pair<double, double>> means;
	std::vector<double> xpeaks_sorted;
	double pedestalPeak;
	int febNr, chi2_red;
	bool fitted;
	std::string pedHisto_name, png_name;
	TString tstring_pedHisto_name, tstring_febName;

	TIter next(fInPed->GetListOfKeys());
  	TKey* key;

  	while ((key = (TKey*)next())) {
		// Extract the histogram name, of the form "hadc_febNr_channerNr", and FEB number
   		pedHisto_name = key->GetName();
		tstring_pedHisto_name = pedHisto_name.c_str();
		TObjArray* tokens = tstring_pedHisto_name.Tokenize("_");
		TString tstring_febName = ((TObjString*)tokens->At(1))->GetString();
    		febNr = tstring_febName.Atoi();
		
		//Select the FEBs to analyze
		if(febNr > 107 && febNr != 113 && febNr < 139 && febNr != 126){
    			std::cout << "Analysing the pedestal: " << pedHisto_name << std::endl;
			fitted = false;
			chi2_red = chi2_red_pedestal;

			while(!fitted) {
				// Clone and analyze the pedestal histogram
				TH1I* histoPedestal = (TH1I*)fInPed->Get(tstring_pedHisto_name)->Clone();
				rebinAndSmooth(histoPedestal, pedestal_rebin, pedestal_smooth);

                		// Get the pedestal peak and add it to the sorted peaks
				pedestalPeak = histoPedestal->GetXaxis()->GetBinCenter(histoPedestal->GetMaximumBin());
				xpeaks_sorted.push_back(pedestalPeak);

                		// Call the "getPeaksMean" function to obtain mean and sigma values, and superposed gaussian function
				means = getPeaksMean(histoPedestal, xpeaks_sorted, chi2_red, "ped");

      				if (means.size() == 0) {
					// If no mean value is found, increase the chi-squared requirement
        				chi2_red += chi2_fit_increase;
					std::cout << "Pedestal fit didn't work for " << pedHisto_name << ": no mean value found. Trying again ..." << std::endl;
        				if (chi2_red > chi2_red_limit) {
						// If the chi-squared limit is reached, consider it as failed and set the mean and sigma values to 0
						std::cout << "Pedestal fit didn't work for " << pedHisto_name << ": no mean value found. Filling with (0,0)" << std::endl;
						means.push_back(std::make_pair(0, 0));
          					fitted = true;
        				}
      				} else {
					// If the fit is successful and a peak is found, exit the fitting loop
					fitted = true;
      				}
			
				// Save the histogram as an image (if save_png is true), in the range [0, range_pedestal_histo]
      				png_name = std::string(inputDir.Data()) + "/" + dataset_specific_name + "_" + pedHisto_name + ".png";
      				if (fitted == true) {
        				pedestals.push_back(std::make_pair(means[0].first, means[0].second));
        				if (save_png == true) {
          					TCanvas* canvas = new TCanvas("canvas", "canvas");
          					histoPedestal->GetXaxis()->SetRangeUser(0, range_pedestal_histo);
          					histoPedestal->Draw();
          					canvas->SaveAs(png_name.c_str());
          					delete canvas;
        				}	
      				}    				
				delete histoPedestal;
			}
		}
		delete tokens;
 	}
	return pedestals;
}

/**
 * Retrieves and analyzes gain information across all channels from both signal and pedestal histograms.
 *
 * INPUT PARAMETERS
 * fInSig: Input file containing histograms for signal analysis.
 * fInPed: Input file containing histograms for pedestal analysis.
 * inputDir: Directory containing the input file.
 *
 * WORKFLOW
 * The function first analyzes pedestal histograms to extract their mean and sigma values, and then processes
 * signal histograms to determine gain characteristics of each channel. It employs various fitting procedures
 * to model the peak distributions and extract peaks mean and standard deviation values.
 * The results are used to calculate the gain across multiple histograms (channels).
 *
 * RETURNS
 * A vector of 'gain' structures, each containing the computed gain information for a channel, including
 * calculated slopes, pedestal values, and errors.
 */
std::vector<gain> hits(TFile *fInSig, TFile *fInPed, TString inputDir){
		
	//Compute and store the mean and sigma values, from gaussian fits of the channels pedestal histograms 
	std::vector<std::pair<double, double>> pedestals = pedestal(fInPed, inputDir);
		
	std::vector<std::pair<double, double>> means;
  	std::vector<double> chi2s, slopes, slopes_error, xpeaks_sorted;
  	std::vector<gain> gains;
	std::string sigHisto_name, png_name;
	TString tstring_sigHisto_name, tstring_febNr;
	int febNr, chi2_red, histogramNr;
	bool fitted;

	TIter next(fInSig->GetListOfKeys());
	TKey* key;

	histogramNr = 0;
  	while ((key = (TKey*)next())) {
		gain gain_entry;
		
		// Extract the histogram name, of the form "hadc_febNr_channerNr", and FEB number
		sigHisto_name = key->GetName();
		tstring_sigHisto_name = sigHisto_name.c_str();
		TObjArray* tokens = tstring_sigHisto_name.Tokenize("_");
		tstring_febNr = ((TObjString*)tokens->At(1))->GetString();
    		febNr = tstring_febNr.Atoi();
		
		//Select the FEBs to analyze, as some FEBs do not work properly
		if(febNr > 107 && febNr != 113 && febNr < 139 && febNr != 126){
    			std::cout << "Analysing the signal: " << sigHisto_name << std::endl;
			std::stringstream namestream(sigHisto_name);
			std::string segment;
			std::vector<std::string> segment_list;
			while(std::getline(namestream, segment, '_')){
				segment_list.push_back(segment);		
			}
			gain_entry.feb = segment_list[1];
			gain_entry.channel = segment_list[2];

			fitted = false;
			chi2_red = chi2_red_signal;
			while(!fitted){
				// Clone and analyze the signal histogram
      				TH1I* histoSignal = (TH1I*)fInSig->Get(sigHisto_name.c_str())->Clone();
      				rebinAndSmooth(histoSignal, signal_rebin, signal_smooth);

				// Extract the ordered peak positions from the signal histogram, keeping only the first 'peaks_of_interest' 
				xpeaks_sorted = getPeaks(histoSignal, max_peaks);
				if (xpeaks_sorted.size() > peaks_of_interest){
					xpeaks_sorted.erase(xpeaks_sorted.begin() + peaks_of_interest, xpeaks_sorted.end());
				}

				// Extract the mean and sigma values from recursive gaussian fits of the p.e. peaks
      				means = getPeaksMean(histoSignal, xpeaks_sorted, chi2_red, "sig"); 

				// Initialize the vector containing the indices of the peaks (peak's number)
      				std::vector<int> index(means.size());
      				std::iota(index.begin(), index.end(), 0);
      				std::vector<int> index_error(means.size(), 0);

				// An additional indices vector, without the first element, is created
      				std::pair<double, double> first_point = means[0];
      				auto means_no_first = means;
      				means_no_first.erase(means_no_first.begin());
      				auto tmp_index = index;
      				tmp_index.erase(tmp_index.begin());

				// Compute a linear fit of the p.e. peaks ADC count mean value vs the corresponding p.e. peak number, with and without the first p.e. peak index
				lFit fit = linearFit(means, index);
      				lFit fit_no_first = linearFit(means_no_first, tmp_index);

				// Compare the two, keeping the linear fit with lower reduced chi squared value
      				bool better_no_first = false;
      				if (fit_no_first.red_chi2 < fit.red_chi2) {
					std::cout << "For the histogram: " << sigHisto_name << " the fit_no_first is better." << std::endl;
     					better_no_first = true;
     					means = means_no_first;
     					index = tmp_index;
     					fit = fit_no_first;
      				}

      				png_name = std::string(inputDir.Data()) + "/" + dataset_specific_name + "_" + sigHisto_name + ".png";

				// If no p.e. peak is found, fill the chi2, slopes and slopes_error with -1
      				if (means.size() == 0) {
      					chi2s.push_back(-1);
      					slopes.push_back(-1);
      					slopes_error.push_back(-1);
      					break;
     				}

				// Search for skipped peaks (not found by TSpectrum) and rearrange the indices, then repeat the linear gain fit
				index = rearrangeIndices(means, index, fit);
      				fit = linearFit(means, index);

				// Add the informations of the pedestal to the gain_entry of the current channel, and in the means vector
				gain_entry.pedestal = pedestals[histogramNr].first;
				gain_entry.pedestal_error = pedestals[histogramNr].second; 
      				means.insert(means.begin(), std::make_pair(gain_entry.pedestal, gain_entry.pedestal_error));
				
				for(int peak_nr = 0; peak_nr < index.size(); peak_nr++){ 
					index[peak_nr]++;
				}
				index.insert(index.begin(), 0);

				// Try a final fit comparison, as sometimes the first p.e. peaks are not found. Keep the best linear fit, based on the reduced chi squared value 
      				tmp_index = index;
      				fit = linearFit(means, tmp_index);
      				int tries = 0;
      				while (tries < 3) {
        				for (int peak_nr = 1; peak_nr < means.size(); peak_nr++){ 
						tmp_index[peak_nr]++;
					}
        				lFit fit_rearranged = linearFit(means, tmp_index);
        				
					if (fit_rearranged.red_chi2 < fit.red_chi2) {
          					fit = fit_rearranged;
          					index = tmp_index;
        				}
        				tries++;
      				}

				// Some checks of the fit quality
			   	if (fit.error / fit.slope > 0.3) {
        				std::cout << "The following histogram: " << sigHisto_name << " has a large error" << std::endl;
					std::cout << "Reducing the red_chi2 and fitting again." << std::endl;
					chi2_red--;
        				if (chi2_red < 1) {
          					fitted = true;
        				}
      				} else {
        				fitted = true;
      				}

				// When the gain fit procedure ends, store the obtained red_chi2 and slope (gain) values
				if(fitted == true){
					chi2s.push_back(fit.red_chi2);
					
					if (fit.red_chi2 > 0.3) { 
						std::cout << "The following histogram: " << sigHisto_name << " has red_chi2 > 0.3" << std::endl;
					}
        				
					slopes.push_back(fit.slope);
        				slopes_error.push_back(fit.error);

        				gain_entry.slope = fit.slope;
        				gain_entry.slope_error = fit.error;
        				gain_entry.index = index;
        
        				gains.push_back(gain_entry);
				}		

				// Draws and saves as a png file the graph of the "p.e. peak number vs p.e. peak ADC count mean value", with overlaid gain fit 
         			if (save_png == true) {
					 saveCanvas(histoSignal, means, index, png_name);
				}	
			
			}
		
			histogramNr++;
		}
	
		delete tokens;
	}

	// Draws the channels' distribution of slopes and reduced chi2 values
	slopeChi2Canvas(slopes, chi2s, inputDir);
	
	return gains;
}

/**
 * Main function, defining the entire histogram analysis workflow for both signal and pedestal data.
 *
 * INPUT PARAMETERS
 * fIn: String pointing to the input file location (path and file name).
 *
 * WORKFLOW
 * The function starts by splitting the input histograms into signal and pedestal groups and storing them
 * in separate files. It then processes these histograms to extract gain information using the 'hits' function.
 * A series of plots are then represented, with the slope (gain) and pedestal distributions information.
 * A .csv file with the main analysis information is generated, for further analysis.
 */
void histosAnalysis(const char* fIn){
	//Divide the name of the file and path, passed as argument
	TString inputFilePath(fIn);
	TString inputDir = gSystem->DirName(inputFilePath);
 	TString inputFileName = gSystem->BaseName(inputFilePath);

	//Divide the histograms, generated in the decoding stage, in two files. One for the signal and one for the pedestal.
  	TString signalFileName = inputDir + "/signal_histograms_" + inputFileName;
  	TString pedestalFileName = inputDir + "/pedestal_histograms_" + inputFileName;
	TFile *signalFile = new TFile(signalFileName, "RECREATE");
  	TFile *pedestalFile = new TFile(pedestalFileName, "RECREATE");
  	separateHistograms(fIn, signalFile, pedestalFile);

	// Start the analysis and save the gain results in a vector of structs
	std::vector<gain> hits_vec = hits(signalFile, pedestalFile, inputDir);

	// Draw a series of plots, with the results of the analysis on the gain and its error, pedestal mean and sigma values 
  	TCanvas* c = new TCanvas("c", "c");
  	c->Divide(2, 2);
  	c->cd(1);

  	TH1D* slope_histo = new TH1D("slope", "slope", 400, 0, 400);
  	TH1D* slope_error_histo = new TH1D("slope_error", "slope_error", 400, 0, 80);
  	TH2D* slope_vs_slope_error_histo = new TH2D("Slope vs Slope_error", "Slope vs Slope_error", 400, 0, 400, 400, 0, 80);

  	TH1D* pedestal_histo = new TH1D("pedestal", "pedestal", 400, 0, 400);
  	TH1D* pedestal_error_histo = new TH1D("pedestal_e", "pedestal_e", 400, 0, 80);
  	TH2D* pedestal_vs_pedestal_error_histo = new TH2D("pedestal vs pedestal_error", "pedestal vs pedestal_error", 400, 0, 400, 400, 0, 80);

  	TH2D* slope_vs_pedestal_histo = new TH2D("Slope vs pedestal", "Slope vs pedestal", 400, 0, 400, 400, 0, 400);

  	for (int i = 0; i < hits_vec.size(); i++) {
    		slope_histo->Fill(hits_vec[i].slope);
    		slope_error_histo->Fill(hits_vec[i].slope_error);
    		slope_vs_slope_error_histo->Fill(hits_vec[i].slope, hits_vec[i].slope_error);
  	}
	
	TCanvas* canvas1 = new TCanvas("canvas1", "canvas1");
  	canvas1->Divide(2, 2);
  	canvas1->cd(4);
  	slope_histo->Draw();
  	canvas1->cd(1);
  	slope_error_histo->Draw();
  	canvas1->cd(2);
  	slope_vs_slope_error_histo->Draw("colz");

  	for (int i = 0; i < hits_vec.size(); i++) {
    		pedestal_histo->Fill(hits_vec[i].pedestal);
    		pedestal_error_histo->Fill(hits_vec[i].pedestal_error);
    		pedestal_vs_pedestal_error_histo->Fill(hits_vec[i].pedestal, hits_vec[i].pedestal_error);
  	}

  	TCanvas* c2 = new TCanvas("c2", "c2");
  	c2->Divide(2, 2);
  	c2->cd(4);
  	pedestal_histo->Draw();
  	c2->cd(1);
  	pedestal_error_histo->Draw();
  	c2->cd(2);
  	pedestal_vs_pedestal_error_histo->Draw("colz");

  	for (int i = 0; i < hits_vec.size(); i++) slope_vs_pedestal_histo->Fill(hits_vec[i].slope, hits_vec[i].pedestal);
  	TCanvas* c3 = new TCanvas("c3", "c3");
  	c3->cd(1);
  	slope_vs_pedestal_histo->Draw("colz");

	
	// Save the calibration information in a .csv file
  	std::ofstream outfile;
  	outfile.open("topCRT_calibration.csv");
  	outfile << "channel,mac5,localchannel,gain,pedestal,gainflags\n";
  	outfile << "bigint,integer,integer,real,real,integer\n";

  	int globalIndex = 107 * 32;
  	for (int i = 0; i < hits_vec.size(); i++) {
    		outfile << globalIndex + i << "," << hits_vec[i].feb << ","
            		<< hits_vec[i].channel << "," << hits_vec[i].slope << ","
            		<< hits_vec[i].pedestal << "," << 0 << "\n";
  	}

	// Final operations
  	signalFile->Close();
  	pedestalFile->Close();
}
