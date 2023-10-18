//THIS MACRO DOES AN ANALYSIS ON THE HISTOGRAMS CONTAINED AT THE FOLLOWING PATH /pnfs/icarus/scratch/users/vpia/cal/run_9989/v09_77_00/ana/test_poms_cal/out/calibrationHistos/ THAT CONTAINS THE MERGED HISTOGRAMS OBTAINED AFTER STAGE0 OF OUR CAMPAIGN

//You will then need to add to the code the custom hadd and the access to the list of broken channels, and so the gain dummy values for those channels
//Change the png name so that takes into account the campaign from which is produced
//Brazilian plot
//Segmentation error possibly from loop on max peaks (make a check for this kind of error). Just check it for the seen histogram.

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
int chi2_red_pedestal = 10;
int chi2_red_signal = 2;
int pedestal_rebin = 5;
int pedestal_smooth = 0;
int signal_rebin = 10;
int signal_smooth = 2;
int adc_last_pe = 900;
int range_pedestal_histo = 500;
int max_peaks = 15;
int peaks_of_interest = 5;
double sigmaXPeakSearch = 1;
double peaks_treshold = 0.10;
std::string dataset_specific_name = "afterCut";

//Represents gain information, including pedestal values, slope, peaks, and associated errors.
struct gain {
  double pedestal;
  double pedestal_e;

  double slope;
  double slope_e;

  std::vector<double> peaks;
  std::vector<double> peaks_e;

  std::vector<int> index;

  std::string feb;
  std::string channel;
};

//Represents the result of a linear fit, containing the slope, error, and reduced chi-squared value.
struct lFit {
  double slope;
  double err;
  double red_chi2;
};

//Creates 2 separate files for the pedestal and signal histograms
void separateHistograms(const char* fIn, TFile *fInSig, TFile *fInPed) {
    TFile *inputFile = TFile::Open(fIn, "READ");

    TDirectory *mainDir = (TDirectory*)inputFile->Get("CRTCalibrationAnalysis");
    TList *list = mainDir->GetListOfKeys();
		bool hasSubDirs = false;

    TIter iter(list);
    TKey *key;
    while ((key = (TKey*)iter())) {
        if (strcmp(key->GetClassName(), "TDirectoryFile") == 0) {
            hasSubDirs = true;
            break;
        }
    }

    for (int feb = 1; feb <= 231; ++feb) {
        for (int ch = 0; ch <= 31; ++ch) {
            TString signalHistName = Form("hadc_%d_%d_signal", feb, ch);
            TString pedestalHistName = Form("hadc_%d_%d_pedestal", feb, ch);
            TH1 *signalHist = nullptr;
            TH1 *pedestalHist = nullptr;

            if (hasSubDirs) {

                TIter dirIter(list);
                while ((key = (TKey*)dirIter())) {
                    if (strcmp(key->GetClassName(), "TDirectoryFile") != 0) {
                        continue;
                    }
                    TDirectory *subDir = (TDirectory*)mainDir->Get(key->GetName());
                    signalHist = (TH1*)subDir->Get(signalHistName);
                    pedestalHist = (TH1*)subDir->Get(pedestalHistName);

                    if (signalHist || pedestalHist) {
                        break;
                    }
                }
            } else {
                signalHist = (TH1*)mainDir->Get(signalHistName);
                pedestalHist = (TH1*)mainDir->Get(pedestalHistName);
            }

            if (signalHist) {
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

//Smoothing and rebinning of the histograms, eliminating also the under/overflow
void rebinAndSmooth(TH1I* histo, int rebin, int smooth) {
	for ( int nbin = histo->GetNbinsX(); nbin > 0; nbin--){
		if( histo->GetBinContent(nbin) != 0) {
			histo->SetBinContent(nbin, 0);
			break;
		}
	}

	histo->SetBinContent(1,0);
	histo->Rebin(rebin);
	histo->Smooth(smooth);

	return;
}

//Returns the mean and sigma value for the peaks of the passed histogram. The info are used for the fit
std::vector<std::pair<double, double>> getPeaksMean(TH1I* histo ,  std::vector<double> xpeaks_sorted, int requested_red_chi2, std::string histo_type) {
	std::vector<std::pair<double, double>> means;
	for(int npeak = 0; npeak < xpeaks_sorted.size(); npeak++){
    TAxis* xaxis = histo->GetXaxis();
    int low_binx;
    int up_binx;

    double left_distance;
    if (npeak == 0) {
      left_distance = xpeaks_sorted[npeak];
    } else {
      left_distance = xpeaks_sorted[npeak] - xpeaks_sorted[npeak - 1];
    }

    double right_distance;
    if (npeak == xpeaks_sorted.size() - 1) {
      right_distance = adc_last_pe - xpeaks_sorted[npeak];
    } else {
      right_distance = xpeaks_sorted[npeak + 1] - xpeaks_sorted[npeak];
    }

		int den_fit_range;
		if(histo_type == "ped") den_fit_range = 1;
		if(histo_type == "sig") den_fit_range = 2;

		double half_mean_distance;
		if(left_distance < right_distance){
      half_mean_distance = left_distance / den_fit_range;
    } else {
      half_mean_distance = right_distance / den_fit_range;
    }

		low_binx = xaxis->FindBin(xpeaks_sorted[npeak] - half_mean_distance);
    up_binx  = xaxis->FindBin(xpeaks_sorted[npeak] + half_mean_distance);

		bool fitted = false;
    double red_chi2;
    int tries = 0;
    while (!fitted) {
      double low_x = xaxis->GetBinCenter(low_binx);
      double up_x  = xaxis->GetBinCenter(up_binx);
      TF1* gaussian_fit = new TF1("gfit", "gaus", low_x, up_x);
      gaussian_fit->SetLineColor(2);
      histo->Fit(gaussian_fit, "QR+");

      double mean = gaussian_fit->GetParameter(1);
      double sigma = gaussian_fit->GetParameter(2);
      red_chi2 = gaussian_fit->GetChisquare() / gaussian_fit->GetNDF();
	  	if (red_chi2 < requested_red_chi2 &&
          mean > xpeaks_sorted[npeak] - half_mean_distance &&
          mean < xpeaks_sorted[npeak] + half_mean_distance && 
          mean + sigma < adc_last_pe) {
        means.push_back(std::make_pair(mean, sigma));
        fitted = true;
      }

			low_binx++;
			up_binx--;

			if (up_binx - low_binx < 2) fitted = true;

			if (tries == 10) fitted = true;
			tries++;

			delete gaussian_fit;
		} 
	}
	return means;
}

//Extracts peak positions from a histogram. Returns: x-values of the detected peaks. TSpectrum maybe obsolete.
std::vector<double> getPeaks(TH1I* histo, int npeaks) {
  TSpectrum* spectra = new TSpectrum(npeaks);
  histo->GetXaxis()->SetRangeUser(0, adc_last_pe);

  int nfound = spectra->Search(histo, sigmaXPeakSearch, "", peaks_treshold);
  double* xpeaks = spectra->GetPositionX();

  std::vector<double> xpeaks_sorted;
  for (int peak = 0; peak < nfound; peak++) {
    xpeaks_sorted.push_back(xpeaks[peak]);
  }
  std::sort(xpeaks_sorted.begin(), xpeaks_sorted.end());

  delete spectra;
  return xpeaks_sorted;
}

//Fills a Tgraph with peaks data and error bars
TGraphErrors fillGraph(std::vector<std::pair<double, double>> means, std::vector<int> indices) {
  TGraphErrors graph(means.size());
  for (int peak_mean_nr = 0; peak_mean_nr < means.size(); peak_mean_nr++) {
    graph.SetPoint(peak_mean_nr, indices[peak_mean_nr], means[peak_mean_nr].first);
    graph.SetPointError(peak_mean_nr, 0, means[peak_mean_nr].second); //Second argument: error on x
  }

  return graph;
}

//Performs a linear fit on given data points.
lFit linearFit(std::vector<std::pair<double, double>> means, std::vector<int> indices) {
  TF1* linear_fit = new TF1("lfit", "pol1");
  auto graph = fillGraph(means, indices);
  linear_fit->SetLineColor(2);
  graph.Fit(linear_fit, "Q");

  lFit my_fit;
  my_fit.slope = linear_fit->GetParameter(1);
  my_fit.err = linear_fit->GetParError(1);
  my_fit.red_chi2 = linear_fit->GetChisquare() / linear_fit->GetNDF();

  delete linear_fit;
  return my_fit;
}

//Adjusts the indices of detected peaks in an attempt to improve the results of a linear fit (in some fits a peak is missed and the indices must be shifted)
std::vector<int> rearrangeIndex(std::vector<std::pair<double, double>> means,
                                std::vector<int> indices, lFit fit) {
  std::vector<int> tmp_indices = indices;
  int tries = 1;
  while (tries < 3) {
    if (ceil(means[1].first - means[0].first) > (tries + 0.33) * (means[2].first - means[1].first)) {
      for (int peak_nr = 1; peak_nr < means.size(); peak_nr++) tmp_indices[peak_nr]++;
    }
    tries++;
  }
  indices = tmp_indices;

  int peak = 1;
  while (peak < means.size()) {
    tries = 1;
    bool skip = false;
    while (tries < 3) {
      if (ceil(means[peak + 1].first - means[peak].first) >
          (tries + 0.33) * (means[peak].first - means[peak - 1].first)) {
        for (int peak_nr = peak + 1; peak_nr < means.size(); peak_nr++) 
          tmp_indices[peak_nr]++;
      } else {
        skip = true;
        break;
      }

      if (skip == false) {
        lFit fit_rearranged = linearFit(means, tmp_indices);
        if (fit_rearranged.red_chi2 < fit.red_chi2) {
          fit = fit_rearranged;
          indices = tmp_indices;
        }
        tries++;
      }
    }
    peak++;
  }
  return indices;
}

//Saves the canvas for the linear fit
void saveCanvas(TH1I* histo, std::vector<std::pair<double, double>> means,
                std::vector<int> indices, std::string png_name) {
  TCanvas* canvas = new TCanvas("canvas", "canvas");
  canvas->Divide(1, 2);
  canvas->cd(1);
  histo->Draw("");

	// Prepare vectors for error band
    canvas->cd(2);
		std::vector<double> x, y, ey;
		for (int i = 0; i < means.size(); ++i) {
    	x.push_back(indices[i]);
    	y.push_back(means[i].first);
    	ey.push_back(means[i].second);  // Replace with your actual error values
		}

		// Create the error graph
		TGraphErrors grError(x.size(), &x[0], &y[0], nullptr, &ey[0]);

		// Customize the error graph to look like a band
		grError.SetFillColorAlpha(kBlue, 0.1);
		grError.SetFillStyle(3002); 
    grError.SetLineColor(kRed);
    grError.SetLineWidth(2);

    // Draw the error graph first to make it a background
    grError.Draw("A3");  // "A3" means draw axes and the fill area
    canvas->Update();  // Update the canvas

    // Then, overlay it with your line graphs
    auto graph = fillGraph(means, indices);
    graph.Draw("LP same");  // "LP" means draw line and markers, "same" means overlay

    canvas->Update();  // Update the canvas again
  	TF1* lfit = new TF1("lfit", "pol1");
    lfit->SetLineColor(2);
    graph.Fit(lfit, "Q");
    graph.Draw("LP same");
		canvas->Update();  

  	canvas->SaveAs(png_name.c_str());
  	delete canvas;
}

//Creates two canvas to show the distributions of slopes and reduced chi2 values
void slopeChi2Canvas(std::vector<double> slopes, std::vector<double> chi2s) {
	int nbins = 500;
	int max_range_slopes = 250;
	int max_range_chi2s = 1;

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

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->Divide(1, 2);
  canvas->cd(1);
  red_chi2s_histo->Draw();
  canvas->cd(2);
  slopes_histo->Draw();
}

//Restrieves gain information from the histograms with pedestal of the pedestal
std::vector<std::pair<double, double>> pedestal(TFile *fInPed, TString inputDir) {
  
	std::vector<std::pair<double, double>> pedestals;

	TIter next(fInPed->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)next())) {
    std::string pedHisto_name = key->GetName();
		TString tstr_pedHisto_name = pedHisto_name.c_str();
		TObjArray* tokens = tstr_pedHisto_name.Tokenize("_");
		TString febStr = ((TObjString*)tokens->At(1))->GetString();
    int febNr = febStr.Atoi();

		if(febNr > 107 && febNr != 113 && febNr < 139 && febNr != 126){
    std::cout << "Analysing the pedestal: " << pedHisto_name << std::endl;

		bool fitted = false;
		int chi2_red = chi2_red_pedestal;
		while(!fitted) {
			TH1I* hPed = (TH1I*)fInPed->Get(pedHisto_name.c_str())->Clone();
			rebinAndSmooth(hPed, pedestal_rebin, pedestal_smooth);

			double pedestalPeak = hPed->GetXaxis()->GetBinCenter(hPed->GetMaximumBin());

			std::vector<double> xpeaks_sorted;
			xpeaks_sorted.push_back(pedestalPeak);

			std::vector<std::pair<double, double>> means = getPeaksMean(hPed, xpeaks_sorted, chi2_red, "ped");

      if (means.size() == 0) {
        chi2_red += 10;
        if (chi2_red > 500) {
          means.push_back(std::make_pair(0, 0));
          fitted = true;
        }
      } else {
        fitted = true;
      }
			
      std::string png_name = std::string(inputDir.Data()) + "/" + dataset_specific_name + "_" + pedHisto_name + ".png";
      if (fitted == true) {
        pedestals.push_back(std::make_pair(means[0].first, means[0].second));

        if (save_png == true) {
          TCanvas* canvas = new TCanvas("canvas", "canvas");
          hPed->GetXaxis()->SetRangeUser(0, range_pedestal_histo);
          hPed->Draw();
          canvas->SaveAs(png_name.c_str());
          delete canvas;
        }
      }
      delete hPed;
		}
	}
	delete tokens;
 }
	return pedestals;

}

//Retrieves gain information for all hits
std::vector<gain> hits(TFile *fInSig, TFile *fInPed, TString inputDir){
  
	std::vector<std::pair<double, double>> pedestals = pedestal(fInPed, inputDir);
  std::vector<double> chi2s;
  std::vector<double> slopes;
  std::vector<double> slopes_e;
  std::vector<gain> gains;
	std::map<int, std::vector<double>> indexToSigmas;
/*	// Declare two histograms, one for index and one for sigma
	TH1D* indexHistogram = new TH1D("indexHistogram", "Index Distribution;Index;Frequency", 6, 0, 6); 
	TH1D* sigmaHistogram = new TH1D("sigmaHistogram", "Sigma Distribution;Sigma;Frequency", 100, 0, 500);
	TH2D* indexVsSigma = new TH2D("indexVsSigma", "Index vs Sigma;Index;Sigma", 7, 0, 7, 100, 0, 300);
*/

	TIter next(fInSig->GetListOfKeys());
	TKey* key;

	int hist_number = 0;
  while ((key = (TKey*)next())) {
		gain gain_entry;
		
		std::string sigHisto_name = key->GetName();
		TString tstr_sigHisto_name = sigHisto_name.c_str();
		TObjArray* tokens = tstr_sigHisto_name.Tokenize("_");
		TString febStr = ((TObjString*)tokens->At(1))->GetString();
    int febNr = febStr.Atoi();

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

		bool fitted = false;
		int chi2_red = chi2_red_signal;
		while(!fitted){
      TH1I* hSig = (TH1I*)fInSig->Get(sigHisto_name.c_str())->Clone();
      rebinAndSmooth(hSig, signal_rebin, signal_smooth);

			std::vector<double> xpeaks_sorted = getPeaks(hSig, max_peaks);
			if (xpeaks_sorted.size() > peaks_of_interest){
				xpeaks_sorted.erase(xpeaks_sorted.begin() + peaks_of_interest, xpeaks_sorted.end());
			}

      std::vector<std::pair<double, double>> means = getPeaksMean(hSig, xpeaks_sorted, chi2_red, "hit"); 

      std::vector<int> index(means.size());
      std::iota(index.begin(), index.end(), 0);
      std::vector<int> index_e(means.size(), 0);

      std::pair<double, double> first_point = means[0];
      auto means_no_first = means;
      means_no_first.erase(means_no_first.begin());
      auto tmp_index = index;
      tmp_index.erase(tmp_index.begin());

			lFit fit = linearFit(means, index);
      lFit fit_no_first = linearFit(means_no_first, tmp_index);

      bool better_no_first = false;
      if (fit_no_first.red_chi2 < fit.red_chi2) {
     		better_no_first = true;
     		means = means_no_first;
     		index = tmp_index;
     		fit = fit_no_first;
      }

      std::string png_name = std::string(inputDir.Data()) + "/" + dataset_specific_name + "_" + sigHisto_name + ".png";

      if (means.size() == 0) {
      	chi2s.push_back(-1);
      	slopes.push_back(-1);
      	slopes_e.push_back(-1);
      	break;
     	}
			
			index = rearrangeIndex(means, index, fit);
      fit = linearFit(means, index);

			gain_entry.pedestal = pedestals[hist_number].first;
			gain_entry.pedestal_e = pedestals[hist_number].second; 
      means.insert(means.begin(), std::make_pair(gain_entry.pedestal, gain_entry.pedestal_e));
			for(int peak_nr = 0; peak_nr < index.size(); peak_nr++) index[peak_nr];
			index.insert(index.begin(), 0);

      tmp_index = index;
      fit = linearFit(means, tmp_index);
      int tries = 0;
      while (tries < 3) {
        for (int peak_nr = 1; peak_nr < means.size(); peak_nr++) tmp_index[peak_nr]++;
        lFit fit_rearranged = linearFit(means, tmp_index);
        if (fit_rearranged.red_chi2 < fit.red_chi2) {
          fit = fit_rearranged;
          index = tmp_index;
        }
        tries++;
      }

			std::cout << "For " << sigHisto_name << " The fit error/slope is: " << std::to_string(fit.err/fit.slope) << " and the reduced chi2: " << std::to_string(fit.red_chi2) << std::endl;



			for (int peak_mean_nr = 0; peak_mean_nr < means.size(); peak_mean_nr++) {
			  indexToSigmas[index[peak_mean_nr]].push_back(means[peak_mean_nr].second);
			}

/*
    //GRAFICO DEVIAZIONI STANDARD
  	for (int peak_mean_nr = 0; peak_mean_nr < means.size(); peak_mean_nr++) {
        indexHistogram->Fill(index[peak_mean_nr]);
   		  sigmaHistogram->Fill(means[peak_mean_nr].second);
				indexVsSigma->Fill(index[peak_mean_nr], means[peak_mean_nr].second);
    }
*/
			fitted = true;
/*    	if (fit.err / fit.slope > 0.3) {
        chi2_red--;
        if (chi2_red < 1) {
          fitted = true;
        }
      } else {
        fitted = true;
      }
*/
//		if(fitted == true){
			chi2s.push_back(fit.red_chi2);
			if (fit.red_chi2 > 0.3) std::cout << "Good one: " << sigHisto_name << std::endl;
			
        slopes.push_back(fit.slope);
        slopes_e.push_back(fit.err);

        gain_entry.slope = fit.slope;
        gain_entry.slope_e = fit.err;
        gain_entry.index = index;
        
        gains.push_back(gain_entry);
//				}		
      	 
         if (save_png == true) saveCanvas(hSig, means, index, png_name);	
		}
		hist_number++;
	}
	delete tokens;
	}
	slopeChi2Canvas(slopes, chi2s);

//-------------------------------------------------------------------

	std::map<int, double> indexToAvgSigma;
	for (auto& entry : indexToSigmas) {
  	double sum = std::accumulate(entry.second.begin(), entry.second.end(), 0.0);
 		double avg = sum / entry.second.size();
  	indexToAvgSigma[entry.first] = avg;
	}


	TGraph graph(indexToAvgSigma.size());
	int pointIndex = 0;
	for (auto& entry : indexToAvgSigma) {
		if(pointIndex > 6) break;
	  graph.SetPoint(pointIndex, entry.first, entry.second);
	  pointIndex++;
	}

	TF1* fit_linear = new TF1("fit_linear", "pol1");
	TF1* fit_quad = new TF1("fit_quad", "pol2");
	TF1* fit_cubic = new TF1("fit_cubic", "pol3");

	graph.Fit(fit_linear, "Q");
	double chi2_linear = fit_linear->GetChisquare() / fit_linear->GetNDF();

	graph.Fit(fit_quad, "Q");
	double chi2_quad = fit_quad->GetChisquare() / fit_quad->GetNDF();

	graph.Fit(fit_cubic, "Q");
	double chi2_cubic = fit_cubic->GetChisquare() / fit_cubic->GetNDF();

	TF1* best_fit;

	// Decide which fit is better based on reduced chi2 values
	if(chi2_linear <= chi2_quad && chi2_linear <= chi2_cubic) {
    best_fit = fit_linear;
    std::cout << "Best fit is linear with chi2: " << chi2_linear << std::endl;
	} else if(chi2_quad <= chi2_linear && chi2_quad <= chi2_cubic) {
    best_fit = fit_quad;
    std::cout << "Best fit is quadratic with chi2: " << chi2_quad << std::endl;
	} else {
    best_fit = fit_cubic;
    std::cout << "Best fit is cubic with chi2: " << chi2_cubic << std::endl;
	}

	// Draw the graph and the best fit
	TCanvas* canvas = new TCanvas("canvas_best_fit", "Best Fit");
	graph.Draw("AP");  // Draw graph with axis and points
	best_fit->SetLineColor(kRed); // Set the color of the best fit line to red
	best_fit->SetLineColor(kRed);  // Set the color of the best fit line to red
	best_fit->SetLineStyle(2);    // Set line style to dashed
	best_fit->SetLineWidth(1);    // Set line width to 1 (default is usually 3)
	best_fit->Draw("same");  // Draw the best fit on the same canvas
	canvas->SaveAs("best_fit.png");  // Save the canvas
	delete canvas;

	// Set graph properties to make points clearly visible
	graph.SetMarkerStyle(20);  // Set the marker style (20 is a full circle)
	graph.SetMarkerSize(1.5);  // Set the marker size to make it bigger
	graph.SetMarkerColor(kBlue);  // Set the marker color to blue

	// Draw the graph with points only
	TCanvas* canvas2 = new TCanvas("canvas_points_only", "Points Only");
	graph.Draw("AP");  // Draw graph with axis and points
	canvas2->SaveAs("points_only.png");  // Save the canvas
	delete canvas2;
/*
	TCanvas* c = new TCanvas("c", "c");
	c->Divide(1, 2); // 1 column, 2 rows

	c->cd(1);
	indexHistogram->Draw();
	c->cd(2); 
	sigmaHistogram->Draw();

	c->SaveAs("Index_and_Sigma_Distributions.png");
	
	TCanvas* c2 = new TCanvas("canvas_indexVsSigma", "Index vs Sigma");
	indexVsSigma->Draw("COLZ");
	c2->SaveAs("Index_vs_Sigma.png");
*/
//	delete indexHistogram;
//	delete sigmaHistogram;
//	delete c;
//  delete c2;
//	delete indexVsSigma;
//--------------------------------------------------------------------------------

	return gains;
}

//The arguments are: the file name and path
void histosAnalysis(const char* fIn){

  TString inputFilePath(fIn);
  TString inputDir = gSystem->DirName(inputFilePath);
  TString inputFileName = gSystem->BaseName(inputFilePath);
  TString signalFileName = inputDir + "/signal_histograms_" + inputFileName;
  TString pedestalFileName = inputDir + "/pedestal_histograms_" + inputFileName;

  TFile *signalFile = new TFile(signalFileName, "RECREATE");
  TFile *pedestalFile = new TFile(pedestalFileName, "RECREATE");

  separateHistograms(fIn, signalFile, pedestalFile);

	std::vector<gain> hits_vec = hits(signalFile, pedestalFile, inputDir);


	//Code to adjust
  TCanvas* c = new TCanvas("c", "c");
  c->Divide(2, 2);
  c->cd(1);

  TH1D* slope_histo = new TH1D("slope", "slope", 400, 0, 400);
  TH1D* slope_e_histo = new TH1D("slope_e", "slope_e", 400, 0, 80);
  TH2D* slope_vs_slope_e_histo =
      new TH2D("Slope vs Slope_e", "Slope vs Slope_e", 400, 0, 400, 400, 0, 80);

  TH1D* pedestal_histo = new TH1D("pedestal", "pedestal", 400, 0, 400);
  TH1D* pedestal_e_histo = new TH1D("pedestal_e", "pedestal_e", 400, 0, 80);
  TH2D* pedestal_vs_pedestal_e_histo =
      new TH2D("pedestal vs pedestal_e", "pedestal vs pedestal_e", 400, 0, 400,
               400, 0, 80);

  TH2D* slope_vs_pedestal_histo = new TH2D(
      "Slope vs pedestal", "Slope vs pedestal", 400, 0, 400, 400, 0, 400);

  for (int i = 0; i < hits_vec.size(); i++) {
    slope_histo->Fill(hits_vec[i].slope);
    slope_e_histo->Fill(hits_vec[i].slope_e);
    slope_vs_slope_e_histo->Fill(hits_vec[i].slope, hits_vec[i].slope_e);
  }

  TCanvas* canvas1 = new TCanvas("canvas1", "canvas1");
  canvas1->Divide(2, 2);
  canvas1->cd(4);
  slope_histo->Draw();
  canvas1->cd(1);
  slope_e_histo->Draw();
  canvas1->cd(2);
  slope_vs_slope_e_histo->Draw("colz");

  for (int i = 0; i < hits_vec.size(); i++) {
    pedestal_histo->Fill(hits_vec[i].pedestal);
    pedestal_e_histo->Fill(hits_vec[i].pedestal_e);
    pedestal_vs_pedestal_e_histo->Fill(hits_vec[i].pedestal,
                                       hits_vec[i].pedestal_e);
  }

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->Divide(2, 2);
  c2->cd(4);
  pedestal_histo->Draw();
  c2->cd(1);
  pedestal_e_histo->Draw();
  c2->cd(2);
  pedestal_vs_pedestal_e_histo->Draw("colz");

  for (int i = 0; i < hits_vec.size(); i++)
    slope_vs_pedestal_histo->Fill(hits_vec[i].slope, hits_vec[i].pedestal);
  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->cd(1);
  slope_vs_pedestal_histo->Draw("colz");

  
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


  signalFile->Close();
  pedestalFile->Close();
}
