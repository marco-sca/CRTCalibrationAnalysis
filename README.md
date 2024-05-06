# Calibration of the Top Cosmic Ray Tagger in ICARUS

## Abstract
The ICARUS detector, a Liquid Argon Time Projection Chamber (LArTPC), is taking data at Fermilab on the Booster's neutrino beam.
Being located at ground level, the TPC is exposed to a high flux of cosmic muons, whose signal constitutes a background in the reconstruction of neutrino interactions. The rejection of the signal induced by the passage of cosmic rays in the TPC is possible thanks to the Cosmic Ray Tagger (CRT), an external subdetector divided into three subsystems: Top, Side and Bottom CRT. The following work focuses on the calibration of the Top CRT, consisting of 123 modules that act as hodoscopes using plastic scintillator bars. The readout of the light generated by the passage of a cosmic ray is carried out through Silicon Photomultipliers (SiPMs). The developed code has the following objectives:
- Estimate the pedestal (noise) and SiPM gain values of each Top CRT channel;
- Develop the analysis and implement it in icaruscode , integrating it with the experiment pipeline.
The analysis campaign, in the CRT raw data decoding stage, runs on the files produced through the Production Operations Management System (POMS). This software allows to launch, modify and monitor large scale campaigns of data processing jobs and for an effective calibration, the analysis of at least 50k events is needed (~10 mln CRT hits).
<!-- Aggiungi una brevissima descrizione di cosa fanno, a cosa servono, i due file di codice della repository -->

# 1. INTRODUCTION
## 1.1 The SBN Program at Fermilab
The Short-Baseline Neutrino (SBN) program at Fermilab is a research initiative aimed at investigating the potential existence of sterile neutrinos at the eV mass-scale. Sterile neutrinos, if they do exist, are challenging to observe directly because they don't interact with regular matter through the weak nuclear force. However, their presence could lead to new oscillations among the known neutrino flavors. The search for light sterile neutrinos at SBN is motivated by a set of anomalous results in past neutrino data, most significantly from the LSND and MiniBooNE experiments.

<a id="fig1"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/5f106f3e-7a8d-423a-9063-c1132bef80ef" alt="SBN Program detectors" width="500">
  <br>
  <p><strong>Figure 1.</strong> Illustration of the three SBN detectors along the Booster Neutrino Beam.</p>
</div>

To conduct this research, as shown in the [Figure 1](#fig1), three Liquid Argon Time Projection Chambers (LArTPCs) are strategically positioned along the Booster Neutrino Beamline (BNB). These detectors include SBND, which serves as the near detector and is located 110 meters from the neutrino source, MicroBooNE, and ICARUS, which acts as the far detector. These last two detectors are positioned 470 meters and 600 meters away from the source, respectively. ICARUS has been refurbished and upgraded to maximize its performance within the SBN program. After its overhauling at CERN, ICARUS was shipped to Fermilab in 2017 and since 2022 has been actively collecting data. The joint work of SBND and ICARUS creates a world-leading sterile neutrino search experiment that can cover the parameters allowed by past anomalies at $\geq 5\sigma$ significance.

At Fermilab the ICARUS detector is exposed to a significant influx of cosmic particles during the brief neutrino beam-spill periods (1.6 µs for BNB and 9.6 µs for NuMI) as well as during the approximately 1 millisecond TPC drift time, which falls outside the beam-spill. To minimize the impact of cosmic ray-induced events, ICARUS is equipped with a Cosmic Ray Tagger (CRT) system that ensures comprehensive coverage of the detector in all directions.
 
Here I will outline the tasks accomplished during my 2-month internship at Fermilab, as part of the "2023 Summer Students Italian program at the Fermi National Accelerator Laboratory and at other US Laboratories". After providing a brief overview of the SBN program and the theoretical and phenomenological aspects related to the search for sterile neutrinos, I will delve into the details of the ICARUS detector and the Cosmic Ray Tagger (CRT) system. My primary focus will be on the top section of the cosmic rays detector, known as the Top CRT, and I will elaborate on the calibration analysis process, which has been the central focus of my work during this internship.

<details>
  <summary>ICARUS and the Cosmic Ray Tagger</summary>

## 1.2 The SBN Far Detector: ICARUS
The ICARUS-T600 detector with an active mass of 476 tons of liquid argon has been the first large-scale operating LArTPC detector. ICARUS (Imaging Cosmic And Rare Underground Signals) consists of two adjacent modules of $3.6 m \times 3.9 m \times 19.9 m$ filled with a total mass of 760 tons of liquid argon, purified by removing the electronegative impurities. Each module is composed of two LAr-TPCs, separated by a common cathode made of a stainless steel frame structure supporting punched stainless-steel sheets. The anode and the cathode planes have a maximum drift lenght of 1.5 m, corresponding to $\sim 0.96$ ms drift time at the nominal 500 V/cm electric drift field. The anode plane is composed of three parallel wire planes 3 mm apart and oriented at different angles: the first with horizontal wires and the other two at $\pm 60^\circ$ from the horizontal direction (see [Figure 2](#fig2)). The optical system is composed of PMTs located behind the anodic wire planes, to collect the scintillation light used to generate the global event trigger.

<a id="fig2"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/1d65c106-01f8-4683-9440-cd96b1b76843" alt="ICARUS overhauling" width="400">
  <br>
  <p><strong>Figure 2.</strong>Picture of the ICARUS TPC during the CERN overhauling. Cathode (left), field shaping electrodes (top and bottom) and PMTs (right) are visible.</p>
</div>

The detector has been operating for 3 years (2011-2013) in the Gran Sasso Laboratory in Italy (LNGS). After that, in 2014 the ICARUS detector was transported to CERN and underwent a significant overhauling.  The two ICARUS modules have then been transported to Fermilab in July 2017 and ICARUS was installed in the SBN far detector building in August 2018.

## 1.3 Cosmic Background
The ICARUS-T600 detector was initially designed to operate in the low muon cosmic background of the Gran Sasso laboratory. The conditions at FNAL are completely different: placed just below the surface the detector is subject to a significant cosmic ray background and this may induce several additional and uncorrelated triggers during the $\sim 1$ ms drift time. Simulations showed that the expected rate of cosmics depositing more than 100 MeV within the T600 active volume is of $\sim 11$ kHz. Cosmic particles entering the detector during the $1.6 \mu s$  BNB neutrino beam-spill interact in the liquid argon generating scintillation light and an event trigger, the so-called \textit{in time activity}. The \textit{out of time} cosmic activity corresponds to cosmic muons crossing the detector during the $\sim 1 ms$ TPC drift time. On average $\sim 11$ cosmic tracks are expected over the full T600 volume during the drift window, generating a background that has to be disentangled from the neutrino event tracks. One of the most important sources of background to the $\nu_e$ appearance analysis is due to electromagnetic showers induced by $\gamma$ produced by cosmic particles propagating through the detector and in the surrounding materials. By showering withing the active liquid argon volume, the cosmogenic photon can mimic a genuine $\nu_e$ CC interaction.  Without systems in place to mitigate cosmic rays, the detector would be unable to effectively conduct any meaningful search. In order to mitigate the cosmogenic induced background, the ICARUS T600 detector is indeed surrounded with an external Cosmic Ray Tagger system (CRT) below a 3 m concrete overburden (6 m water equivalent). The CRT system is described in the following sections.

## 1.4 The CRT
The CRT system serves as an external subdetector located outside the cryostats, and its primary purpose is to identify charged particles that pass through or come close to the active volume of the TPC. With both the PMT and CRT systems offering an expected time resolution of a few nanoseconds, their synchronization and synergy allows for the determination of the direction of detected particles using the Photodetection system (PMT). This allows discrimination between events coming from the outside the detector from those generated inside and therefore rejecting cosmic ray induced triggers. Through a precise timing calibration effort, it becomes possible to filter out events in which the initial trigger was triggered by an identified cosmic particle entering the detector.

The CRT system encompasses an area of approximately 1100 square meters and is divided into three distinct subsystems: the \textit{Top CRT}, \textit{Side CRT}, and \textit{Bottom CRT}. These subsystems complement each other, ensuring complete coverage ($4\pi$) of the active LAr volume and enabling the identification of nearly 95$\%$ of passing through cosmics. In [Figure 3](#fig3) a representation of the Top and Side CRT sub-systems from the beam perspective.

<a id="fig3"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/9a0106e1-15ca-41b3-9c15-33fa8a369cbd" alt="Top and Side CRT representation">
  <br>
  <p><strong>Figure 3.</strong> Representation of the Top and Side CRT sub-systems.</p>
</div>

### 1.4.1 The Top CRT
The Top CRT is designed to capture around 80$\%$ of the cosmic muons that enter the ICARUS LArTPC. It consists of 123 modules, with 84 modules placed on the top horizontal plane and 39 modules covering the upper perimeter of the TPC (vertical rims). You can view an image of the Top CRT in [Figure 5](#fig5) taken from the ground floor of the Far Detector Building at Fermilab, before the concrete overburden was installed. These modules function as hodoscopes and are composed of two perpendicular layers, each containing 8 scintillator bars, which are 23 cm wide. These scintillator bars are enclosed in aluminum boxes measuring 1.86 meters $\times$ 1.86 meters, as depicted in [Figure 4](#fig4) (left picture) below.

<a id="fig4"></a>
<div align="center">
  <!-- Image 1 -->
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/d911f8c7-532e-4339-9466-efcdc4a51681" alt="Top CRT module" width="300">
  <!-- Image 2 -->
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/e7cae507-fa10-4163-a0b6-f1d013a516bf" alt="Channels" width="300">
  <!-- Image 3 -->
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/6826f2c2-2b27-4b86-840d-c76eeb6c312c" alt="Silicon Photomultiplier" width="300">
  <br>
  <p><strong>Figure 4.</strong> <i> Left:</i> Sketch of a Top CRT module and its components. <i>Center:</i> Representation of the scintillator bar with the two fibers embedded along the longitudinal direction of the bar. <i>Right:</i> Picture of the SiPM connection scheme to the fiber.</p>
</div>

In the top layer, the scintillator bars are 10 mm thick, while in the bottom layer, they are 15 mm thick. Each scintillator strip in the Top CRT has two WLS fibers embedded along the length of the bar, positioned 6 cm from each side, as shown in [Figure 4](#fig4) (center). These fibers are read-out from only one end, with the opposite end mirrored to enhance the light yield. A Hamamatsu S13360-1350CS SiPM is used for light detection, and the coupling of the WLS to SiPM is illustrated in [Figure 4](#fig4) (right). The system has a crosstalk probability of approximately 3$\%$ and a photon detection efficiency of around 40$\%$ at 450 nm \cite{Poppi:phd}.

<a id="fig5"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/cd890a96-2053-4e8f-bdd9-41e8e85c098d" alt="Top CRT installation" width="200">
  <br>
  <p> <strong>Figure 5.</strong> Picture of the fully installed Top CRT, before the OB installation. </p>
</div>

The SiPMs in each module are read out and biased by their respective Front End Board ([Figure 6](#fig6)), of the same type of those used for the Side CRT. The analog input signal is processed by a 32-channel ASIC (CITIROC\footnote{The Cherenkov Imaging Telescope Integrated Read Out Chip (CITIROC) is a 32 channel fully analog front-end ASIC dedicated to read-out of SiPMs}\cite{citroc}). These 32 signals are directed to an XILINX Spartan-6 FPGA chip, which handles basic input coincidence and triggering logic. Communication between the board and the host computer is facilitated through the Ethernet protocol.

<a id="fig6"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/50371b77-fd29-4f42-ad00-2bc59c449c49" alt="Front End Board picture" width="300">
  <br>
  <p> <strong>Figure 6.</strong> The Front End Board and its internal components. </p>
</div>

The primary function of the CRT modules is to accurately determine the precise position where muons cross through them. In the case of the Top CRT modules, they employ an XY scintillator layer configuration, enabling the creation of 64 coincidences of crossing strips (referred to as "sectors") within each module. You can see an example of a possible coincidence sector in [Figure 7](#fig7) when a cosmic muon passes through.

Each of the 32 channels is equipped with a CITIROC ASIC, which includes a charge amplifier with an adjustable gain and a dynamic range of 1 to 2000 photo-electrons (p.e.). The 32 trigger signals, denoted as C0 to C31, are in LVCMOS logic with a 3.3 V active state. These signals are directed to an FPGA, where they are combined using an AND logic operation to create coincidence signals for each of the two fibers from the same scintillator bar (the logic pairs the signals from even-odd channels, for example, C0$\And$C1,C2$\And$C3, and so on) so that if both fibers have detected light signals at the same time it indicates that a particle (such as a muon) has crossed that specific sector of the scintillator \cite{Poppi:phd}.

<a id="fig7"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/711673a9-c33a-443c-85a8-fb20a97a5177" alt="Coincidence sector" width="300">
  <br>
  <p> <strong>Figure 7.</strong> A possible coincidence sector at the passage of a cosmic muon. </p>
</div>

</details>

# 2. Calibration of the Top CRT
During my work at Fermilab I focused mainly on the implementation of a C++ macro for the calibration of the Top CRT channels. The analysis code is written in ROOT framework\cite{rootsite} and was implemented in icaruscode\footnote{\url{https://github.com/SBNSoftware/icaruscode}}\footnote{Version: v77$\_$00$\_$00}, integrating it with the experiment pipeline. 

The primary objective of the calibration analysis is to collect a large volume of CRT data over time to construct comprehensive response spectrums from the silicon photomultipliers (SiPMs). These SiPMs are responsible for detecting and converting light signals from the scintillator modules into electronic signals. To make sense of the data, we must translate the readings from the front-end boards (FEBs), initially in generic Analog to Digital Conversion (ADC) units, into measurements of photon energies. To achieve this conversion, the code seeks out for \textit{photopeaks} within the spectra. These photopeaks are characterized by their distinctive ADC values linked to the mean values of the peaks. Additionally, a calculated parameter called "peak number" is used to correlate and plot photopeaks against their respective ADC values. By applying a linear regression to multiple data points consisting of photopeak-ADC pairs on a graph, we can determine a conversion factor for linking ADC units to the number of detected photoelectrons (p.e.) in that particular channel. This calibration process is crucial for accurately interpreting the SiPM data in terms of the energy of detected particles or photons.
In the following section I present what the calibration analysis code does.

<details>
<Summary>Analysis code</Summary>
  
## 2.1 The calibration analysis
Goal of the analysis is to estimate the pedestal and gain values of each Top CRT channel. Those can be obtained by fitting the integrated ADC charge spectrum of each channel, exploiting the feature that at each trigger the FEB stores the ADC value of each of the 32 channels.

Before running the code, there is a decoding stage, where the raw data from each FEB are selected and converted into a readable format (decoding). The most relevant information of the CRT data product\footnote{\url{https://github.com/SBNSoftware/sbnobj/blob/develop/sbnobj/ICARUS/CRT/CRTData.hh}} that are used in the analysis are:
- the Front End Board MAC5 address (whose variable name is \textit{fMac5}), as mapped in [Figure 8](#fig8);
- the ADC values of all 32 FEB channels (\textit{fAdc[32]});
- the flag (\textit{fFlags}), that represents the CRT hit status.

<a id="fig8"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/6d3beaf3-0c1e-4d8c-a348-b48478cb1916" alt="Mac5 map" width="400">
  <br>
  <p> <strong>Figure 8.</strong> Map of the MAC5 addresses of the FEBs/modules of the Top CRT. </p>
</div>

The flag variable is an integer and can take the values: 3 if it is related to a regular CRT signal hit, 7 or 9 if it was a special reset hit. These reset hits are special events associated with a global trigger signal or a PPS (Pulse Per Second) signal regulated by the FEBs and generated by a White Rabbit system \cite{Poppi:phd}.

After the decoding stage, the data entries are used to construct integrated ADC spectra for all the channels. This process involves iterating through the dataset and, for each entry corresponding to a CRT hit caused by a cosmic particle, recording the ADC value for each channel. This information is then used to generate histograms for each of the 32 channels within the 231 modules/FEBs of the Top CRT. Each histogram is named "hadc$\_$channel$\#\_$feb$\#$," where "channel$\#$' and "feb$\#$" are replaced with the specific channel and FEB numbers, respectively. An example of spectrum is showed in [Figure 9](#fig9).
Additionally, the same dataset is employed to create histograms representing the inherent electronic noise in each channel, commonly referred to as the \textit{pedestal}. Furthermore, histograms are generated to display the ADC distribution of signal hits as well as histograms displaying the sum of the ADC counts of all the channels, for each module/FEB ([Figure 10](#fig10)).

<a id="fig9"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/948a78ef-72d1-4deb-bc36-33d44beefdfe" alt="channel spectrum example" width="400">
  <br>
  <p> <strong>Figure 9.</strong> Example of spectrum for a 15 mm scintillator channel (Top Layer) zoomed in the range 0 – 1100 ADC Counts. The pedestal and signal peaks are visible. </p>
</div>

<a id="fig10"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/8717f603-5af6-4f21-b659-7e8a5fd469aa" alt="Example all signal sum on channels" width="400">
  <br>
  <p> <strong>Figure 10.</strong> Histogram obtained by the sum of the ADC count values of all the 32 channels of the FEB with mac5 address 136. </p>
</div>

### 2.1.1 Signal and pedestal selection
To generate the channels' signal spectra I selected the data entries\footnote{Data from calibration run 9989 of 6/26/2023 ($\sim$ 23 hours)}\footnote{In order to increase the statistic of CRT Hits, mainly for calibration purposes, the acquisition window for the Top CRT was set to $\pm$ 25 ms w.r.t. trigger timestamp. The window is extremely larger than the $\sim$ 2 ms drift window of the TPC, this is why the CRT hits sample of the 50 ms window is used only for calibration purposes, the normal data flow of event reconstruction uses a software reduced window of $\pm$ 3 ms.\cite{Poppi:phd}} with $fFlags = 3$ (the module recorded signal generated by a cosmic particle) excluding all ADC \textless 250 and took the highest 2 ADC counts in each layer of scintillator bars. These values were used to fill the histograms of the corresponding channels.

Initially the pedestal distributions were derived by analyzing the ADC values recorded in the channels of each layer with $fFlags = 3$, with the exclusion of the top 6 highest values per layer. In this way I had the ADC spectrum of a channel when it did not participate in the CRT hit channel coincidence (also referred to as non-triggering channel logic) \footnote{CRT triggering coincidence: Signal hits have at least 4 channels above the threshold, due to the internal trigger logic}. As can be observed in [Figure 11](#fig11) and in [Figure 12](#fig12) the distribution of the pedestal is larger then the average distance between the photoelectron peaks. This behaviour is not suitable to correctly estimate the waveform baseline for the pedestal.

<a id="fig11"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/dbd0822f-d124-46cc-8d4f-90df984853e4" alt="Noise spectrum from non-triggering channels" width="400">
  <br>
  <p> <strong>Figure 11.</strong> ADC spectrum of a top layer channel when it did not participate in the CRT hit channel coincidence. </p>
</div>

<a id="fig12"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/093ca8ea-63ea-4619-963c-726d9efeef9e" alt="Full top layer channel's spectrum" width="400">
  <br>
  <p> <strong>Figure 12.</strong> ADC spectrum of a top layer channel. </p>
</div>

I have then explored a different extraction method for the pedestal selection. By definition, random triggers of the CRT FEB should result in random values of each channel around its pedestal. Using the same dataset I exploited the T1 and T0 special reset events which behave as an external random triggers. Those hits are generated by an external uncorrelated source (Pulse per second signal or PMT trigger), so that the ADC value of all 32 channels are most likely electronic noise and a new sub-sample with a reduced statistic was obtained ([Figures 13](#fig13) and [Figure 14](#fig14)). The new distribution was considered to be more suitable for the pedestal evaluation.

<a id="fig13"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/24dedd25-9843-488b-8645-21d8dd6ba52b" alt="Example of pedestal from reset hits, log scale" width="400">
  <br>
  <p> <strong>Figure 13.</strong> Pedestal distribution for a Top Layer channel obtained from the reset hits, with a lower statistic and with the y-axis in log scale. </p>
</div>

<a id="fig14"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/d05d4e66-b9c1-4c12-ae7f-a3e58845edd9" alt="Example of pedestal from reset hits, no log scale" width="400">
  <br>
  <p> <strong>Figure 14.</strong> Pedestal distribution for a Top Layer channel obtained from the reset hits, with a lower statistic. </p>
</div>

A problem was observed when digitizing the special reset events of the T0/T1 counters: not all reset events were correctly identified and flagged as special events, but they were treated as regular signal hits (65$\%$ of the times the flag is correct \cite{Poppi:phd}). To solve this issue the sum of all the 32 ADC values for each hit can be used, in order to separate T1/T0 reset hits from signal ones. In [Figure 15](#fig15) we can see that the sum of the signal given by reset hits (red peak on the left) is superposed on a similar peak related to the sum of the ADC values of signal hits (in blue). In the next calibration analysis a cut for signal sum values over 7000 ADC will be tried for a better selection of reset hits.

<a id="fig15"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/fbe4c4db-bc37-4011-b881-be7e6be2c7c9" alt="Signal sum" width="400">
  <br>
  <p> <strong>Figure 15.</strong> Superposition of the distributions of adc values sum on all 32 channels of a FEB for pedestal obtained with reset hits correctly flagged (red) and signal hits (blue). </p>
</div>

Another viable option for the pedestal estimation, for future calibrations, is to look at the signal distribution of broken channels, where no signal is detected above the pedestal ([Figure 16](#fig16)).

<a id="fig16"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/febf01c9-3c9c-458f-9c6b-c018a17a7d39" alt="Broken channel spectrum" width="400">
  <br>
  <p> <strong>Figure 16.</strong> Example of a top layer broken channel's extracted signal. </p>
</div>

In conclusion,  the distribution for the pedesal is still very large, even if it is better if compared with the «Non triggering channels» extraction method of the data. Is also possible to notice the presence of two peaks in the distribution (circled in red in [Figure 17](#fig17)), where the left peak is generated by electronic noise when there’s a signal hit in other channels of the same FEB \cite{Poppi:phd} and the right peak could be a p.e. peak covered by the pedestal or SiPM intrinsic electronic noise. In an attempt to obtain a "clean" pedestal further investigation is required.

<a id="fig17"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/a70cf1f3-c2c4-45dd-bdb2-2d80358b4a7c" alt="Broken channel with strange peaks" width="400">
  <br>
  <p> <strong>Figure 17.</strong> Pedestal distribution with highlighted peaks, due to the presence of electronic noise.</p>
</div>

## 2.1.2 Analysis algorithm
Multiple gaussian fits are performed on the pedestals obtained through the reset hits, optimizing the fit range until the reduced $\chi^2$ is smaller than 10 or until there are no more bins in the selected range (see [Figure 18](#fig18)). The mean value extracted from the fit is then stored for the channel's gain evaluation.

<a id="fig18"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/bd45ae07-d2de-4121-a371-d424ae4106d3" alt="Pedestal Fit" width="400">
  <br>
  <p> <strong>Figure 18.</strong> Pedestal distribution for a Bottom Layer channel obtained from the reset hits (blue) and superposed gaussian fit (red) with $\chi^2$\textless 10. </p>
</div>

For the signal a similar procedure is followed: the ROOT function TSpectrum\cite{root} is used to search for the first 5 peaks in the hits spectrum, quantized photoelectron peaks are fitted recursively using a gaussian distribution, adjusting the fitted range of the histogram in order to minimize the reduced $\chi^2$. For each fit, the minimum distance between the previous and following peaks is used as the range and is recursively reduced until $\chi^2$ \textless 2 or until there are no more bins in the selected range. In [Figure 19](#fig19) a distribution of the charge spectrum for a bottom layer's channel is shown with overlayed the recursive gaussian fit of the first 5 photoelectron peaks. The mean and standard deviation values of the peaks, extracted from the fit, are then stored.

<a id="fig19"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/a097c2cc-64e7-4a4f-895e-e44ed9fb9ece" alt="Recursive Fit" width="400">
  <br>
  <p> <strong>Figure 19.</strong> Recursive single photoelectron peaks fitting with a gaussian distribution for a bottom layer's channel when participating in the CRT triggering coincidence. The signal has a cut for ADC counts \textgreater 250 and the left peak corresponds to 2 p.e. </p>
</div>

The gain estimation process relies on mean ADC values of detected peaks, plotted against their corresponding peak numbers. The gain for a specific channel is determined as the slope of a linear fit applied to this distribution. However, a peak is sometimes skipped by TSpectrum and some peaks are misidentified with others, introducing errors in the fit. To address this, an index rearrangement function was introduced to adjust the order of peak indices, exploiting a minimization of the reduced chi-squared ($\chi^2$) value of the fit. This function iteratively analyzes the spacing between adjacent peaks using this information  to adjusts their positions. The result of this work is presented in [Figure 20](#fig20), where an example of gain fit for a top layer's channel is shown. The colored band represents the growing sigma value of the fitted photo-electron peaks and the y-intercept is the ADC count mean value of the pedestal peak.

<a id="fig20"></a>
<div align="center">
  <img src="https://github.com/marco-sca/CRTCalibrationAnalysis/assets/140084724/47251686-4c99-437f-ba25-33e3fb81ce96" alt="Gain Fit" width="400">
  <br>
  <p> <strong>Figure 20.</strong> Distribution of the p.e. peaks mean value versus the corresponding p.e. number with the superposed linear fit (in red) used to evaluate the gain from the slope. The peak with index 1 is skipped. The blue band shows how the standard deviation of the fitted gaussians grows with the peak number. </p>
</div>

Following the calibration of pedestal and gains for all the Top CRT channels, the conversion of ADC counts to photo-electrons can be obtained by:
\begin{equation}
    n_{p.e.} = \frac{ADC_i - Ped_i}{G_i}
\end{equation}
where $n_{p.e.}$ is the resulting number of photo-electrons, $ADC_i$ is the ADC value of the i-th channel and $Ped_i$ and $G_i$ are, respectively, its pedestal and its gain as evaluated from the calibration. As future work, the average amount of light ("light yield") produced by the particles when they pass through scintillator bars will be determined. With the gain value and an adequate statistic we can obtain the distribution of the p.e. for each bar, fit and search for the peak (whose value represents the most probable number of p.e. produced per event) that is the average light yield for each channel.
</details>

<details>
<summary>Software and code framework</summary>
  
## 2.1.3 POMS
The calibration analysis campaign was run on the Production Operations Management System (POMS) \cite{poms} that allows to launch, modify and monitor large scale campaigns of data processing jobs. This was needed given the large scale of the analysis work: it was estimated that for an effective calibration of the Top CRT, at least 50 thousand events are needed, corresponding to $\sim$ 10 million CRT hits \cite{Poppi:phd}. POMS provides a web service interface that enables automated jobs submission on distributed resources according to customers’ requests and subsequent monitoring and recovery of failed submissions. Part of the calibration work included understanding the procedure to submit a POMS campaign gauged on my needs. Only the decoding stage was executed as a campaign stage in POMS and produced a substantial number of histograms. However, due to the nature of the decoding stage, which processes data file by file, each file containing information on approximately 50 PMT triggered events (around 10 CRT hits from cosmic rays within the data acquisition window), the resulting histograms had relatively few entries. Therefore, before running the calibration analysis code, I had to develop a script to merge a large number of ROOT files, enabling the creation of histograms with an higher number of entries.

## 2.1.4 LArSoft
<!-- To complete -->
</details>

# BIBLIOGRAPHY
<a id="references"></a> [1]
<a id="references"></a> [2]
