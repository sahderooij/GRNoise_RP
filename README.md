# Reproduction Package for: 

### *Strong Reduction of Quasiparticle Fluctuations in a Superconductor due to Decoupling of the Quasiparticle Number and Lifetime*

*[doi to paper]*

This is a package of data and code to reproduce the results of the paper by S. A. H. de Rooij et al. (2021). The code in this package is also on [github](https://github.com/sahderooij/GRNoise_RP). The python code is a scaled down fork of [this repository](https://github.com/sahderooij/MKID-models). This document will explain the content of the package and the steps needed to reproduce the results. 

## Software prerequisites

- Python 3 environment with jupyter, numpy, scipy and matplotlib
- MATLAB (2020+) with curve fitting toolbox

## Package Contents

- **figures** -- raw figure files, which are presented in the paper
- **data** -- all the needed raw and processed data needed to reproduce the figures.
  - **A#/B#/C#** are subfolder containing the data for different devices as listed in the supplemental material. These have the following subfolders:
    - **Noise_vs_T/TD_2D** contains the raw Time Domain (TD) data, for different Kinetic Inductance Detectors (KIDs, i.e. the superconducting resonators) on the chip. Both read power and temperature are varied.
    - **NoiseTDanalyse** contains the file **TDresults.mat**, that stores the Power Spectral Densities (PSDs) calculated from the data in Noise_vs_T
    - **S21/2D**, contains the raw $S_{21}$ measurements (.dat-files) for different KIDs, with varying read power and temperature. It also contains .csv-files with extracted quality factors, responsivities and various other parameters, at different temperatures. The .png-files show summaries of $S_{21}$ analysis
    - *KIDparam.txt* is used in the $S_{21}$ analysis script and contains a table of design resonance frequencies, design Q-factors and lengths of the Al strips, for each KID.
- **notebooks** -- Jupyter notebooks that generate and save the figures
- **python_code** -- all necessary scripts for data analysis and plotting
  - Packages:
    
    - **kidata** is a package that deals with data from experiments. It contains the modules:
        - **io**, to find, read and write data
      - **calc**, which uses (multiple) raw data files to calculate Kinetic Inductance Detector (KID, which are the superconducting resonators) parameters, such as the kinetic induction fraction, phonon escape time and the quasiparticle lifetime from a Lorentzian fit to a PSD.
        - **filters**, to manipulate PSDs. For instance, to remove amplifier noise
      - **noise**, which processes the raw TD noise files: filter pulses and calculate PSDs, which are then stored in **NoiseTDanalyse/TDresults.mat**
        - **plot**, which contain (rather elaborate) plotting function to plot the data in various ways
    
  - Modules:
    
    - **kidcalc** implements superconductor theory that governs the behavior of MKIDs (mostly [BSC](https://link.aps.org/doi/10.1103/PhysRev.108.1175), [Mattis-Bardeen](https://link.aps.org/doi/10.1103/PhysRev.111.412) and [Kaplan et al.](https://link.aps.org/doi/10.1103/PhysRevB.14.4854))
      - **trapmodels** is used for trapping models, based on modified [Rothwarf-Taylor](https://link.aps.org/doi/10.1103/PhysRevLett.19.27) equations, to explain the observed reduction in noise level

    - Ddata:
    
      - **Ddata** contains calculated gap energies over temperature for different $T_c$, to speed up the calculation of $\Delta(T)$.
- **MATLAB_code** contains the code to process the $S_{21}$ measurement raw data and write the .csv-files with the extracted parameters, such as Q-factors and responsivities. 
  - *S21analysisV5.m*, the script to process the files
  - **subroutines** needed to run *S21analysisV5.m*

## Reproducing the Results

The Jupyter notebooks in the **notebooks** directory can be directly executed to obtain the figures presented in both the main text of the paper and the supplemental information. These notebooks mostly use the plot functions, in the *kidata.plot* module. 

In order for these to work, the raw data needs to be processed. The processed data files are present in the package, so that the notebooks work right away, but in order to reproduce these files two analyses must be done:

1. The PSDs need to be calculated from the raw TD noise files. This is done via the function *kidata.noise.do_TDanalysis()*, which takes the chipnumber (e.g. 'A1A2') as only required argument. This must be done for each chip seperately. The calculated PSDs are stored in **data/NoiseTDanalyse/TDresults.mat**
2. The KID parameters from the $S_{21}$ measurement need to be calculated. This is done with the MATLAB script *S21analysisV5.m* in the **MATLAB_code** folder. In order to run this properly, the *ChipInfo* variable should be set correctly at the top of the .m-file. This includes: Chipnumber in the path, $T_c$ and Al thickness and width. The scripts loads the *KIDparam.txt*- file, which is in the data folder of the chip, for the Al length for each KID. Furthermore, to calculate the correct responsivities, the *stopoffset* variable should be set such that the high temperature data points, where the $S_{21}$-fits go wrong, are not included in the analysis. This can be done by setting it to 0 first, after which the threshold can be determined from the generated .png-images and the script can be run again. The results are stored in .cvs-files in the same folder (e.g. **A1A2/S21/2D/KID2_99dBm_Tdep.csv**)







