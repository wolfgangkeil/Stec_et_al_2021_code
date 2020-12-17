# Code accompanying Stec et al. An Epigenetic Priming Mechanism Mediated by Nutrient Sensing Regulates Transcriptional Output during C. elegans development, Current Biology, 2020

This repository contains whole-animal fluorescence traces obtain from long-term live imaging as well as MATLAB® scripts to replicate the analyses of pulsatile gene expression, performed with the long-term live imaging data in Stec et al. Current Biology (2020)

If you use parts of this code please cite: Stec et al., Current Biology (2020), https://doi.org/10.1016/j.cub.2020.11.060
For questions, please contact wolfgang.keil@curie.fr

## Data
All data is contained in the subfolder 'data/'. Subfolders within 'data/' are named according to the C. elegans strain names (e.g. GR1395, HML474, see Key Resources Table in STAR Methods of Stec et al  Curr Biol 2020).  The 'data/' folder also contains .txt files (e.g. 'HML474_list.txt', 'all_WT_list.txt') which specify lists of individual animals that can be analyzed or plotted together. 

Within the strain folders, for each animal, information about experimental conditions and the molting times as scored manually from the imaging data, is contained in .txt files, e.g. 'data/HML474/HML474_02-Feb-2018_2_2.txt'. These files are parsed by the function 'read_single_worm_molting_data' (see below). The .mat files with the same filenames contain the fluorescence time traces as well as information about background and parameter values of the Gaussian peaks fitted to the background-subtracted fluorescence time traces.  

## Plotting individual fluorescence time courses for pulsatile gene expression together with Gaussian fits
Start by executing  'plot_worm_overview('data/', 'HML474','HML474_list.txt', 'GFP', 15,'constant')' as an example. You will see a plot showing all the fluorescence time traces
overlaid with the corresponding Gaussian fits.  By changing the strain and and list name, you can plot different genetic backgrounds or fluorescent reporter genes.

## Analysis of dynamic features (onset, peak phase, duration etc.) in each larval stage
Start by executing  'plot_peak_statistics_comparison({'data/', 'data/'}, {'GR1395', 'HML274'},{'GR1395_list.txt','HML274_list.txt'},{'WT','blmp-1(0)'}, 'mlt-10');' as an example. This will replicate figure 5J of Stec et al. 2020. By changing the strain and and list names, you can compare different genetic backgrounds for the same reporter gene to each other.

## Installation

Download by 

$ git clone  https://github.com/wolfgangkeil/Stec_et_al_2020_code.git


## License
Copyright (c) [2020] [Wolfgang Keil]

This repositry contains free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

All code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details at <https://www.gnu.org/licenses/>.
