# Stec_et_al_2020_code
This repository contains whole-animal fluorescence traces obtain from long-term live imaging as well as MATLAB® scripts to replicate the analyses performed in Stec et al. Current Biology (2020)

If you use this code please cite: Stec et al., Current Biology (2020)
For questions, please contact wolfgang.keil@curie.fr

## Peak fitting
Executing statistical_analyses_Z1Z4lineages.m will load all early somatic gonad cell division timing 
data gathered at 25˚C, provide some statistics and scatter plots.
Within statistical_analyses_Z1Z4lineages.m, you may change the variable list_file to
'WT_20degrees_list.txt' to analyze and plot data gathered at 20˚C.
 
The subfolders
fluorescence_data/mlt-10 & 
fluorescence_data/zk180.5
contain the files with whole animal fluorescence of individual animals and .txt files specifying overall developmental progression of the animals (molts etc.)

## Analysis of stage durations and scaling of fluorescence time traces

## Analysis of dynamic features (onset, peak phase, duration etc.) in each larval stage


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
