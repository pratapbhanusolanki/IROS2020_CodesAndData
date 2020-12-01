# IROS2020_CodesAndData

## Description

This repository contains the complete code and data files to repeat the simulation results and generate figures that are presented in our paper titled "Active Alignment Control-based LED Communication for Underwater Robots" which is presented at, and will be published in the proceedings of, 2020 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), the submitted manuscript is included here for reference. The dependency npy-matlab-master, that is used to convert numpy files to matlab variables, is included with the original licence. 

## Instructions

The Data and scripts in this repo can be used to regenerate Figures, 7,8,9 and 12. 
The instructions and details to simulate the scenarios/read the experimental data and replicate the figures are as follows:
1. Figure 7: Run [StaticComparison_ES_TE_Fig7.m](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/StaticComparison_ES_TE_Fig7.m) using Matlab to simulate a run for the extremum seeking (ES) and the triangular exploration algorithms, and generate a figure similar to Figure 7. Due to randomness in the simulation, the plot may not exactly look the same but the qualitative illustration would be the similar.  

2. Figure 8: Run [MovingOptimaComparison_ES_TE_Fig8.m](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/MovingOptimaComparison_ES_TE_Fig8.m) simulate the algorithms with a moving optimum scenario and regenerate a figure similar to Figure 8.  

3. Figure 9: This figure reflects the statistical data from multiple simulation runs over a range of frequencies. The data from the simulation is stored in the file [DataFrequencyResponseFig9.mat](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/DataFrequencyResponseFig9.mat), which can be plotted by running the file [PlotDataFig9.m](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/PlotDataFig9.m) on Matlab. The complete simulation takes hours to run, the code of the simulation can be found in [SimulationMultipleMotionFrequenciesFig9.m](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/SimulationMultipleMotionFrequenciesFig9.m).  

4. Figure 12: The experiment was run in Python and hence the data is stored in .pkl and .npz/.npy fprmat. To regenrate the figure 12, please follow these instructions: 

	  a. Run [BER_calculation_Fig12.py](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/BER_calculation_Fig12.py) using python to generate [BitTransmissionDataIROSFig12.mat](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/BitTransmissionDataIROS_Fig12.mat) (This step can be skipped as the file BitTransmissionDataIROS_Fig12.mat is also included).
    
	  b. Using MATLAB run [data_analyzer_Fig12.m](https://github.com/pratapbhanusolanki/IROS2020_CodesAndData/blob/main/data_analyzer_Fig12.m), this will generate the figure corresponding to Figure 12 in the paper. 
