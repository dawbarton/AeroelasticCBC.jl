Data files for paper "Analysis of self-excited flutter oscillations with control-based continuation."

Data packaged on 24-08-2022 by Kyoung Hyun Lee (jz18526@bristol.ac.uk)

Contents
========

This package contains all the experimental data used in the paper, along with the
files necessary to generate the figures for the paper.

Description of the folders
	- CBC-master: Matlab files used in the CBC experiments
	- Linear_system_Identification: linear system identification code and the free decay data sets for model-1 and model-2.
	- Model1: CBC results for model-1.
	- Model2: CBC results for model-2.

Description of the files that generate figures
	1. Linear_system_Identification
		- Linear_system_ID.m: linear system identification using "system identification" toolbox for model-1. 
		- Linear_system_ID2.m: linear system identification using "system identification" toolbox for model-2
	
	2. Model1
		- save_data_plot.m: generates a bifurcation diagram, and the estimated flutter speed and frequency from the CBC results
	2. Model2
		- save_data_plot.m: generates a bifurcation diagram, and the estimated flutter speed and frequency from the CBC results


Descriptions of the experiments and the data

- CBC tests variables saved in Model1 and Model2 folder

	amplitude - noninvasive control target iteratively calculated (see cbc-master/search_non_invasive.m for detail)

	Wind_speed: wind speed

	actual_frequency: measured frequency of the LCO

	Array "data" was saved after noninvasive control was achieved for saving time-series

Note that bifurcation diagrams are plotted using (Wind_speed, amplitude). 

description of "data"
	data(1,:)- heave response in voltage (sensitivity = 18/100 m/V)
	data(2,:)- control target in voltage 
	data(3,:)- the force of the actuator in voltage 
	data(4,:)- pitch angle in rad 
	data(5,:)- phase signal of the CBC

Measurement files in Model1
=====
CBC test results:
CBC_stable_v14_9.mat
CBC_stable_v15_6.mat
CBC_stable_v16_5.mat
CBC_stable_v17_3.mat
CBC_unstable_v14_9.mat
CBC_unstable_v15_6.mat
CBC_unstable_v16_5.mat
CBC_unstable_v17_2.mat

Measurement files in Model2
=====
CBC test results:
CBC_v20_7.mat
CBC_v21_6.mat
CBC_v22_4.mat
CBC_v23_3.mat
CBC_v24_1.mat
CBC_v25_0.mat
CBC_v25_8.mat
CBC_stable_v20_7.mat
CBC_stable_v21_5.mat
CBC_unstable_v22_4.mat
CBC_unstable_v23_2.mat


The Matlab code that the experiment was running is (initialize.m, search_non_invasive.m).

All data were collected with a low-cost real-time control board based on the
BeagleBone Black. All hardware schematics and software are open source and
available from <http://github.com/~db9052/rtc> and
<http://github.com/~db9052/cbc>.








