close all
clear all
clc

% SetUp - add the libraries
addpath C:\Users\Jisuk\Documents\controlcont\rtc\matlab
%Load control-based continuation helper utilities.
addpath C:\Users\Jisuk\Documents\MATLAB\CBC_Flutter\cbc-master
addpath C:\Users\Jisuk\Documents\MATLAB\CBC_Flutter\PPCBC_FlutterRig\RTC_Code\flutter_phase_mean

%Create a control interface.
rtc = flutter_interface();

rtc.par.forcing_freq =1; % frequency
rtc.par.forcing_amp = 0; %amplitude of additional input force
rtc.par.x1_control = 0; %control on/off

rtc.opt.samples = 5e4 ; % recording sampling points 
rtc.datafields.stream_fields = {'aksim_angle','phi_bis', 'x1', 'x1_target','Fshaker' 'out','mean_h'}; % data that is stored, maximum 8 fields

% % constant for the laser
% a=47;
% c=-135;

% Get the indices of the fundamental harmonic
idx_sin = rtc.fourier.idx_fund(1);
idx_cos = rtc.fourier.idx_fund(2);

% Tolerance for the Picard iteration in the mean A_0 of the target A_0+A_1 sin(2 pi f t)
PicardTol=10^-4;
PicardTol_Freq=0.009;

% CONSTANTS in the beagle bone - they can be updated directly in the
% interface
rtc.opt.frequency=0.05;
%tolerance adopted for checking the steady condition
rtc.opt.x1_coeffs_var_tol_rel=0.1;
error_tol=0.01;
Wind_sped=16.4;
