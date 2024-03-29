%
% 2D single-echo/multi-echo GRE sequence for C13 imaging
% 
% Authors: Xiaoxi Liu, Di Cui,  
% Xiaoxi.Liu@ucsf.edu, Di.Cui@ucsf.edu, 03/2024
% 
clc; clear all;
% add paths for the pulseq and toppe functions
% change this to proper directories as needed
addpath '../PulCeq-main/matlab'
addpath(genpath('../pulseq-master'));
addpath(genpath('../toppe-main'));
export_filename = '../seq/C13_2DGRE_spiral';

%% Set system limits (update these values with configuration of the scanner)
sys = mr.opts('maxGrad', 50/sqrt(2), 'gradUnit','mT/m', ...
'maxSlew', 200/sqrt(2), 'slewUnit', 'T/m/s', ...
'rfDeadTime', 100e-6, ...
'rfRingdownTime', 60e-6, ...
'rfRasterTime', 4e-6, ...
'adcDeadTime', 20e-6, ...
'adcRasterTime', 2e-6, ...
'gradRasterTime', 4e-6, ...
'blockDurationRaster', 4e-6, ...
'gamma',10708000, ...
'B0', 3.0);

%% Set sequence parameters
para = [];  
para.nmet = 3; % number of metabolites
para.freqshift = [0, 395, -322]; % Hz, frequency shift for each metabolite
para.res = [5 10 15]*1e-3; % m, resolution for each metabolites
para.nx = 38;  % matrix size
para.sth = 21e-3;  % m, slice thickness
para.slc_sep = 21e-3;   % m, slice seperation, min = slice thickness
para.nz = 3;  % number of slice
para.necho = [15, 10, 10];  % echo number for each metabolite
para.deltaTE = 27e-3; % s, echo spacing
para.nT = 30; % number of time points
para.TR = [430, 285, 285]*1e-3; % s, Repeat Time for each metabolite
para.FA = [20, 30, 30]; % degree, flip angle for each metabolite


%% path for the RF and/or ReadOut gradient
RF_path = '../data/spsp_pulse_C1pyr.mat';
% RO_path = '../data/grad_ssspiral.mat'; % load customized readout gradient

%% build sequence
% write_2DGRE_C13(para, RF_path, RO_path, sys, export_filename);  % load customized readout gradient
write_2DGRE_C13(para, RF_path, [], sys, export_filename); % built-in spiral generation

%% write into GE file, set these values based on the configuration of scanner
sysGE = toppe.systemspecs('maxGrad', sys.maxGrad/sys.gamma*100, ...   % G/cm
    'maxSlew', sys.maxSlew/sys.gamma/10, ...                        % G/cm/ms
    'maxView', para.nmet*para.nz*para.nT, ...           % maximum is 2048
    'maxRF', 0.2, ...
    'rfDeadTime', 100, ... % us
    'rfRingdownTime', 60, ... % us
    'adcDeadTime', 20, ... % us  % 40
    'psd_rf_wait', 148, ... % RF/gradient delay (us)
    'psd_grd_wait', 156, ... % ADC/gradient delay (us)
    'gamma', sys.gamma); 

seq2ge([export_filename '.seq'], sysGE, [export_filename '.tar']);

