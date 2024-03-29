% write_MS_3DbSSFP.m
%
% MS-3DbSSFP sequence for C13 imaging
% 
% writen by Xiaoxi Liu, Di Cui,  
% Xiaoxi.Liu@ucsf.edu, Di.Cui@ucsf.edu, 03/24/2024
% 
clc; clear all;
addpath '/home/xliu18/Documents/pulseq/PulCeq-main/matlab'
addpath(genpath('/home/xliu18/Documents/pulseq/pulseq-master'));
addpath(genpath('/home/xliu18/Documents/pulseq/toppe-main'));

export_filename = 'C13_MS_3DbSSFP';
RF_path = '../Data/';
RO_path = '../Data/readout_grad.mat';

%% Set system limits
maxGrad = 50;
maxSR = 200;
sys = mr.opts('maxGrad', maxGrad, 'gradUnit','mT/m', ...
'maxSlew', maxSR, 'slewUnit', 'T/m/s', ...
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
para.nmet = 2; % number of metabolites
para.met_name = ['pyr'; 'lac']; % order of metabolites acquisition
para.freqshift = [0, 395]; % Hz, frequency shift for each metabolite
para.res = [10, 15]*1e-3; % m, resolution for each metabolites
para.nx = 45;  % matrix size
para.sth = 20e-3;  % m, slice thickness
para.nz = 16;  % number of slice
para.nT = 1; % number of time points
para.FA = [40, 60]; % degree, flip angle for each metabolite
para.interleave = 4;
para.cat_FA = [0.05, 0.12, 0.38, 0.69, 0.93, 1]; 
para.tem_res = 5; % s, temporal resolution

% write_C13_MS_bSSFP(para, sys, export_filename, RF_path, RO_path);  % load customized readout gradient
write_C13_MS_bSSFP(para, sys, export_filename, RF_path);  % built-in spiral generation

%%
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
