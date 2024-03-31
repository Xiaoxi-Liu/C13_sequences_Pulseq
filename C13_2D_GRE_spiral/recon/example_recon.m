clear all, clc

% addpath
addpath('../orchestra-sdk-2.1-1.matlab'); % add orchestra path
addpath('../utils'); % http://web.stanford.edu/class/ee369c/mfiles/gridkb.m
rawdata_path = '../Exam257/Series3/';  % add rawdatapath
load('../ktraj_C1_0322_BicCo2.mat');  % load ktrajectory info.: k & dcf

rawdata = read_raw_ScanArchive(rawdata_path);

ncoil = size(rawdata, 2);
nfov = 32;
slice = 5;
nmet = 3;
np = size(rawdata, 3);
necho = [10, 5, 5];
nT = 20;
maxecho = max(necho);
 
data = rawdata;
ims = zeros(nfov,nfov,ncoil,np);
for i = 1:np
    ims(:,:,:,i) = gridkb(data(:,:,i),k,dcf,nfov,1.2,4.5);
end
ims = reshape(ims,[nfov nfov ncoil slice*sum(necho) nT]);
ims1 = zeros(nfov, nfov, ncoil, maxecho, slice, nT, nmet);
num = 0;
for i = 1:nmet
    ims1(:,:,:,1:necho(i),:,:,i) = reshape(ims(:,:,:,num+1:sum(necho(1:i))*slice,:),[nfov nfov ncoil necho(i) slice nT]);
    num = sum(necho(1:i))*slice;
end
% AUC image
ims_auc = sum(ims1, 6);
% coil combine
ims_sos = squeeze(sqrt(sum(abs(ims_auc).^2,3)));
figure; imagesc(abs(reshape(ims_sos(:,:,1,floor(slice/2)+1,:),[nfov nfov*nmet]))); axis image;


