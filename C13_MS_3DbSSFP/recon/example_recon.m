clear all, clc

% addpath
addpath('../orchestra-sdk-2.1-1.matlab'); % add orchestra path
rawdata_path = '../Exam258/Series23/';  % add rawdatapath
load('../ktraj_C1_0324_3DbSSFP.mat'); % load ktrajectory info.: k & dcf

rawdata = read_raw_ScanArchive(rawdata_path);

ncoil = size(rawdata, 2);
nfov = 32;
slice = 16;
nmet = 3;
np = size(rawdata, 3);
nshot = 5;
nT = 1;
kz_list = cal_kz(slice, 1)+floor(slice/2)+1;
k = reshape(k, [length(k)/nshot nshot]);
dcf = reshape(dcf, [length(dcf)/nshot nshot]);
k = repmat(k, [1 ceil(np/nshot)]); 
dcf = repmat(dcf, [1 ceil(np/nshot)]); 
data = rawdata;
ims = zeros(nfov,nfov,ncoil,np);
for i = 1:np
    ims(:,:,:,i) = gridkb_batch(data(:,:,i),k(:,i),dcf(:,i),nfov,1.2,4.5);
end
ims = reshape(ims,[nfov nfov ncoil nshot slice nmet nT]);
ims(:,:,:,:,kz_list,:,:) = ims;
ims1 = fftc(ims,5);
ims1 = sum(ims1, 4);
% coil combine
ims_sos = squeeze(sqrt(sum(abs(ims1).^2,3)));
% AUC image
ims_auc = sum(ims_sos, 5);
figure; imagesc(abs(reshape(ims_auc(:,:,floor(slice/2),:),[nfov nfov*nmet]))); axis image;


function kz = cal_kz(nslice, type)
% type: 1-center out; 2-outside in; 3-last to first
    if ~exist('type', 'var')
	    type = 1;
    end
    kz = zeros(nslice, 1);
    switch type
        case 1
            kz(1:2:end) = (floor(nslice/2+1):nslice) - floor(nslice/2)-1;
            kz(2:2:end) = (floor(nslice/2):-1:1) - floor(nslice/2)-1;
        case 2
            if mod(nslice, 2) == 1
                kz(1:2:end) = (nslice:-1:floor(nslice/2+1)) - floor(nslice/2)-1;
                kz(2:2:end) = (1:floor(nslice/2))- floor(nslice/2)-1;
            else
                kz(1:2:end) = (nslice:-1:floor(nslice/2+1)) - floor(nslice/2);
                kz(2:2:end) = (1:floor(nslice/2))- floor(nslice/2);
            end
        case 3
            kz = (nslice:-1:1) - floor(nslice/2)-1;
        otherwise
            error('type out of range');
    end
end