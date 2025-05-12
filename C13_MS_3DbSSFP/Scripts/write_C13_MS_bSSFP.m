function write_C13_MS_bSSFP(para, sys, export_filename, RF_path, RO_path)
%
% Metabolite-specific 3D bSSFP sequence for C13 imaging
% 
% Authors: Xiaoxi Liu, Di Cui,  
% Xiaoxi.Liu@ucsf.edu, Di.Cui@ucsf.edu, 03/24/2024
% 

for i = 1:para.nmet
    path = dir([RF_path '*' para.met_name(i,:) '*.mat']); 
    load(path.name);
    rf2p = rf2pulseq(rf, dwell, sys.rfRasterTime, sys);
    rf_bssfp = mr.makeArbitraryRf(rf2p, 2*pi*abs(sum(rf2p*sys.rfRasterTime)), 'system', sys);
    rf_ex{i} = rf_bssfp;
    para.TR(i) = TR;
end

fov = para.res*para.nx;
fovmin = min(fov);
if ~exist('RO_path', 'var')
    [k,dcf,t,ind,out,grad] = design_spiral_pulseq(fovmin,para.nx,para.interleave,sys.gradRasterTime*1e6,false,sys.maxGrad,sys.maxSlew,true,false,true);
    grad = grad.*sys.gamma; 
    gradn = grad(:,:,1)+1i*grad(:,:,2);
else
    load(RO_path);
    grad = grad*res/min(para.res);
    gradn = grad.*sys.gamma;
end

adcTime = sys.gradRasterTime*size(k,2)/para.interleave;
adcSamples = size(k,2)/para.interleave;
adc = mr.makeAdc(adcSamples,'Duration',adcTime,'delay',sys.adcDeadTime+1.2e-5,'system', sys);
for n = 1:para.interleave
    gx{n} = mr.makeArbitraryGrad('x', squeeze(real(gradn(:,n))), 'delay', sys.adcDeadTime, 'system', sys);
    gy{n} = mr.makeArbitraryGrad('y', squeeze(imag(gradn(:,n))), 'delay', sys.adcDeadTime, 'system', sys);
end
gz = mr.makeTrapezoid('z', sys, 'Area', 1/para.sth, 'delay', 8e-4);

delayTR = zeros(1, para.nmet);
for n = 1:para.nmet
    delay = ceil((para.TR(n) - mr.calcDuration(gx) - mr.calcDuration(gz) - mr.calcDuration(rf_ex{n}))/sys.gradRasterTime)*sys.gradRasterTime;
    delayTR(n) = delay - mod(delay, 0.4e-5);
end
assert(all(delayTR > 0));

delay1 = ceil((para.tem_res-(length(para.cat_FA)*2+para.nz*para.interleave)*sum(para.TR))/sys.gradRasterTime/para.nmet)*sys.gradRasterTime;
delay1 = delay1 - mod(delay1, 0.4e-5);
assert(delay1 > 0);

seq = mr.Sequence(sys);
rf_ex_met = rf_ex;
gz_met = gz;
gz_repmet = gz;
gz_repmet.amplitude = -1*gz.amplitude;
gz_repmet.delay = mr.calcDuration(adc);
kz_order = cal_kz(para.nz,1);

for ntp = 1:para.nT
    for nmet = 1:para.nmet
        gz_met.delay = mr.calcDuration(rf_ex{nmet}) + gz.delay;
        rf_ex_met{nmet}.freqOffset = para.freqshift(nmet);
        adc.freqOffset = para.freqshift(nmet); 
        for ncat = 1:length(para.cat_FA)
            rf_ex_met{nmet}.signal = para.cat_FA(ncat)*para.FA(nmet)/90*rf_ex{nmet}.signal;
            rf_ex_met{nmet}.phaseOffset = mod(ncat,2)*pi;
            seq.addBlock(rf_ex_met{nmet}, mr.scaleGrad(gz_met,0), mr.makeLabel('SET', 'TRID', 1+(nmet-1)*(para.interleave+2)));
            seq.addBlock(mr.scaleGrad(gx{1},0),mr.scaleGrad(gy{1},0));
            seq.addBlock(mr.makeDelay(delayTR(nmet)));
        end
        res = mod(ncat,2);
        for nslc = 1:para.nz
            for ni = 1:para.interleave
                rf_ex_met{nmet}.phaseOffset = mod((nslc-1)*para.interleave+ni+res,2)*pi;
                adc.phaseOffset = mod((nslc-1)*para.interleave+ni+res,2)*pi;
                seq.addBlock(rf_ex_met{nmet}, mr.scaleGrad(gz_met,kz_order(nslc)/floor(para.nz/2)), mr.makeLabel('SET', 'TRID', ni+1+(nmet-1)*(para.interleave+2)));
                seq.addBlock(mr.scaleGrad(gx{ni},fovmin/fov(nmet)),mr.scaleGrad(gy{ni},fovmin/fov(nmet)),adc,mr.scaleGrad(gz_repmet,kz_order(nslc)/floor(para.nz/2)));
                seq.addBlock(mr.makeDelay(delayTR(nmet)));
            end
        end
        res = mod((nslc-1)*para.interleave+ni+res,2);
        for ncat = 1:length(para.cat_FA)
            rf_ex_met{nmet}.signal = para.cat_FA(length(para.cat_FA)-ncat+1)*para.FA(nmet)/90*rf_ex{nmet}.signal;
            rf_ex_met{nmet}.phaseOffset = mod(ncat+res,2)*pi;
            if (ncat == length(para.cat_FA))
                seq.addBlock(rf_ex_met{nmet}, mr.scaleGrad(gz_met,0), mr.makeLabel('SET', 'TRID', nmet*(para.interleave+2)));
                seq.addBlock(mr.scaleGrad(gx{1},0),mr.scaleGrad(gy{1},0));
                seq.addBlock(mr.makeDelay(delayTR(nmet)+delay1));
            else
                seq.addBlock(rf_ex_met{nmet}, mr.scaleGrad(gz_met,0), mr.makeLabel('SET', 'TRID', 1+(nmet-1)*(para.interleave+2)));
                seq.addBlock(mr.scaleGrad(gx{1},0),mr.scaleGrad(gy{1},0));
                seq.addBlock(mr.makeDelay(delayTR(nmet)));
            end
        end
    end
end

[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.setDefinition('FOV', fov);
seq.setDefinition('Name', '3D-bSSFP');
seq.write([export_filename '.seq']);       

end

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
