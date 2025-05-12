function write_2DGRE_C13(para, RF_path, RO_path, sys, export_filename)
%
% 2D single-echo/multi-echo GRE sequence for C13 imaging
% 
% Authors: Xiaoxi Liu, Di Cui,  
% Xiaoxi.Liu@ucsf.edu, Di.Cui@ucsf.edu, 03/2024
% 
fov = para.res*para.nx;
fovmin = min(fov);

%% load RF waveforms
load(RF_path);

rf2p = rf2pulseq(rf, dwell, sys.rfRasterTime, sys);
rf_ex = mr.makeArbitraryRf(rf2p, 2*pi*abs(sum(rf2p*sys.rfRasterTime)), 'system', sys);

Gss = g.*sys.gamma.*maxslc/para.sth; % T/m -> Hz/m
gzRF = mr.makeArbitraryGrad('z', squeeze(Gss), 'delay', sys.rfDeadTime, 'system', sys);

%% create spiral waveform
if isempty(RO_path)
    [k,dcf,t,ind,out,grad] = design_spiral_pulseq(fovmin,para.nx,1,sys.gradRasterTime*1e6,false,sys.maxGrad,sys.maxSlew,true,true,true);
else
    load(RO_path);
    grad = grad*res/min(para.res);
end
grad = grad.*sys.gamma; % T/m -> Hz/m

%% create ADC event and gradients
adcTime = sys.gradRasterTime*size(k,2);
adcSamples = size(k,2);
adc = mr.makeAdc(adcSamples,'Duration',adcTime,'delay',sys.adcDeadTime+1.2e-5,'system', sys);
gx = mr.makeArbitraryGrad('x', squeeze(grad(:,1)), 'delay', sys.adcDeadTime, 'system', sys);
gy = mr.makeArbitraryGrad('y', squeeze(grad(:,2)), 'delay', sys.adcDeadTime, 'system', sys);
gzSpoil = mr.makeTrapezoid('z', 'Area', 10/para.sth, 'system', sys);

%% Timing calculation
if max(para.necho)>1
    delayTE = ceil((para.deltaTE - mr.calcDuration(gx))/sys.gradRasterTime)*sys.gradRasterTime;
    assert(all(delayTE > 0));
else 
    delayTE = 0;
end

delayTR = zeros(1, para.nmet);
for n = 1:para.nmet
    delay = ceil((para.TR(n) - mr.calcDuration(gx)*para.necho(n) - mr.calcDuration(gzSpoil) - ...
        mr.calcDuration(gzRF) - delayTE*(para.necho(n)-1))/sys.gradRasterTime)*sys.gradRasterTime;
    delayTR(n) = delay - mod(delay, 4e-6);
end
assert(all(delayTR > 0));

%% Assemble sequence
seq = mr.Sequence(sys);
rf_phase = 0;
rf_inc = 0;
rfSpoilingInc = 117; 
rf_ex_met = rf_ex;
phase = 0*rf_ex.signal;
for i = 1:length(rf_ex.signal)
    phase(i) = dwell*sum(Gss(1,1:i));
end

for ntp = 1:para.nT
    for nmet = 1:para.nmet
        for nslc = 1:para.nz
            % excitation
            rf_ex_met.phaseOffset = rf_phase/180*pi;
            rf_ex_met.signal = abs(rf_ex.signal)*para.FA(nmet)/90.*exp(1i*(phase*2*pi*ceil((nslc-para.nz/2-1))*para.slc_sep+angle(rf_ex.signal)));
            adc.phaseOffset = rf_phase/180*pi;
            rf_ex_met.freqOffset = para.freqshift(nmet);
            adc.freqOffset = para.freqshift(nmet);
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            seq.addBlock(rf_ex_met, gzRF, mr.makeLabel('SET', 'TRID', nslc+(nmet-1)*para.nz));

            % readout
            for echoes = 1:para.necho(nmet)
                seq.addBlock(mr.scaleGrad(gx,fovmin/fov(nmet)),mr.scaleGrad(gy,fovmin/fov(nmet)),adc);
                if echoes ~= para.necho(nmet)
                    seq.addBlock(mr.makeDelay(delayTE));
                end
            end

            seq.addBlock(gzSpoil);
            seq.addBlock(mr.makeDelay(delayTR(nmet)));
        end
    end
end

fprintf('\n');

%% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Write .seq file
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'C1_2D_ME-GRE');
seq.write([export_filename '.seq']);       

end
