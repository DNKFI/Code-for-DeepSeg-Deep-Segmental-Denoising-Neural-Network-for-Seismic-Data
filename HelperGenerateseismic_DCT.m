function [targets1,predictors1] = HelperGenerateseismic_DCT(seismic,noisy, WindowLength,win,overlap,NumFeatures,NumSegments)

ffTLength    = WindowLength;

noisySeismic   = noisy;
cleanSTFT = stft_dct(seismic,'Window',win,'OverlapLength',overlap,'FFTLength',ffTLength,'Center',false);
cleanSTFT = (cleanSTFT  - mean(cleanSTFT(:))) / max(abs(cleanSTFT(:)));
cleanSTFT1 = (cleanSTFT(1:NumFeatures,:));
noisySTFT = stft_dct(noisySeismic,'Window',win,'OverlapLength',overlap,'FFTLength',ffTLength,'Center',false);
noisySTFT  = (noisySTFT   - mean(noisySTFT (:))) / max(abs(noisySTFT (:)));
noisySTFT1 = (noisySTFT(1:NumFeatures,:));
noisySTFTAugmented1   = [zeros(NumFeatures,ceil(NumSegments/2-1)) noisySTFT1 zeros(NumFeatures,ceil(NumSegments/2-1))]; 
STFTSegments1 = zeros( NumFeatures, NumSegments , size(noisySTFTAugmented1,2) - NumSegments + 1);
for index     = 1 : size(noisySTFTAugmented1,2) - NumSegments + 1
    STFTSegments1(:,:,index) = noisySTFTAugmented1(:,index:index+NumSegments-1);
end
targets1    = cleanSTFT1;
predictors1 = STFTSegments1;