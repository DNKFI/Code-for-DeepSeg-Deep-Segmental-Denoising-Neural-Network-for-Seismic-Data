
close all
SNR=0;
Nt=20; % number of traces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load training_data.mat
y =R+z; %y=noisy traces, R=traces, z=noise
SNR_before_denoising=round(mean(10*log10(var(R)./var(y-R))))
windowLength = 128;
win          = ones(windowLength,1);%hamming(windowLength);
overlap      = round(.90* windowLength);%.90
ffTLength    = windowLength;
numFeatures  = ffTLength;%ffTLength/2 + 1;
numSegments  = 15;

Dshot=R;
Dshot_cn=y;
predictors =[];targets =[];cs=200;
seq=1:size(Dshot,1);
for nn=1:max(1,ceil(size(Dshot_cn,2)/cs))
    Temp1=num2cell(Dshot(seq,((nn-1)*cs+1):min(size(Dshot,2),cs*nn)),1).';
    T1=tall(Temp1);
    Temp2=num2cell(Dshot_cn(seq,((nn-1)*cs+1):min(size(Dshot,2),cs*nn)),1).';
    T2=tall(Temp2);
    [temp_targets,temp_predictors] = cellfun(@(x,y)HelperGenerateseismic_DCT(x,y,...
        windowLength,win,overlap,numFeatures,numSegments), T1,T2,"UniformOutput",false);
    [temp_targets,temp_predictors] = gather(temp_targets,temp_predictors);
    targets=[targets;temp_targets];
    predictors=[predictors;temp_predictors];
end
% end
%%
predictors    = cat(3,predictors{:});
targets       = cat(2,targets{:});

noisyMean     = mean(predictors(:));
noisyStd      = std(predictors(:));
predictors(:) = (predictors(:) - noisyMean)/noisyStd;


cleanMean     = mean(targets(:));
cleanStd      = std(targets(:));
targets(:)    = (targets(:) - cleanMean)/cleanStd;


predictors = reshape(predictors,size(predictors,1),size(predictors,2),1,size(predictors,3));
targets    = reshape(targets,1,1,size(targets,1),size(targets,2));

inds                = randperm(size(predictors,4));
L                   = round(.90* size(predictors,4));
trainPredictors     = predictors(:,:,:,inds(1:L));
trainTargets        = targets(:,:,:,inds(1:L));
validatePredictors  = predictors(:,:,:,inds(L+1:end));
validateTargets     = targets(:,:,:,inds(L+1:end));

%%  Denoising with Convolutional Layers
%%{
layers2 = [imageInputLayer([numFeatures,numSegments],'Name','ip')

            convolution2dLayer(3,8,"Stride",1,"Padding","same",'Name','Cl')
            leakyReluLayer('Name','r1')
            batchNormalizationLayer('Name','B1')
            
            
            convolution2dLayer(3,16,"Stride",2,"Padding","same",'Name','C2')
            leakyReluLayer('Name','r2')
            batchNormalizationLayer('Name','B2')
            
            
            convolution2dLayer(3,16,"Stride",1,"Padding","same",'Name','C3')
            leakyReluLayer('Name','r3')
            batchNormalizationLayer('Name','B3')
            
            
            convolution2dLayer(3,32,"Stride",2,"Padding","same",'Name','C4')
            leakyReluLayer('Name','r4')
            batchNormalizationLayer('Name','B4')
            
            
            convolution2dLayer(3,32,"Stride",1,"Padding","same",'Name','C5')
            leakyReluLayer('Name','r5')
            batchNormalizationLayer('Name','B5')
            
            
            convolution2dLayer(3,64,"Stride",2,"Padding","same",'Name','C6')
            leakyReluLayer('Name','r6')
            batchNormalizationLayer('Name','B6')
            
            
            convolution2dLayer(3,64,"Stride",1,"Padding","same",'Name','C7')
            leakyReluLayer('Name','r7')
            batchNormalizationLayer('Name','B7')
            
            
            convolution2dLayer(3,128,"Stride",2,"Padding","same",'Name','C8')
            leakyReluLayer('Name','r8')
            batchNormalizationLayer('Name','B8')
            
            
            convolution2dLayer(3,128,"Stride",1,"Padding","same",'Name','C8a')
            leakyReluLayer('Name','r8a')
            batchNormalizationLayer('Name','B8a')
            
            %%%
            convolution2dLayer(3,256,"Stride",2,"Padding","same",'Name','C8b')
            leakyReluLayer('Name','r8b')
            batchNormalizationLayer('Name','B8b')
            
            
            convolution2dLayer(3,256,"Stride",1,"Padding","same",'Name','C8c')
            leakyReluLayer('Name','r8c')
            batchNormalizationLayer('Name','B8c')
            
            
            convolution2dLayer(3,512,"Stride",2,"Padding","same",'Name','C8d')
            leakyReluLayer('Name','r8d')
            batchNormalizationLayer('Name','B8d')
            
            
            transposedConv2dLayer(3,512,'Stride',1,"Cropping","same",'Name','Tla')
            leakyReluLayer('Name','r9a')
            batchNormalizationLayer('Name','B9a')
            
            
            transposedConv2dLayer(3,256,'Stride',[2 1],"Cropping","same",'Name','T2a')
            leakyReluLayer('Name','rl0a')
            additionLayer(2,'Name','add6')
            batchNormalizationLayer('Name','Bl0a')
            
            
            
            
            transposedConv2dLayer(3,256,'Stride',1,"Cropping","same",'Name','T2b')
            leakyReluLayer('Name','rl0b')
            batchNormalizationLayer('Name','Bl0b')
            
            
            
            transposedConv2dLayer(3,128,'Stride',[2 1],"Cropping","same",'Name','Tlc')
            leakyReluLayer('Name','r9c')
            additionLayer(2,'Name','add5')
            batchNormalizationLayer('Name','B9c')
            
            
            
            %%%
            transposedConv2dLayer(3,128,'Stride',1,"Cropping","same",'Name','Tl')
            leakyReluLayer('Name','r9')
            batchNormalizationLayer('Name','B9')
            
            
            
            transposedConv2dLayer(3,64,'Stride',2,"Cropping","same",'Name','T2')
            leakyReluLayer('Name','rl0')
            additionLayer(2,'Name','add4')
            batchNormalizationLayer('Name','Bl0')
            
            
            
            %
            transposedConv2dLayer(3,64,'Stride',1,"Cropping","same",'Name','T3')
            leakyReluLayer('Name','rl1')
            batchNormalizationLayer('Name','Bl1')
            
            
            transposedConv2dLayer(3,32,'Stride',2,"Cropping","same",'Name','T4')
            leakyReluLayer('Name','rl2')
            additionLayer(2,'Name','add3')
            batchNormalizationLayer('Name','Bl2')
            
            
            
            transposedConv2dLayer(3,32,'Stride',1,"Cropping","same",'Name','T5')
            leakyReluLayer('Name','rl3')
            batchNormalizationLayer('Name','Bl3')
            
            
            transposedConv2dLayer(3,16,'Stride',2,"Cropping","same",'Name','T6')
            leakyReluLayer('Name','rl4')
            additionLayer(2,'Name','add2')
            batchNormalizationLayer('Name','Bl4')
            
            
            
            
            transposedConv2dLayer(3,16,'Stride',1,"Cropping","same",'Name','T7')
            leakyReluLayer('Name','rl5')
            additionLayer(2,'Name','add1')
            batchNormalizationLayer('Name','Bl5')
            
            
            
            
            transposedConv2dLayer(3,8,'Stride',2,"Cropping","same",'Name','T8')
            leakyReluLayer('Name','rl6')
            batchNormalizationLayer('Name','Bl6')
            
            
            convolution2dLayer(1,1,"Stride",[1 1],"Padding","same",'Name','C')
            fullyConnectedLayer(numFeatures,'Name','F')
            regressionLayer('Name','op')
            ];
            layers2 = layerGraph(layers2);
            layers2  = connectLayers(layers2 ,'r2','add1/in2');
            layers2  = connectLayers(layers2 ,'r3','add2/in2');
            layers2  = connectLayers(layers2,'r5','add3/in2');
            layers2  = connectLayers(layers2 ,'r7','add4/in2');
            layers2  = connectLayers(layers2,'r8a','add5/in2');
            layers2  = connectLayers(layers2 ,'r8c','add6/in2');

miniBatchSize = 512;
LearnRate = 1e-2;
options2 = trainingOptions("adam", ...
    "MaxEpochs",50, ...
    'ValidationPatience',5,...
    "InitialLearnRate",LearnRate,...
    'Plots','training-progress',...
    "MiniBatchSize",miniBatchSize, ...
    "Shuffle","once", ...
    "Verbose",0, ...
    'ExecutionEnvironment','auto',...
    "VerboseFrequency",floor(size(trainPredictors,4)/miniBatchSize), ...
    "ValidationFrequency",floor(size(trainPredictors,4)/miniBatchSize),...
    "LearnRateSchedule","piecewise",...
    "LearnRateDropFactor",0.9,...
    "LearnRateDropPeriod",1,...
    "ValidationData",{validatePredictors,permute(validateTargets,[1 2 3 4])});


denoiseNetDCT = trainNetwork(trainPredictors,permute(trainTargets,[1 2 3 4]),layers2,options2);
%% Test the Denoising Networks
load testing_data
y =R+z;
for ii=1:size(y,2)
    tr=ii;

    noisysignal   = y(:,tr);
    noisySTFT  = stft_dct(noisysignal,'Window',win,'OverlapLength',overlap,'FFTLength',ffTLength,'Center',false);
    TEMP = (noisySTFT);
    noisySTFT  = (noisySTFT(1:numFeatures,:));
    noisySTFT  = [zeros(numFeatures,ceil(numSegments/2-1)) noisySTFT zeros(numFeatures,ceil(numSegments/2-1))];
    predictors = zeros( numFeatures, numSegments , size(noisySTFT,2) - numSegments + 1);
    for index = 1 : size(noisySTFT,2) - numSegments + 1
        predictors(:,:,index) = noisySTFT(:,index:index + numSegments - 1);
    end
    predictors(:) =(predictors(:) -mean(predictors(:) ))/max(abs(predictors(:)));
    predictors(:) = (predictors(:) - noisyMean) / noisyStd;

    predictors = reshape(predictors, [numFeatures,numSegments,1,size(predictors,3)]);

    STFTFullyConvolutional = predict(denoiseNetDCT, predictors).';
    denoised = real(istft_dct(squeeze(STFTFullyConvolutional),  ...
        'Window',win,'OverlapLength',overlap, ...
        'FFTLength',ffTLength,'ConjugateSymmetric',false,'Center',false));

    GT(:,ii)=     R(1:numel(denoised),ii);

    DE(:,ii)=    real(denoised);


end
SNR_after_denoising=mean(10*log10(var(GT)./var(DE-GT)))
%%

figure
II=1;
ax(1) =subplot(331);
plot((GT(:,II))/max(abs(GT(:,II))));hold on
title('(a)')
ylabel({'Normalized'; 'Amplitude'})
ax(2) =subplot(332);
plot(y(:,II)/max(abs(y(:,II))));hold on
title('(b)')
ax(8) =subplot(333);
plot(DE(:,II)/max(abs(DE(:,II))));hold on
title('(c)')
linkaxes(ax);
II=2;
ax(3) =subplot(334);
plot((GT(:,II))/max(abs(GT(:,II))));hold on
title('(d)')
ylabel({'Normalized'; 'Amplitude'})
ax(4) =subplot(335);
plot(y(:,II)/max(abs(y(:,II))));hold on
title('(e)')
ax(9) =subplot(336);
plot(DE(:,II)/max(abs(DE(:,II))));hold on
title('(f)')
linkaxes(ax);
II=3;
ax(5) =subplot(337);
plot((GT(:,II))/max(abs(GT(:,II))));hold on
title('(g)')
ylabel({'Normalized'; 'Amplitude'})
xlabel('Sample no.')
ax(6) =subplot(338);
plot(y(:,II)/max(abs(y(:,II))));hold on
title('(h)')
xlabel('Sample no.')
ax(7) =subplot(339);
plot(DE(:,II)/max(abs(DE(:,II))));hold on
title('(i)')
linkaxes(ax);
xlabel('Sample no.')
bx=gca;
set(bx,'FontSize',14,'fontname','times');
set(gcf, 'Position',  [100, 100, 1000, 600])