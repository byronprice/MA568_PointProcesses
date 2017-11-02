% SimpleSpikingNetworkModel.m

%  Poisson GLM network model with Hebbian and STDP-like plasticity
%   illustrative example for talk in December
%   basically a leaky integrate-and-fire neuron network with plasticity
%   rules
%   simulates a Poisson process with history dependence and dependence on 
%    the other neuron's recent spiking history

simTime = 20;
sampleFreq = 1000;
timeMultiplier = 1000;
excitProb = 0.8;inhibProb = 1-excitProb;
histLen = 30;

x = 0:(histLen);y = sin(2*pi*x/histLen+3*pi/2);
xx = 0:floor(histLen/2);yy = sin(2*pi*xx/(2*histLen)+pi/2);
y(floor(histLen/2)+1:end) = yy;y(y<0) = y(y<0)*4;
subZeroSum = sum(abs(y(y<0)));y(y>0) = (y(y>0)./sum(y(y>0)))*subZeroSum;
excitHistFilter = y';

tempHistLen = 20;
x = 0:(tempHistLen);y = sin(2*pi*x/tempHistLen+3*pi/2);
xx = 0:floor(tempHistLen/2);yy = sin(2*pi*xx/(2*tempHistLen)+pi/2);
y(floor(tempHistLen/2)+1:end) = yy;y(y<0) = y(y<0)*3;
subZeroSum = sum(abs(y(y<0)));y(y>0) = (y(y>0)./sum(y(y>0)))*subZeroSum;
inhibHistFilter = y';

inhibHistFilter = [inhibHistFilter;zeros(histLen-tempHistLen,1)];

histLen = histLen+1;

% post-pre spike timing subtraction (pre-before-post is greater than zero)
plasticityLen = 100;
positiveTime = 1:plasticityLen;negativeTime = -plasticityLen:0;
time = [negativeTime,positiveTime];
negTimeLen = length(negativeTime);posTimeLen = length(positiveTime);
hebbTau = 12;
hebbianFilter = [0;exp(-(positiveTime(1:end-1)-1)./hebbTau)'];

stdpTau = 10;
stdpFilterExcitToExcit = [-exp((negativeTime)./stdpTau),1.1*exp(-(positiveTime-1)./stdpTau)]';
stdpFilterExcitToInhib = [exp((negativeTime)./stdpTau),-1.1*exp(-(positiveTime-1)./stdpTau)]';

x = 0:(plasticityLen*2);
halfWay = plasticityLen/2;
gaussKernel = exp(-((x-halfWay).*(x-halfWay))./(2*35*35));
stdpFilterInhibToExcit = (cos(2*pi*(x-halfWay)/125).*gaussKernel)';

connTau = 2;
inhibConnFilter = -exp(-(positiveTime-1)./connTau)';
excitConnFilter = 0.5*exp(-(positiveTime-1)./connTau)';


numNeurons = 50;spatialFootprintSize = [720,1280];
A = zeros(prod(spatialFootprintSize),numNeurons);
videoFrameRate = 60;
% temporalComponent = zeros(numNeurons,simTime*videoFrameRate);

for ii=1:numNeurons
    centerPos = [(spatialFootprintSize(1)-100)*rand+100,(spatialFootprintSize(2)-100)*rand+100];
    tempFootprint = zeros(spatialFootprintSize);
    for jj=1:spatialFootprintSize(1)
        for kk=1:spatialFootprintSize(2)
            dist = sqrt((centerPos(1)-jj).^2+(centerPos(2)-kk).^2);
            if dist < 30
               tempFootprint(jj,kk) = 1;
            end
        end
    end
    A(:,ii) = tempFootprint(:);
end
% to get back to the footprint
%  reshape(A(:,ii),[spatialFootprintSize(1),spatialFootprintSize(2)])';

% simple Hebbian rule, random connectivity to start
% baseRateExcit = 2/timeMultiplier;baseRateInhib = 10/timeMultiplier;
% pointProcessSpikes = zeros(numNeurons,simTime*timeMultiplier);
% pointProcessRate = zeros(numNeurons,simTime*timeMultiplier);
% excitNeurons = sort(randperm(numNeurons,round(0.8*numNeurons)));
% inhibNeurons = find(~ismember(1:numNeurons,excitNeurons));
% numInhibs = length(inhibNeurons);numExcits = length(excitNeurons);
% 
% Adjacency = zeros(simTime*timeMultiplier,4,4);
% 
% 
% currentSpikeHist = zeros(numNeurons,plasticityLen);
% for ii=1:numNeurons
%     if ismember(ii,excitNeurons) == 1
%         currentSpikeHist(ii,:) = poissrnd(baseRateExcit,[1,plasticityLen]);
%     else
%         currentSpikeHist(ii,:) = poissrnd(baseRateInhib,[1,plasticityLen]);
%     end
% end
% 
% inhibToExcitConnMat = rand([numInhibs,numExcits])*0.1;
% ExcitToExcitConnMat = rand([numExcits,numExcits])*0.1;
% ExcitToInhibConnMat = rand([numExcits,numInhibs])*0.1;
% inhibToInhibConnMat = rand([numInhibs,numInhibs])*0.1;
% 
% learningRate = 0.1;
% for ii=1:simTime*timeMultiplier
%     excitConnResults = currentSpikeHist(excitNeurons,:)*excitConnFilter;
%     inhibConnResults = currentSpikeHist(inhibNeurons,:)*inhibConnFilter;
%     
%     count = 1;
%     for jj=excitNeurons
%         inhibWeights = inhibToExcitConnMat(:,count);
%         excitWeights = ExcitToExcitConnMat(:,count);
%         
%         innerSum = sum(inhibWeights.*inhibConnResults)+...
%             sum(excitWeights.*excitConnResults)+...
%             currentSpikeHist(jj,1:histLen)*excitHistFilter;
% %          add correlated impulse
% %         if ii == timeMultiplier*10
% %             innerSum = innerSum+2;
% %         end
%         pointProcessRate(jj,ii) = exp(innerSum)*baseRateExcit;
%         pointProcessSpikes(jj,ii) = poissrnd(pointProcessRate(jj,ii));
%         currentSpikeHist(jj,:) = [pointProcessSpikes(jj,ii),currentSpikeHist(jj,1:end-1)];
%         count = count+1;
%     end
%     count = 1;
%     for jj=inhibNeurons
%         inhibWeights = inhibToInhibConnMat(:,count);
%         excitWeights = ExcitToInhibConnMat(:,count);
%         
%         innerSum = sum(inhibWeights.*inhibConnResults)+...
%             sum(excitWeights.*excitConnResults)+...
%             currentSpikeHist(jj,1:histLen)*inhibHistFilter;
% %         add correlated impulse
% %         if ii==timeMultiplier*10
% %             innerSum = innerSum-2;
% %         end
%         pointProcessRate(jj,ii) = exp(innerSum)*baseRateInhib;
%         pointProcessSpikes(jj,ii) = poissrnd(pointProcessRate(jj,ii));
%         currentSpikeHist(jj,:) = [pointProcessSpikes(jj,ii),currentSpikeHist(jj,1:end-1)];
%         count = count+1;
%     end
%     
%     changeConn = currentSpikeHist*hebbianFilter;
%     
%     count = 1;
%     for jj=excitNeurons
%         if pointProcessSpikes(jj,ii) == 1
%             tempChangeConn = changeConn;tempChangeConn(jj) = 0;
%             inhibToExcitConnMat(:,count) = inhibToExcitConnMat(:,count)+...
%                 tempChangeConn(inhibNeurons)*learningRate;
%             ExcitToExcitConnMat(:,count) = ExcitToExcitConnMat(:,count)+...
%                 tempChangeConn(excitNeurons)*learningRate;
%         end
%         count = count+1;
%     end
%     
%     count = 1;
%     for jj=inhibNeurons
%         if pointProcessSpikes(jj,ii) == 1
%             tempChangeConn = changeConn;tempChangeConn(jj) = 0;
% %             inhibToInhibConnMat(:,count) = inhibToInhibConnMat(:,count)+...
% %                 tempChangeConn(inhibNeurons)*learningRate;
%             ExcitToInhibConnMat(:,count) = ExcitToInhibConnMat(:,count)+...
%                 tempChangeConn(excitNeurons)*learningRate;
%         end
%         count = count+1;
%     end
%     
%       % excite 1, excite 2, excite 3, inhib 1
%     Adjacency(ii,1,2) = ExcitToExcitConnMat(1,2);
%     Adjacency(ii,2,1) = ExcitToExcitConnMat(2,1);
%     Adjacency(ii,1,3) = ExcitToExcitConnMat(1,3);
%     Adjacency(ii,3,1) = ExcitToExcitConnMat(3,1);
%     Adjacency(ii,1,4) = ExcitToInhibConnMat(1,1);
%     Adjacency(ii,4,1) = inhibToExcitConnMat(1,1);
%     Adjacency(ii,2,3) = ExcitToExcitConnMat(2,3);
%     Adjacency(ii,3,2) = ExcitToExcitConnMat(3,2);
%     Adjacency(ii,2,4) = ExcitToInhibConnMat(2,1);
%     Adjacency(ii,4,2) = inhibToExcitConnMat(1,2);
%     Adjacency(ii,3,4) = ExcitToInhibConnMat(3,1);
%     Adjacency(ii,4,3) = inhibToExcitConnMat(1,3);
% end
% 
% downSampleFactor = round(sampleFreq/50);
% tic;
% for ii=1:downSampleFactor:size(pointProcessRate,2)
%     frame = A*pointProcessRate(:,ii)*timeMultiplier;
%     frame = reshape(frame,[spatialFootprintSize(1),spatialFootprintSize(2)]);
%     imagesc(frame);colormap('bone');caxis([0 25]);
%     pause(1/50);
% end
% toc;


%% STDP rule, random connectivity to start
baseRateExcit = 2/timeMultiplier;baseRateInhib = 10/timeMultiplier;
pointProcessSpikes = zeros(numNeurons,simTime*timeMultiplier);
pointProcessRate = zeros(numNeurons,simTime*timeMultiplier);
excitNeurons = sort(randperm(numNeurons,round(0.8*numNeurons)));
inhibNeurons = find(~ismember(1:numNeurons,excitNeurons));
numInhibs = length(inhibNeurons);numExcits = length(excitNeurons);

Adjacency = zeros(simTime*timeMultiplier,4,4);


currentSpikeHist = zeros(numNeurons,plasticityLen*2+1);
for ii=1:numNeurons
    if ismember(ii,excitNeurons) == 1
        currentSpikeHist(ii,:) = poissrnd(baseRateExcit,[1,plasticityLen*2+1]);
    else
        currentSpikeHist(ii,:) = poissrnd(baseRateInhib,[1,plasticityLen*2+1]);
    end
end

inhibToExcitConnMat = rand([numInhibs,numExcits])/numInhibs;
ExcitToExcitConnMat = rand([numExcits,numExcits])/numExcits;
ExcitToExcitConnMat(find(diag(ones(numExcits,1)))) = 0;
ExcitToInhibConnMat = rand([numExcits,numInhibs])/numExcits;
inhibToInhibConnMat = rand([numInhibs,numInhibs])/numInhibs;
inhibToInhibConnMat(find(diag(ones(numInhibs,1)))) = 0;

learningRate = 0.1;
for ii=1:simTime*timeMultiplier
    excitConnResults = currentSpikeHist(excitNeurons,1:plasticityLen)*excitConnFilter;
    inhibConnResults = currentSpikeHist(inhibNeurons,1:plasticityLen)*inhibConnFilter;
    
    count = 1;
    for jj=excitNeurons
        inhibWeights = inhibToExcitConnMat(:,count);
        excitWeights = ExcitToExcitConnMat(:,count);
        
        innerSum = sum(inhibWeights.*inhibConnResults)+...
            sum(excitWeights.*excitConnResults)+...
            currentSpikeHist(jj,1:histLen)*excitHistFilter;
%          add correlated impulse
%         if ii == timeMultiplier*10
%             innerSum = innerSum+2;
%         end
        pointProcessRate(jj,ii) = exp(innerSum)*baseRateExcit;
        pointProcessSpikes(jj,ii) = poissrnd(pointProcessRate(jj,ii));
        currentSpikeHist(jj,:) = [pointProcessSpikes(jj,ii),currentSpikeHist(jj,1:end-1)];
        count = count+1;
    end
    count = 1;
    for jj=inhibNeurons
        inhibWeights = inhibToInhibConnMat(:,count);
        excitWeights = ExcitToInhibConnMat(:,count);
        
        innerSum = sum(inhibWeights.*inhibConnResults)+...
            sum(excitWeights.*excitConnResults)+...
            currentSpikeHist(jj,1:histLen)*inhibHistFilter;
%         add correlated impulse
%         if ii==timeMultiplier*10
%             innerSum = innerSum-2;
%         end
        pointProcessRate(jj,ii) = exp(innerSum)*baseRateInhib;
        pointProcessSpikes(jj,ii) = poissrnd(pointProcessRate(jj,ii));
        currentSpikeHist(jj,:) = [pointProcessSpikes(jj,ii),currentSpikeHist(jj,1:end-1)];
        count = count+1;
    end
    
    changeConnEE = currentSpikeHist*stdpFilterExcitToExcit;
    changeConnEI = currentSpikeHist*stdpFilterExcitToInhib;
    changeConnIE = currentSpikeHist*stdpFilterInhibToExcit;
    count = 1;
    for jj=excitNeurons
        if currentSpikeHist(jj,plasticityLen+1)==1
            tempChangeConn = changeConnIE;
            inhibToExcitConnMat(:,count) = max(inhibToExcitConnMat(:,count)+...
                tempChangeConn(inhibNeurons)*learningRate,0);
            
            tempChangeConn = changeConnEE;tempChangeConn(jj) = 0;
            ExcitToExcitConnMat(:,count) = max(ExcitToExcitConnMat(:,count)+...
                tempChangeConn(excitNeurons)*learningRate,0);
        end
        count = count+1;
    end
    
    count = 1;
    for jj=inhibNeurons
        if currentSpikeHist(jj,plasticityLen+1) == 1
            tempChangeConn = changeConnEI;
%             inhibToInhibConnMat(:,count) = max(inhibToInhibConnMat(:,count)+...
%                 tempChangeConn(inhibNeurons)*learningRate,0);
            ExcitToInhibConnMat(:,count) = max(ExcitToInhibConnMat(:,count)+...
                tempChangeConn(excitNeurons)*learningRate,0);
        end
        count = count+1;
    end
    % excite 1, excite 2, excite 3, inhib 1
    Adjacency(ii,1,2) = ExcitToExcitConnMat(1,2);
    Adjacency(ii,2,1) = ExcitToExcitConnMat(2,1);
    Adjacency(ii,1,3) = ExcitToExcitConnMat(1,3);
    Adjacency(ii,3,1) = ExcitToExcitConnMat(3,1);
    Adjacency(ii,1,4) = ExcitToInhibConnMat(1,1);
    Adjacency(ii,4,1) = inhibToExcitConnMat(1,1);
    Adjacency(ii,2,3) = ExcitToExcitConnMat(2,3);
    Adjacency(ii,3,2) = ExcitToExcitConnMat(3,2);
    Adjacency(ii,2,4) = ExcitToInhibConnMat(2,1);
    Adjacency(ii,4,2) = inhibToExcitConnMat(1,2);
    Adjacency(ii,3,4) = ExcitToInhibConnMat(3,1);
    Adjacency(ii,4,3) = inhibToExcitConnMat(1,3);
end

downSampleFactor = round(sampleFreq/50);
newProcessRate = zeros(numNeurons,length(pointProcessRate)/downSampleFactor);
for ii=1:numNeurons
   temp = pointProcessRate(ii,:);
   newProcessRate(ii,:) = decimate(pointProcessRate(ii,:),downSampleFactor);
end
    
tic;
for ii=1:size(newProcessRate,2)
    frame = A*newProcessRate(:,ii)*timeMultiplier;
    frame = reshape(frame,[spatialFootprintSize(1),spatialFootprintSize(2)]);
    imagesc(frame);colormap('bone');caxis([0 25]);
    pause(1/50);
end
toc;
