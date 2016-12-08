function [Results] = Byron_FitGLM(Process,X,testInds,modelName)
%Byron_FitGLM
%   Detailed explanation goes here
myGreen = [0 0.5 0];
myRed = [0.5 0 0];
myBlue = [0 0 0.5];

forLegend = {'Celtics','Opponent'};
forTitle = {'Shot Attempts','Hit Shots','Fouls'};

numCoeffs = size(X{1},2)+1;
numProcesses = length(testInds);
Results = zeros(numProcesses,numCoeffs,4);
for ii=1:numProcesses
    [b,dev,stats] = glmfit(X{testInds(8)},Process{testInds(ii)}','poisson');
    for jj=1:numCoeffs
        Results(ii,jj,:) = [b(jj),stats.se(jj),stats.p(jj),dev];
    end
end

minB = min(min(exp(Results(:,:,1))));
maxB = max(max(exp(Results(:,:,1))));
expError = Results(:,:,2).*exp(Results(:,:,1));

signifResults = double(Results(:,:,3)<0.05);
signifResults(signifResults>0) = -0.5;
signifResults(signifResults==0) = -50;
figure();
for ii=1:length(testInds)
   if ii <= 3
%         subplot(1,3,mod(ii,3)+1);e = errorbar(0:(numCoeffs-1),squeeze(exp(Results(testInds(ii),:,1))),...
%             2.*squeeze(expError(testInds(ii),:)),'+','LineWidth',2);axis([-1 numCoeffs minB-1 maxB+1]);
%         e.Color = myBlue;hold on;e.LineStyle = 'none';
        %plot(0:(numCoeffs-1),signifResults(testInds(ii),:),'*','Color',myBlue);
   elseif ii <= 6
       subplot(1,3,ii-3);e = errorbar(0:(numCoeffs-1),squeeze(exp(Results(ii,:,1))),...
            2.*squeeze(expError(ii,:)),'s','LineWidth',2);axis([-1 numCoeffs minB-1 maxB+1]);
        e.Color = myGreen;hold on;e.LineStyle = 'none';
       % plot(0:(numCoeffs-1),signifResults(testInds(ii),:),'*','Color',myGreen);
   elseif ii <= 9
       subplot(1,3,ii-6);e = errorbar(0:(numCoeffs-1),squeeze(exp(Results(ii,:,1))),...
            2.*squeeze(expError(ii,:)),'d','LineWidth',2);axis([-1 numCoeffs minB-1 maxB+1]);
        e.Color = myRed;hold on;e.LineStyle = 'none';
       % plot(0:(numCoeffs-1),signifResults(testInds(ii),:),'*','Color',myRed);
   end
end
subplot(1,3,1);legend(forLegend);title(['Coefficients for ',modelName,': ',forTitle{1}]);
 xlabel('Coefficient #');ylabel('Modulation Factor [exp(Coefficient)]');
subplot(1,3,2);legend(forLegend);title(['Coefficients for ',modelName,': ',forTitle{2}]);
xlabel('Coefficient #');ylabel('Modulation Factor [exp(Coefficient)]');
subplot(1,3,3);legend(forLegend);title(['Coefficients for ',modelName,': ',forTitle{3}]);
xlabel('Coefficient #');ylabel('Modulation Factor [exp(Coefficient)]');

for ii=1:length(testInds)
   if ii <= 6 && ii > 3
       for jj=0:(numCoeffs-1)
           if jj==0
               subplot(1,3,ii-3);plot(jj,signifResults(ii,jj+1),'*','Color',myGreen);
           elseif Results(ii,jj+1,1) > 0 
                subplot(1,3,ii-3);plot(jj,signifResults(ii,jj+1),'^','Color',myGreen);
           elseif Results(ii,jj+1,1) < 0 
                subplot(1,3,ii-3);plot(jj,signifResults(ii,jj+1),'v','Color',myGreen);
           end
       end
       set(gca,'FontSize',12);
   elseif ii <= 9 && ii > 6
       for jj=0:(numCoeffs-1)
           if jj==0
               subplot(1,3,ii-6);plot(jj,signifResults(ii,jj+1)-0.1,'*','Color',myRed);
           elseif Results(ii,jj+1,1) > 0 
                subplot(1,3,ii-6);plot(jj,signifResults(ii,jj+1)-0.1,'^','Color',myRed);
           elseif Results(ii,jj+1,1) < 0
                subplot(1,3,ii-6);plot(jj,signifResults(ii,jj+1)-0.1,'v','Color',myRed);
           end
       end
       set(gca,'FontSize',12);
   end
end
end

