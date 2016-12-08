% Uri_Project.m
myGreen = [0 0.5 0];
myRed = [0.5 0 0];
myBlue = [0 0 0.5];
load('NBASampleData.mat');
totalTime = length(T);

quarterTimes = zeros(4,2);
for ii=1:4
    quarterTimes(ii,1) = 1+60*12*(ii-1);
    quarterTimes(ii,2) = 60*12*ii;
end
halfTimes = zeros(2,2);
for ii=1:2
   halfTimes(ii,1) = 1+60*24*(ii-1);
   halfTimes(ii,2) = 60*24*ii;
end
gameTime = max(T);

% remove outliers ... three games had > 10 minutes with no shots taken
indeces = [1:115200-1,115200+gameTime:152641-1,152641+gameTime:184321-1,...
    184321+gameTime:totalTime];

T = T(indeces);
shots = shots(:,indeces);oshots = oshots(:,indeces);
fouls = fouls(indeces);ofouls = ofouls(indeces);
totalTime = length(T);numGames = totalTime/gameTime;

%  create cell array with each point process inside
%Order goes all shots, all hit shots, all missed shots, all 2's, all 3's,
%  all hit 2's, all hit 3's, repeat same order for celtics, repeat same
%  order for opponents (that's 21), add all fouls, celtic fouls, opponent
%  fouls for 24
    
% shot probability array
%   all hits, all 2's, all 3's, repeat for celts, repeat for opponents

ProcessNames = {'All Shot Attempts','All Hit Shots','All Missed Shots',...
    'All 2 Attempts','All Hit 2''s','All Miss 2''s','All 3 Attempts',...
    'All Hit 3''s','All Miss 3''s','Celtic Shot Attempts',...
    'Celtic Hit Shots','Celtic Missed Shots','Celtic 2 Attempts','Celtic Hit 2''s',...
    'Celtic Miss 2''s','Celtic 3 Attempts','Celtic Hit 3''s','Celtic Miss 3''s',...
    'Opponent Shot Attempts','Opponent Hit Shots','Opponent Missed Shots',...
    'Opponent 2 Attempts','Opponent Hit 2''s','Opponent Miss 2''s',...
    'Opponent 3 Attempts','Opponent Hit 3''s','Opponent Miss 3''s',...
    'All Fouls','Celtic Fouls','Opponent Fouls'};

Process = cell(30,1);
ShotProbs = cell(9,1);

% all shots
Process{1} = sum(shots+oshots,1);
Process{10} = sum(shots,1);Process{19} = sum(oshots,1);

% hit shots
inds = [1,3];
Process{2} = sum(shots(inds,:)+oshots(inds,:),1);
Process{11} = sum(shots(inds,:),1);
Process{20} = sum(oshots(inds,:),1);

% missed shots
inds = [2,4];
Process{3} = sum(shots(inds,:)+oshots(inds,:),1);
Process{12} = sum(shots(inds,:),1);
Process{21} = sum(oshots(inds,:),1);

% all 2s
inds = [1,2];
Process{4} = sum(shots(inds,:)+oshots(inds,:),1);
Process{13} = sum(shots(inds,:),1);
Process{22} = sum(oshots(inds,:),1);

% all 3s
inds = [3,4];
Process{7} = sum(shots(inds,:)+oshots(inds,:),1);
Process{16} = sum(shots(inds,:),1);
Process{25} = sum(oshots(inds,:),1);

% hit 2s
inds = 1;
Process{5} = sum(shots(inds,:)+oshots(inds,:),1);
Process{14} = shots(inds,:);Process{23} = oshots(inds,:);

% hit 3s
inds = 3;
Process{8} = sum(shots(inds,:)+oshots(inds,:),1);
Process{17} = shots(inds,:);Process{26} = oshots(inds,:);

% miss 2s
inds = 2;
Process{6} = sum(shots(inds,:)+oshots(inds,:),1);
Process{15} = shots(inds,:);Process{24} = oshots(inds,:);

% miss 3s
inds = 4;
Process{9} = sum(shots(inds,:)+oshots(inds,:),1);
Process{18} = shots(inds,:);Process{27} = oshots(inds,:);

Process{28} = fouls+ofouls;Process{29} = fouls;Process{30} = ofouls;
% for ii=1:30
%     Process{ii} = Process{ii}>0;
% end

figure();plot(T(1:gameTime)./60,Process{10}(1:gameTime),'Color',myGreen);hold on;
plot(T(1:gameTime)./60,Process{19}(1:gameTime),'Color',myRed);legend({'Celtics','Opponent'});
title('Shots Taken in One Game');xlabel('Game Time (Minutes)');
ylabel('Shot Taken (1) or Not Taken (0)');

% probability of making shots
for ii=1:3
    ShotProbs{1+3*(ii-1)} = sum(Process{2+9*(ii-1)})/sum(Process{1+9*(ii-1)});
    ShotProbs{2+3*(ii-1)} = sum(Process{5+9*(ii-1)})/sum(Process{4+9*(ii-1)});
    ShotProbs{3+3*(ii-1)} = sum(Process{8+9*(ii-1)})/sum(Process{7+9*(ii-1)});
end

figure();
hb = bar([ShotProbs{1},ShotProbs{4},ShotProbs{7};ShotProbs{2},...
    ShotProbs{5},ShotProbs{8};ShotProbs{3},ShotProbs{6},ShotProbs{9}]);hold on;
set(hb(1),'facecolor',myBlue);set(hb(2),'facecolor',myGreen);set(hb(3),'facecolor',myRed);
axis([0 4 0 2]);
set(gca,'XTick',[1,2,3]);
set(gca,'XTickLabel',{'Any Shot','Two','Three'});
legend({'Both','Celtics','Opponents'},'FontSize',20);


scatter([1-0.2,2-0.2,3-0.2],[(ShotProbs{2}*2*sum(Process{5})+...
    ShotProbs{3}*3*sum(Process{8}))/(sum(Process{5})+sum(Process{8})),...
    ShotProbs{2}*2,ShotProbs{3}*3],50,myBlue);
scatter([1,2,3],[(ShotProbs{5}*2*sum(Process{14})+...
    ShotProbs{6}*3*sum(Process{17}))/(sum(Process{14})+sum(Process{17})),...
    ShotProbs{5}*2,ShotProbs{6}*3],50,myGreen);
scatter([1+0.2,2+0.2,3+0.2],[(ShotProbs{8}*2*sum(Process{23})+...
    ShotProbs{9}*3*sum(Process{26}))/(sum(Process{23})+sum(Process{26})),...
    ShotProbs{8}*2,ShotProbs{9}*3],50,myRed);hold off;


ISI = cell(30,1);
for ii=1:30
    ISI{ii} = zeros(sum(Process{ii}),1);
    count = 1;
    isiCount = 0;
    for jj=1:totalTime
        if Process{ii}(jj) == 0
            isiCount = isiCount+1;
        else
            ISI{ii}(count) = isiCount;
            count = count+1;
            isiCount = 0;
        end
        if T(jj) == gameTime
            isiCount = 0;
        end
    end
end

figure();
% subplot(3,3,1);histogram(ISI{1}(:),'FaceColor',myBlue);
% title('Inter-Shot Interval For All Shot Attempts');xlabel('Inter-Shot Interval (seconds)');
% ylabel('Count');axis([0 500 0 1000]);
% 
% subplot(3,3,2);histogram(ISI{2}(:),'FaceColor',myBlue);
% title('Inter-Shot Interval For All Hit Shots');xlabel('Inter-Shot Interval (seconds)');
% ylabel('Count');axis([0 500 0 1000]);
% 
% subplot(3,3,3);histogram(ISI{28}(:),'FaceColor',myBlue);
% title('Inter-Foul Interval For All Fouls');xlabel('Inter-Foul Interval (seconds)');
% ylabel('Count');axis([0 500 0 1000]);

subplot(2,3,1);histogram(ISI{10}(:),'FaceColor',myGreen);
title('Inter-Shot Interval Histogram For Celtic Shot Attempts');xlabel('Inter-Shot Interval (seconds)');
ylabel('Count');axis([0 500 0 1000]);

subplot(2,3,2);histogram(ISI{11}(:),'FaceColor',myGreen);
title('Inter-Shot Interval Histogram For Celtic Hit Shots');xlabel('Inter-Shot Interval (seconds)');
ylabel('Count');axis([0 500 0 1000]);

subplot(2,3,3);histogram(ISI{29}(:),'FaceColor',myGreen);
title('Inter-Foul Interval Histogram For Celtic Fouls');xlabel('Inter-Foul Interval (seconds)');
ylabel('Count');axis([0 500 0 1000]);

subplot(2,3,4);histogram(ISI{19}(:),'FaceColor',myRed);
title('Inter-Shot Interval Histogram For Opponent Shot Attempts');xlabel('Inter-Shot Interval (seconds)');
ylabel('Count');axis([0 500 0 1000]);

subplot(2,3,5);histogram(ISI{20}(:),'FaceColor',myRed);
title('Inter-Shot Interval Histogram For Opponent Hit Shots');xlabel('Inter-Shot Interval (seconds)');
ylabel('Count');axis([0 500 0 1000]);

subplot(2,3,6);histogram(ISI{21}(:),'FaceColor',myRed);
title('Inter-Foul Interval Histogram For Opponent Fouls');xlabel('Inter-Foul Interval (seconds)');
ylabel('Count');axis([0 500 0 1000]);

% point process analysis ... dependence on time into game, on quarter, on
% half, and on time remaining in quarter
quarter = zeros(totalTime,4);
for ii=1:4
    indeces = find(T>=quarterTimes(ii,1) & T<=quarterTimes(ii,2));
    quarter(indeces,ii) = 1;
end

% quarter = zeros(totalTime,4);
% for ii=1:totalTime
%     for jj=1:4
%         if T(ii) >= quarterTimes(jj,1) && T(ii) <= quarterTimes(jj,2)
%            quarter(ii,jj) = 1; 
%         end
%     end
% end

half = zeros(totalTime,2);
for ii=1:2
    indeces = find(T>=halfTimes(ii,1) & T<=halfTimes(ii,2));
    half(indeces,ii) = 1;
end

quarterRemaining = zeros(totalTime,1);
for ii=1:totalTime
    if T(ii)==quarterTimes(1,1) || T(ii)==quarterTimes(2,1) || T(ii)==quarterTimes(3,1) || T(ii)==quarterTimes(4,1)
        value = 720;
    end
    quarterRemaining(ii) = value;
    value = value-1;
end

testInds = [1,2,28,10,11,29,19,20,30];

% display('Game-Clock Model');
% [time_Results] = Byron_FitGLM_Inhomogeneous(Process,T,testInds,'Game-Clock Model');
% 
% display('Which-Quarter Model: LAMBDA = exp(b0+b1*I(Q1)+b2*I(Q2)+b3*I(Q3)+b4*I(Q4))');
% [quarter_Results] = Byron_FitGLM_Inhomogeneous(Process,quarter,testInds,'Which-Quarter Model');
% 
% display('Which-Half Model');
% [half_Results] = Byron_FitGLM_Inhomogeneous(Process,half,testInds,'Which-Half Model');
% 
% display('Time-Remaining-In-Quarter Model');
% [remain_Results] = Byron_FitGLM_Inhomogeneous(Process,quarterRemaining,testInds,'Time-Remaining-In-Quarter Model');

%display('Full Model');
%[full_Results] = Byron_FitGLM_Inhomogeneous(Process,[T,quarter,quarterRemaining],testInds,'Full Inhomogeneous Model');

% history-dependent models
ProcessHistory = cell(30,1);
for ii=1:30
    ProcessHistory{ii} = zeros(totalTime,30);
    for jj=1:30
        tempProcess = Process{ii}';
        zeroVec = zeros(jj,1);
        tempProcess = [zeroVec;tempProcess];
        tempProcess = tempProcess(1:end-jj);
        ProcessHistory{ii}(:,jj) = tempProcess;
    end
end

% devianceNested = time_Results(testInds,1,4);devianceFull = full_Results(testInds,1,4);
% k = size(full_Results,2)-size(time_Results,2);
% diffDevs = devianceNested-devianceFull;
% p_values = 1-chi2cdf(diffDevs,k);

% display('History-Dependence Model');
% [history_Results] = Byron_FitGLM(Process,ProcessHistory,testInds,'History-Dependence Model');

% display('Celtic-Shot-Attempt-Dependence Model');
% [celticAttempt_Results] = Byron_FitGLM(Process,ProcessHistory,testInds,'Celtic-Shot-Attempt-Dependence Model');

% display('Celtic-Hit-Shot-Dependence Model');
% [celticHit_Results] = Byron_FitGLM(Process,ProcessHistory,testInds,'Celtic-Hit-Shot-Dependence Model');
% 
% display('Opponent-Shot-Attempt-Dependence Model');
% [oppAttempt_Results] = Byron_FitGLM(Process,ProcessHistory,testInds,'Opponent-Shot-Attempt-Dependence Model');
% 
% [oppHit_Results] = Byron_FitGLM(Process,ProcessHistory,testInds,'Opponent-Hit-Shot-Dependence Model');


% nested models for Celtics hit shots
AIC = zeros(7,1);modelParams = zeros(7,1);

focusProcess = Process{11}';

[b,dev,~] = glmfit(T,focusProcess,'poisson');
modelParams(1) = length(b);
AIC(1) = dev+2*modelParams(1);

missHistory = ProcessHistory{12};
[b,dev,~] = glmfit([T,missHistory],focusProcess,'poisson');
modelParams(2) = length(b);
AIC(2) = dev+2*modelParams(2);

hitHistory = ProcessHistory{11};
[b,dev,~] = glmfit([T,missHistory,hitHistory],focusProcess,'poisson');
modelParams(3) = length(b);
AIC(3) = dev+2*modelParams(3);

oppMissHist = ProcessHistory{21};
[b,dev,~] = glmfit([T,missHistory,hitHistory,oppMissHist],focusProcess,'poisson');
modelParams(4) = length(b);
AIC(4) = dev+2*modelParams(4);

oppHitHist = ProcessHistory{20};
[b,dev,~] = glmfit([T,missHistory,hitHistory,oppMissHist,oppHitHist],focusProcess,'poisson');
modelParams(5) = length(b);
AIC(5) = dev+2*modelParams(5);

selfFouls = ProcessHistory{29};
[b,dev,~] = glmfit([T,missHistory,hitHistory,oppMissHist,oppHitHist,selfFouls],focusProcess,'poisson');
modelParams(6) = length(b);
AIC(6) = dev+2*modelParams(6);

oppFouls = ProcessHistory{30};
[b,dev,stats] = glmfit([T,missHistory,hitHistory,oppMissHist,oppHitHist,selfFouls,oppFouls],focusProcess,'poisson');
modelParams(7) = length(b);
AIC(7) = dev+2*modelParams(7);

figure();plot(modelParams,AIC,'Color',myGreen,'LineWidth',2);title('AIC for Nested Models of Celtic-Hit-Shots Process');
xlabel('Number of Model Parameters');ylabel('AIC');set(gca,'FontSize',12);
set(gca,'XTicks',modelParams);
set(gca,'XTickLabel',{'Game Clock','Celtic Miss Hist',...
    'Celtic Hit Hist','Opp Miss Hist','Opp Hit Hist','Celtic Foul Hist','Opp Foul Hist'});

% two-sided KS test for shot attempts and made shots
display('Attempted Shots');
[h,p,ks2stat] = kstest2(ISI{10},ISI{19});
h
p
ks2stat

display('Hit Shots');
[h,p,ks2stat] = kstest2(ISI{11},ISI{20});
h
p
ks2stat

display('Fouls');
[h,p,ks2stat] = kstest2(ISI{29},ISI{30});
h
p
ks2stat

% time-rescaling theorem for goodness of fit
%[b,dev,stats] = glmfit([],focusProcess,'poisson');
errors = stats.se.*(exp(b).^2);
figure();hold on;
e = errorbar(0,exp(b(1)),2*errors(1),'*b','LineWidth',2);hold on;
e.LineStyle = 'none';
e = errorbar(1,exp(b(2)),2*errors(2),'*k','LineWidth',2);
e.LineStyle = 'none';
e = errorbar(2:31,exp(b(3:32)),2*errors(3:32),'*c','LineWidth',2);
e.LineStyle = 'none';
e = errorbar(32:61,exp(b(33:62)),2*errors(33:62),'*g','LineWidth',2);
e.LineStyle = 'none';
e = errorbar(62:91,exp(b(63:92)),2*errors(63:92),'*m','LineWidth',2);
e.LineStyle = 'none';
e = errorbar(92:121,exp(b(93:122)),2*errors(93:122),'*r','LineWidth',2);
e.LineStyle = 'none';
e = errorbar(122:151,exp(b(123:152)),2*errors(123:152),'*b','LineWidth',2);
e.LineStyle = 'none';
e = errorbar(152:181,exp(b(153:182)),2*errors(153:182),'*k','LineWidth',2);
e.LineStyle = 'none';
legend({'Baseline','Game Clock','Celtic Miss','Celtic Hit','Opp Miss','Opp Hit','Celtic Foul','Opp Foul'});
axis([-1 190 -1 3.5]);

pvals = stats.p;
signifResults = double(pvals<0.05);
signifResults(signifResults==1) = -0.5;
signifResults(signifResults== 0) = -50;

plot(0,signifResults(1),'vb');plot(1,signifResults(2),'vk');
for ii=2:31
    if b(ii+1) < 0
        plot(ii,signifResults(ii+1),'vc')
    elseif b(ii+1) > 0
        plot(ii,signifResults(ii+1),'^c')
    end
end  
for ii=32:61
    if b(ii+1) < 0
        plot(ii,signifResults(ii+1),'vg')
    elseif b(ii+1) > 0
        plot(ii,signifResults(ii+1),'^g')
    end
end 
for ii=62:91
    if b(ii+1) < 0
        plot(ii,signifResults(ii+1),'vm')
    elseif b(ii+1) > 0
        plot(ii,signifResults(ii+1),'^m')
    end
end 
for ii=92:121
    if b(ii+1) < 0
        plot(ii,signifResults(ii+1),'vr')
    elseif b(ii+1) > 0
        plot(ii,signifResults(ii+1),'^r')
    end
end 
for ii=122:151
    if b(ii+1) < 0
        plot(ii,signifResults(ii+1),'vb')
    elseif b(ii+1) > 0
        plot(ii,signifResults(ii+1),'^b')
    end
end 
for ii=152:181
    if b(ii+1) < 0
        plot(ii,signifResults(ii+1),'vk')
    elseif b(ii+1) > 0
        plot(ii,signifResults(ii+1),'^k')
    end
end 

% revist model with fewer coefficients, check AIC
[b,dev,stats] = glmfit([T,missHistory(:,1:15),hitHistory,oppMissHist(:,1:15)...
    ,oppHitHist,selfFouls(:,1:15),oppFouls(:,1:15)],focusProcess,'poisson');

lambdaT = zeros(totalTime,1);
lambdaT = exp(b(1)).*exp(b(2)*T).*exp((b(3:17)'*missHistory(:,1:15)')').*...
    exp((b(18:47)'*hitHistory')').*exp((b(48:62)'*oppMissHist(:,1:15)')').*...
    exp((b(63:92)'*oppHitHist')').*exp((b(93:107)'*selfFouls(:,1:15)')').*...
    exp((b(108:122)'*oppFouls(:,1:15)')');

numHit = sum(focusProcess);
Z = zeros(numHit,1);
indeces = find(focusProcess==1);

lowInd = 1;
for ii=1:numHit
    Z(ii) = sum(lambdaT(lowInd:indeces(ii)));
    lowInd = indeces(ii)+1;
end
figure();histogram(Z);title('Histogram of Rescaled Waiting Times');
ylabel('Count');xlabel('Rescaled Time');

sortedISI = sort(Z);
ISItimes = 0:0.01:max(sortedISI);
eCDF = zeros(length(ISItimes),1);
for ii=1:length(ISItimes)
    eCDF(ii) = sum(sortedISI <= ISItimes(ii));
end
eCDF = eCDF./numHit;

tCDF = 1-exp(-ISItimes');
bounds = (1.36/sqrt(numHit)).*ones(length(ISItimes),1);
figure();plot(eCDF,tCDF,'LineWidth',2);hold on;
xlabel('Empirical CDF');ylabel('Theoretical CDF');
title('KS Plot with 95% Confidence Bounds');
plot(tCDF,tCDF,'k');plot(tCDF,min(tCDF+bounds,1),'r');
plot(tCDF,max(tCDF-bounds,0),'r');hold off;

KSstat = max(abs(eCDF-tCDF));
display(sprintf('KS statistic: %3.3f',KSstat));