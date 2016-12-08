% Script to generate and fit simulated data according to a conditional intensity
% with spline based history dependence
clear; rng(0);

% Define spline parameters 
lastknot = 100; numknots = 5;
%c_pt_times_all = [-10 linspace(0,lastknot+1,numknots+1) lastknot+10];
c_pt_times_all = [-10 0 5 10 30 50 75 101 102];
s = 0.5;  % Define Tension Parameter

% Construct spline matrix
S = zeros(lastknot,numknots);
for i=1:lastknot  
    nearest_c_pt_index = max(find(c_pt_times_all<i));
    nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
    next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
    next2 = c_pt_times_all(nearest_c_pt_index+2);
    u = (i-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    l = (next2-next_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    p=[u^3 u^2 u 1]*[-s 2-s/l s-2 s/l;2*s s/l-3 3-2*s -s/l;-s 0 s 0;0 1 0 0];
    S(i,nearest_c_pt_index-1:nearest_c_pt_index+2) = p;   
end

% Simulate spiking activity 
nsteps = 50000;
spiketrain = zeros(nsteps,1);
theta = [log(.05) 0 -2 .5 0 -.2 .2 0 0 0];
for i=lastknot+1:nsteps,
    lambda(i) = exp(theta*[1; S'*spiketrain(i-1:-1:i-lastknot)]);
    spiketrain(i) = min(poissrnd(lambda(i)),1);
end;

%figure();plot(spiketrain);

% Build design matrix for multiplicative history model 
Hist = [];
for i=1:lastknot,
	Hist = [Hist spiketrain(lastknot-i+1:end-i)];
end;
figure();
for tt = 1:100
    plot(S(tt,:));
    pause(0.1)
end
X = Hist*S;
y=spiketrain(lastknot+1:end);

% Fit point process GLM 
[b dev stats] = glmfit(X,y,'poisson');
[yhat,dylo,dyhi] = glmval(b,S,'log',stats);
ytrue = glmval(theta',S,'log');

figure();
% Plot results
plot(1:lastknot,ytrue,1:lastknot,yhat,1:lastknot,yhat+dyhi,'r',1:lastknot,yhat-dylo,'r');
xlabel('Lag (ms)');
ylabel('Intensity based on a single spike at given lag');
legend('True','Model fit','Error bounds');