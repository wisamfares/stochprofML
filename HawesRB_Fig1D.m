clear
close all
set(groot, 'DefaultTextInterpreter', 'none');
set(groot, 'DefaultAxesTickLabelInterpreter', 'none');
set(groot, 'DefaultLegendInterpreter', 'none');

%% Set simulation parameters

mu1 = 0.5;          %First log mean
mu2 = -2.5;         %Second log mean
F = 0.225;          %Fraction in the high population
sigma = 0.225;      %Log standard deviation
mu_ps = 0.01;       %Reference log mean
mu_mb = 0.1;        %Mini-bulk log mean
mu_batch = 0.1;       %Batch log mean
sigma_ps = 0.005;     %Reference log standard deviation
sigma_batch = 0.01;  %Batch log standard deviation

%% Create synthetic mini-bulk datasets
rng(8);   %Set seed
numsamps = 2000; %Set number of observations
figure(1)
hold on
subplot(2,5,2);
plot(linspace(0.8,1.2,100),logndens(linspace(0.8, ...
    1.2,100),mu_ps,sigma_ps)/sum(logndens(linspace(0.8, ...
    1.2,100),mu_ps,sigma_ps)));
axis([0.8 1.2 0 0.4])
title('Technical effect')
text(0.82,0.35,sprintf('β = %0.2f\nγ = %0.3f',mu_ps,sigma_ps));

subplot(2,5,3);
plot([],[]);
box on;
line([1 1],[0 1]);
axis([0.8 1.2 0 0.4]);
title('Batch effect');
text(0.82,0.35,'ẟ(x)');

subplot(2,5,7);
plot(linspace(0.8,1.2,100),logndens(linspace(0.8, ...
    1.2,100),mu_ps,sigma_ps)/sum(logndens(linspace(0.8, ...
    1.2,100),mu_ps,sigma_ps)));
axis([0.8 1.2 0 0.4])
text(0.82,0.35,sprintf('β = %0.2f\nγ = %0.3f',mu_ps,sigma_ps));

subplot(2,5,8);
plot(linspace(0.5,1.5,100),logndens(linspace(0.5, ...
    1.5,100),mu_batch,sigma_batch)/sum(logndens(linspace(0.5, ...
    1.5,100),mu_batch,sigma_batch)));
axis([0.8 1.2 0 0.4]);
text(0.82,0.35,sprintf('β = %0.2f\nγ = %0.3f',mu_batch,sigma_batch));

numcells = 10;
mbdata1 = zeros(numsamps,1);
for i = 1:numsamps
    isamp = 0;
    num1 = binornd(numcells,F);
    for j = 1:num1
        isamp = isamp + lognrnd(mu1,sigma);
    end
    for j = num1+1:numcells
        isamp = isamp + lognrnd(mu2,sigma);
    end
    mbdata1(i) = isamp;
end
[y,x] = kde(mbdata1,'EvaluationPoints',linspace(0,20,101));
subplot(2,5,4);
plot(x,y/sum(y));
axis([-1 20 0 0.04]);
title('Mini-bulk + residual');

subplot(2,5,1);
%Convolve and normalize
mbdata1blend = conv(y,logndens(linspace(0.8, ...
    1.2,100),mu_ps,sigma_ps),'same');
mbdata1blend = mbdata1blend/sum(mbdata1blend);
plot(x,mbdata1blend);
axis([-1 20 0 0.04]);
ylabel('Probability density');
title('Original dataset');
text(12,0.035,{'Batch 1'});

mbdata2 = zeros(numsamps,1);
for i = 1:numsamps
    isamp = 0;
    num1 = binornd(numcells,F);
    for j = 1:num1
        isamp = isamp + lognrnd(mu1,sigma);
    end
    for j = num1+1:numcells
        isamp = isamp + lognrnd(mu2,sigma);
    end
    mbdata2(i) = isamp;
end
[y,x] = kde(mbdata2,'EvaluationPoints',linspace(0,20,101));
subplot(2,5,9);
plot(x,y/sum(y));
axis([-1 20 0 0.04]);

subplot(2,5,6);
%Double convolve and normalize
mbdata2blend = conv(conv(y,logndens(linspace(0.8,1.2,100),mu_ps,sigma_ps),'same'), ...
    logndens(linspace(0.8,1.2,100),mu_batch,sigma_batch),'same');
mbdata2blend = mbdata2blend/sum(mbdata2blend);
plot(x,mbdata2blend);
axis([-1 20 0 0.04]);
ylabel('Probability density')
text(12,0.035,{'Batch 2'});

[y,x] = kde([mbdata1; mbdata2],'EvaluationPoints',linspace(0,20,101));
subplot(2,5,10);
plot(x,y/sum(y));
axis([-1 20 0 0.04]);
title('Combined corrected')

%% Clean up graphs
for i = 1:10
    if i ~= 5
        subplot(2,5,i);
        h = gca;
        stylegraph(h);
        pbaspect([1 1 1]);
        xlabel('Target abundance');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d=logndens(x,m,s) 
d = (1./(x*s*sqrt(2*pi()))).*exp(-((log(x-m)).^2/(2*s^2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stylegraph(h)
%This function changes the graph style from MATLAB defaults
grid off %turn off grid
h.XColor = [0 0 0]; %change x-axis color
h.YColor = [0 0 0]; %change y-axis color
h.TickDir = 'out'; %tick marks out
h.TickLength = [0.02 0.05]; %increase tick length
h.Box = 'off'; %no box
end
