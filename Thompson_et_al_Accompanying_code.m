clc
clear all
close all

%% This code accompanies the manuscript entitled "Link between the number of virus
%% particles and variants founding new HIV-1 infections depends on the timing of transmission"
%% by Thompson et al. For further information about the paper or this code, please email
%% robin.thompson@chch.ox.ac.uk

%% We request that users cite the original publication when referring to this code or
%% any results generated from it.

%% The intention of this code is not to replicate every figure in the original publication,
%% but is instead to allow readers to experiment with the behaviour of model under different parameter
%% values. However, by changing the parameter values to match those at each point in the paper, a
%% number of the results in the publication above can be straightforwardly reproduced.


% Please note - this code may take some time to run

% First, we plot the distribution of variants in donors for different times since the
% donor became infected. For default parameters, these are the distributions shown in
% the left column of Fig. 3 of the publication.

maxnvariants = 56;
timesToShow = [0.4 2.2 4 6.4 8.4];
nTimes = length(timesToShow);

variantsByTime = zeros(maxnvariants, nTimes);

for tVal = 1:nTimes
    for str = 1:maxnvariants
        variantsByTime(str, tVal) = gampdf(str, 0.417, timesToShow(tVal)/0.563);
    end
    variantsByTime(:,tVal) = variantsByTime(:,tVal)./sum(variantsByTime(:,tVal));
end

for i = 1:length(timesToShow)
    subplot(3,2,i)
    bar([1:maxnvariants], variantsByTime(:,i), 'barWidth', 1)
    xlim([0.5 10.5])
    ylim([0 1])
    if i == 1
        title([strcat({'Time since infection:'}, {' '}, {num2str(timesToShow(i))}, {' '}, {'year'})])
    else
        title([strcat({'Time since infection:'}, {' '}, {num2str(timesToShow(i))}, {' '}, {'years'})])
    end
    ylabel('Proportion within host')
    xlabel('Variant (xth most common)')
end
figtitleRT('Proportion of each variant in donors (cf Fig 3 of our publication)','fontweight','bold');


% Now we characterise the number of particles available for transmission
% within donors in the population.

% First calculate the proportion of donors with each SPVL in the population

ncVals = [10^2 10^2.5 10^3 10^3.5 10^4 10^4.5 10^5 10^5.5 10^6 10^6.5 10^7]; % The different SPVLs of donors
kappa = 1;
alpha = -3.55;
sigma = 0.78/(sqrt(1 - ((2*alpha^2)/(pi*(1 + alpha^2)))));
mu = 4.74 - (2*sigma*alpha)/(sqrt(2*pi*(1 + alpha^2)));

for i = 1:length(ncVals)
    ncVals(i) = round(ncVals(i));
end

Dmax = 25.4;
Dfifty = 3058;
Dk = 0.41;

g = zeros(length(ncVals), 1);
tauc = zeros(length(ncVals), 1);

for j = 1:length(ncVals)
    nc = ncVals(j);
    g(j) = (2/sigma)*normpdf((log10(ncVals(j)) - mu)/sigma)*normcdf(alpha*(log10(ncVals(j)) - mu)/sigma);
    tauc(j) = Dmax*(Dfifty^Dk)/(nc^Dk + (Dfifty^Dk));
end

taup = 0.24;
taua = 0.75;

g = g./sum(g);

% Now scale the distribution of donors with each SPVL, g(V), by the infectious periods
% of the donors, since the expressions in Fraser et al. (2007) represent numbers of individuals of
% each viral load at seroconversion, rather than absolute numbers of individuals in the population

numIndividualTypes = length(g);
totalOfTimeVals = 0;
for i = 1:numIndividualTypes
    totalOfTimeVals = totalOfTimeVals + taup + tauc(i) + taua;
end

for i = 1:numIndividualTypes
    g(i) = g(i)*(taup + tauc(i) + taua)/totalOfTimeVals;
end

g = g./sum(g);

% Plot the proportion of donors of each SPVL (on a log scale) in the population
figure(2)
subplot(1,3,1)
bar(log10(ncVals), g)
ax = gca;
xtickPreCalc = [log10(ncVals(1)): ((log10(ncVals(length(ncVals))) - log10(ncVals(1)))/5): log10(ncVals(length(ncVals)))];
xtick = zeros(length(xtickPreCalc),1);
for i = 1:length(xtickPreCalc)
    xtick(i) = 10^xtickPreCalc(i);
end
ax.XTick = xtickPreCalc;
ax.XTickLabel = xtick;
xlim([1.5 7.5])
sum(g)
box off
title('Proportion of donors with each SPVL')
xlabel('SPVL')
ylabel('Proportion')

% Plot the length of chronic infection for donors of each SPVL (on a log scale) in the population
figure(2)
subplot(1,3,2)
bar(log10(ncVals), tauc)
ax = gca;
xtickPreCalc = [log10(ncVals(1)): ((log10(ncVals(length(ncVals))) - log10(ncVals(1)))/5): log10(ncVals(length(ncVals)))];
xtick = zeros(length(xtickPreCalc),1);
for i = 1:length(xtickPreCalc)
    xtick(i) = 10^xtickPreCalc(i);
end
ax.XTick = xtickPreCalc;
ax.XTickLabel = xtick;
xlim([1.5 7.5])
sum(g)
box off
title('Length of chronic infection for each SPVL')
xlabel('SPVL')
ylabel('Length of chronic infection (years)')

% Now plot viral load profiles (the number of particles available for transmission
% in donors as a function of time since infection)
Vp = 87173000;
Va = 24004000;

np = round(kappa*Vp);
na = round(kappa*Va);
maxInfectiousPeriod = taup + max(tauc)+ taua;
timestep = 0.1;

figure(2)
subplot(1,3,3)
iValsHere = [2:2:10];
logSPVLsHere = log10(ncVals(iValsHere));
x = colormap(hsv(10));
y = [];

for i = 1:length(logSPVLsHere)
    if i == 2
        y(i) = plot([0 taup], [log10(np) log10(np)], 'Color', 'k');
    else
        y(i) = plot([0 taup], [log10(np) log10(np)], 'Color', x(i,:));
    end
    hold on
    plot([taup taup], [log10(np) logSPVLsHere(i)], 'Color', x(i,:))
    hold on
    plot([taup (taup + tauc(iValsHere(i)))], [logSPVLsHere(i) logSPVLsHere(i)], 'Color', x(i,:))
    hold on
    plot([(taup + tauc(iValsHere(i))) (taup + tauc(iValsHere(i)))], [logSPVLsHere(i) log10(na)], 'Color', x(i,:))
    hold on
    plot([(taup + tauc(iValsHere(i))) (taup + tauc(iValsHere(i)) + taua)], [log10(na) log10(na)], 'Color', x(i,:))
    hold on
    plot([(taup + tauc(iValsHere(i)) + taua) (taup + tauc(iValsHere(i)) + taua)], [log10(na) -5000], 'Color', x(i,:))
    hold on
    plot([(taup + tauc(iValsHere(i)) + taua) (taup + max(tauc) + taua)], [-5000 -5000], 'Color', x(i,:))
end
ylim([1 8])
xlabel('Time since infection (years)')
ylabel('log(Viral load)')
title('Viral load profiles of donors in the population')
figtitleRT('Characterising the number of particles available for transmission in donors (cf Fig 2 of our publication)','fontweight','bold');


% Now we examine various population-scale quantities

p =  4.715*1e-8;  % Probablility of transmission per particle
f = 0.029;  % Fraction of the time when transmission possible

% First, the probability of transmitting n particles vs n

probNoTransmissionPerSexAct = 0;
for j = 1:length(ncVals)
    nc = ncVals(j);
    probNoTransmissionPerSexAct = probNoTransmissionPerSexAct + g(j)*((1 - f) + f*((taup/(taup + tauc(j) + taua))*(1 - p)^np + (tauc(j)/(taup + tauc(j) + taua))*(1 - p)^nc + (taua/(taup + tauc(j) + taua))*(1 - p)^na));
end
probTransmissionPerSexAct = 1 - probNoTransmissionPerSexAct;

probTransmitnparticles = zeros(10^7,1);
threshold = 10^(-6);
n = 1;

while sum((probTransmitnparticles)/(probTransmissionPerSexAct)) < (1 - threshold)
    p0 = 0;
    if mod(n,1000) == 0
        n
    end
    for j = 1:length(ncVals)
        nc = ncVals(j);
        if n <= np
            logBitToAddOn = log(taup/(taup + tauc(j) + taua)) + n*log(p) + (np - n)*log(1-p);
            if np > n
                for z = (np - n + 1):np
                    logBitToAddOn = logBitToAddOn + log(z);
                end
                for z = 1:n
                    logBitToAddOn = logBitToAddOn - log(z);
                end
            end
            p0 = p0 + exp(logBitToAddOn)*g(j)*f;
        end
        if n <= nc
            logBitToAddOn = log(tauc(j)/(taup + tauc(j) + taua)) + n*log(p) + (nc - n)*log(1-p);
            if nc > n
                for z = (nc - n + 1):nc
                    logBitToAddOn = logBitToAddOn + log(z);
                end
                for z = 1:n
                    logBitToAddOn = logBitToAddOn - log(z);
                end
            end
            p0 = p0 + exp(logBitToAddOn)*g(j)*f;
        end
        if n <= na
            logBitToAddOn = log(taua/(taup + tauc(j) + taua)) + n*log(p) + (na - n)*log(1-p);
            if na > n
                for z = (na - n + 1):na
                    logBitToAddOn = logBitToAddOn + log(z);
                end
                for z = 1:n
                    logBitToAddOn = logBitToAddOn - log(z);
                end
            end
            p0 = p0 + exp(logBitToAddOn)*g(j)*f;
        end
    end
    probTransmitnparticles(n) = p0;
    n = n + 1;
end
probTransmitnparticles = probTransmitnparticles(1:(n-1));
probTransmitnparticles = probTransmitnparticles./sum(probTransmitnparticles);

% Now, the probability of transmitting N variants vs N
nparticlesConsidered = length(probTransmitnparticles);
maximumTime = max(taup + tauc + taua);

nTimeSteps = 1000;
timeWindowEdges = [0:maximumTime/nTimeSteps:maximumTime];
timeVals = zeros(nTimeSteps,1);
probTransmitNvariantsGivenTimeTAndTransmitNparticles = zeros(nparticlesConsidered, nTimeSteps, nparticlesConsidered);

disp('Calculating the probability of transmitting N variants')
disp('given a distribution of variants and number of particles')
disp('transmitted. This may take some time...')

timeSinceInfectionBeingCalculated = 0;
for timeV = 1:nTimeSteps
    timeVals(timeV) = (timeWindowEdges(timeV) + timeWindowEdges(timeV + 1))/2;
    
    if timeVals(timeV) > timeSinceInfectionBeingCalculated
        timeSinceInfectionBeingCalculated
        timeSinceInfectionBeingCalculated = timeSinceInfectionBeingCalculated + 1;
    end
    
    probDist = zeros(maxnvariants,1);
    for i = 1:maxnvariants
        probDist(i) = gampdf(i, 0.417, timeVals(timeV)/0.563);
    end
    probDist = probDist/sum(probDist);
    
    nSims = 100000;
    
    for simNo = 1:nSims
        variantIIndicator = zeros(maxnvariants,1);
        for nparticles = 1:nparticlesConsidered
            
            variantPicked = 1;
            variantPickedCount = 0;
            r = rand();
            while r > variantPickedCount
                variantPickedCount = variantPickedCount + probDist(variantPicked);
                variantPicked = variantPicked + 1;
            end
            variantIIndicator(variantPicked - 1) = 1;
            
            nvariantsTransferred = sum(variantIIndicator);
            probTransmitNvariantsGivenTimeTAndTransmitNparticles(nvariantsTransferred, timeV, nparticles) = probTransmitNvariantsGivenTimeTAndTransmitNparticles(nvariantsTransferred, timeV, nparticles) + 1;
            
        end
    end
end
probTransmitNvariantsGivenTimeTAndTransmitNparticles = probTransmitNvariantsGivenTimeTAndTransmitNparticles./nSims;
probTransmitMvariants = zeros(min(maxnvariants, nparticlesConsidered),1);

for M = 1:min(maxnvariants, nparticlesConsidered)
    
    for i = 1:length(ncVals)
        
        integralPrimary = 0;
        for j = 1:nTimeSteps
            if timeVals(j) <= taup
                for nVal = M:min(np, nparticlesConsidered)
                    binomialBit = 1;
                    if np > nVal
                        logBinBit = nVal*log(p) + (np - nVal)*log(1-p);
                        for z = (np - nVal + 1):np
                            logBinBit = logBinBit + log(z);
                        end
                        for z = 1:nVal
                            logBinBit = logBinBit - log(z);
                        end
                        binomialBit = exp(logBinBit);
                    end
                    integralPrimary = integralPrimary + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                end
            end
        end
        
        nc = ncVals(i);
        integralChronic = 0;
        
        for j = 1:nTimeSteps
            
            if ((timeVals(j) > taup) && (timeVals(j) <= (taup + tauc(i))))
                for nVal = M:min(nc, nparticlesConsidered)
                    
                    binomialBit = 1;
                    if nc > nVal
                        logBinBit = nVal*log(p) + (nc - nVal)*log(1-p);
                        for z = (nc - nVal + 1):nc
                            logBinBit = logBinBit + log(z);
                        end
                        for z = 1:nVal
                            logBinBit = logBinBit - log(z);
                        end
                        binomialBit = exp(logBinBit);
                    end
                    integralChronic = integralChronic + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                end
            end
        end
        
        integralPreAids = 0;
        
        for j = 1:nTimeSteps
            if ((timeVals(j) > (taup + tauc(i))) && (timeVals(j) <= (taup + tauc(i) + taua)))
                for nVal = M:min(na, nparticlesConsidered)
                    
                    binomialBit = 1;
                    if na > nVal
                        logBinBit = nVal*log(p) + (na - nVal)*log(1-p);
                        for z = (na - nVal + 1):na
                            logBinBit = logBinBit + log(z);
                        end
                        for z = 1:nVal
                            logBinBit = logBinBit - log(z);
                        end
                        binomialBit = exp(logBinBit);
                    end
                    integralPreAids = integralPreAids + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                end
            end
        end
        
        probTransmitMvariants(M) = probTransmitMvariants(M) + g(i)*(integralPrimary + integralChronic + integralPreAids);
        
    end
end
probTransmitMvariants = probTransmitMvariants/sum(probTransmitMvariants);

probMultivariantTransmission = 1 - probTransmitMvariants(1);

figure(3)
bar(probTransmitMvariants, 'barWidth', 1, 'FaceColor', [0.2 0.2 0.5])
hold on
bar(probTransmitnparticles, 'barWidth', 0.5, 'FaceColor', [0 0.7 0.7])
xlim([0.5 16.5])
title('Probability of transmitting n particles (teal)/N variants (blue)')
xlabel('Number (n/N)')
ylabel('Probability')
ylim([0 1])
figtitleRT('The numbers of particles and variants transmitted (cf top left panel of Fig 4 of our publication)','fontweight','bold');

% Now, the joint probability distribution of transmitting N variants and n particles
probTransmitMvariantsAndNparticles = zeros(min(maxnvariants, nparticlesConsidered), nparticlesConsidered);
probTransmitMvariantsAndNparticlesPrimary = zeros(min(maxnvariants, nparticlesConsidered), nparticlesConsidered);
probTransmitMvariantsAndNparticlesChronic = zeros(min(maxnvariants, nparticlesConsidered), nparticlesConsidered);
probTransmitMvariantsAndNparticlesPreAIDS = zeros(min(maxnvariants, nparticlesConsidered), nparticlesConsidered);

for nVal = 1:nparticlesConsidered
    for M = 1:min(maxnvariants, nparticlesConsidered)
        for i = 1:length(ncVals)
            
            integralPrimary = 0;
            for j = 1:nTimeSteps
                if timeVals(j) <= taup
                    if nVal >= M && nVal <= min(np, nparticlesConsidered)
                        binomialBit = 1;
                        if np > nVal
                            logBinBit = nVal*log(p) + (np - nVal)*log(1-p);
                            for z = (np - nVal + 1):np
                                logBinBit = logBinBit + log(z);
                            end
                            for z = 1:nVal
                                logBinBit = logBinBit - log(z);
                            end
                            binomialBit = exp(logBinBit);
                        end
                        integralPrimary = integralPrimary + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                    end
                end
            end
            
            nc = ncVals(i);
            integralChronic = 0;
            
            for j = 1:nTimeSteps
                if ((timeVals(j) > taup) && (timeVals(j) <= (taup + tauc(i))))
                    if (nVal >= M) && (nVal <= min(nc, nparticlesConsidered))
                        
                        binomialBit = 1;
                        if nc > nVal
                            logBinBit = nVal*log(p) + (nc - nVal)*log(1-p);
                            for z = (nc - nVal + 1):nc
                                logBinBit = logBinBit + log(z);
                            end
                            for z = 1:nVal
                                logBinBit = logBinBit - log(z);
                            end
                            binomialBit = exp(logBinBit);
                        end
                        integralChronic = integralChronic + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                    end
                end
            end
            
            integralPreAids = 0;
            
            for j = 1:nTimeSteps
                if ((timeVals(j) > (taup + tauc(i))) && (timeVals(j) <= (taup + tauc(i) + taua)))
                    if (nVal >= M) && (nVal <= min(na, nparticlesConsidered))
                        binomialBit = 1;
                        if na > nVal
                            logBinBit = nVal*log(p) + (na - nVal)*log(1-p);
                            for z = (na - nVal + 1):na
                                logBinBit = logBinBit + log(z);
                            end
                            for z = 1:nVal
                                logBinBit = logBinBit - log(z);
                            end
                            binomialBit = exp(logBinBit);
                        end
                        integralPreAids = integralPreAids + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                    end
                end
            end
            
            probTransmitMvariantsAndNparticles(M,nVal) = probTransmitMvariantsAndNparticles(M,nVal) + g(i)*(integralPrimary + integralChronic + integralPreAids);
            probTransmitMvariantsAndNparticlesPrimary(M,nVal) = probTransmitMvariantsAndNparticlesPrimary(M,nVal) + g(i)*(integralPrimary);
            probTransmitMvariantsAndNparticlesChronic(M,nVal) = probTransmitMvariantsAndNparticlesChronic(M,nVal) + g(i)*(integralChronic);
            probTransmitMvariantsAndNparticlesPreAIDS(M,nVal) = probTransmitMvariantsAndNparticlesPreAIDS(M,nVal) + g(i)*(integralPreAids);
        end
    end
end

probTransmitMvariantsAndNparticles = probTransmitMvariantsAndNparticles./(sum(sum(probTransmitMvariantsAndNparticles)));
probTransmitMvariantsAndNparticlesPrimary = probTransmitMvariantsAndNparticlesPrimary./(sum(sum(probTransmitMvariantsAndNparticles)));
probTransmitMvariantsAndNparticlesChronic = probTransmitMvariantsAndNparticlesChronic./(sum(sum(probTransmitMvariantsAndNparticles)));
probTransmitMvariantsAndNparticlesPreAIDS = probTransmitMvariantsAndNparticlesPreAIDS./(sum(sum(probTransmitMvariantsAndNparticles)));

figure(4)
for i = 1:min(maxnvariants, nparticlesConsidered)
    for j = 1:nparticlesConsidered
        if probTransmitMvariantsAndNparticles(i,j) > 0
            hold on
            plot([j],[i],'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize', 55*sqrt(probTransmitMvariantsAndNparticles(i,j)))
        end
    end
end
title('The joint probability distribution for the numbers of particles and variants founding new infections')
xlabel('Number of particles')
ylabel('Number of variants')
figtitleRT('The numbers of particles and variants transmitted (cf top middle panel of Fig 4 of our publication)','fontweight','bold');


% Now the distributions of particles and variants transmitted when the donor is in early infection (donor infected for less than 2 years)
probTransmitMvariantsAndNparticlesFirstTwoYears = zeros(min(maxnvariants, nparticlesConsidered), nparticlesConsidered);

for nVal = 1:nparticlesConsidered
    for M = 1:min(maxnvariants, nparticlesConsidered)
        for i = 1:length(ncVals)
            
            integralPrimary = 0;
            for j = 1:nTimeSteps
                if ((timeVals(j) > 0) && (timeVals(j) <= 2))
                    if timeVals(j) <= taup
                        if nVal >= M && nVal <= min(np, nparticlesConsidered)
                            
                            binomialBit = 1;
                            if np > nVal
                                logBinBit = nVal*log(p) + (np - nVal)*log(1-p);
                                for z = (np - nVal + 1):np
                                    logBinBit = logBinBit + log(z);
                                end
                                for z = 1:nVal
                                    logBinBit = logBinBit - log(z);
                                end
                                binomialBit = exp(logBinBit);
                            end
                            integralPrimary = integralPrimary + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                        end
                    end
                end
            end
            
            nc = ncVals(i);
            integralChronic = 0;
            
            for j = 1:nTimeSteps
                if min(2, taup) < min(taup + tauc(i), 2)
                    if ((timeVals(j) > min(2, taup)) && (timeVals(j) <= (min(taup + tauc(i), 2)))) 
                        if ((timeVals(j) > taup) && (timeVals(j) <= (taup + tauc(i))))
                            if (nVal >= M) && (nVal <= min(nc, nparticlesConsidered))
                                binomialBit = 1;
                                if nc > nVal
                                    logBinBit = nVal*log(p) + (nc - nVal)*log(1-p);
                                    for z = (nc - nVal + 1):nc
                                        logBinBit = logBinBit + log(z);
                                    end
                                    for z = 1:nVal
                                        logBinBit = logBinBit - log(z);
                                    end
                                    binomialBit = exp(logBinBit);
                                end
                                integralChronic = integralChronic + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                            end
                        end
                    end
                end
            end
            
            integralPreAids = 0;
            for j = 1:nTimeSteps
                if 2 > min(taup + tauc(i) + taua, 2)
                    if ((timeVals(j) > min(taup + tauc(i) + taua, 2)) && (timeVals(j) <= (2)))
                        if ((timeVals(j) > (taup + tauc(i))) && (timeVals(j) <= (taup + tauc(i) + taua)))
                            if (nVal >= M) && (nVal <= min(na, nparticlesConsidered))
                                binomialBit = 1;
                                if na > nVal
                                    logBinBit = nVal*log(p) + (na - nVal)*log(1-p);
                                    for z = (na - nVal + 1):na
                                        logBinBit = logBinBit + log(z);
                                    end
                                    for z = 1:nVal
                                        logBinBit = logBinBit - log(z);
                                    end
                                    binomialBit = exp(logBinBit);
                                end
                                integralPreAids = integralPreAids + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                            end
                        end
                    end
                end
            end
            
            probTransmitMvariantsAndNparticlesFirstTwoYears(M,nVal) = probTransmitMvariantsAndNparticlesFirstTwoYears(M,nVal) + g(i)*(integralPrimary + integralChronic + integralPreAids);
            
        end
    end
end

probTransmitMvariantsAndNparticlesFirstTwoYears = probTransmitMvariantsAndNparticlesFirstTwoYears./(sum(sum(probTransmitMvariantsAndNparticlesFirstTwoYears)));

probTransmitMvariantsFirstTwoYears = sum(transpose(probTransmitMvariantsAndNparticlesFirstTwoYears));
probTransmitnparticlesFirstTwoYears = sum(probTransmitMvariantsAndNparticlesFirstTwoYears);
figure(5)
bar(probTransmitMvariantsFirstTwoYears, 'barWidth', 1, 'FaceColor', [0.2 0.2 0.5])
hold on
bar(probTransmitnparticlesFirstTwoYears, 'barWidth', 0.5, 'FaceColor', [0 0.7 0.7])
xlim([0.5 16.5])
xlabel('Number (n/N)')
ylabel('Probability')
title({'Probability distributions for the number of particles (teal)/variants (blue)'; 'transmitted, for donors in the first two years of infection'})
figtitleRT('The numbers of particles and variants transmitted in early infection (cf top right panel of Fig 4 of our publication)','fontweight','bold');


% Now, the probability that a randomly chosen infection is with M variants and arises from a donor with each SPVL

probTransmitMvariantsAndSPVLx = zeros(min(maxnvariants, nparticlesConsidered), length(g));
for nVal = 1:nparticlesConsidered
    for M = 1:min(maxnvariants, nparticlesConsidered)
        for i = 1:length(ncVals)
            
            integralPrimary = 0;
            
            for j = 1:nTimeSteps
                if timeVals(j) <= taup
                    if nVal >= M && nVal <= min(np, nparticlesConsidered)
                        binomialBit = 1;
                        if np > nVal
                            logBinBit = nVal*log(p) + (np - nVal)*log(1-p);
                            for z = (np - nVal + 1):np
                                logBinBit = logBinBit + log(z);
                            end
                            for z = 1:nVal
                                logBinBit = logBinBit - log(z);
                            end
                            binomialBit = exp(logBinBit);
                        end
                        integralPrimary = integralPrimary + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                    end
                end
            end
            
            nc = ncVals(i);
            integralChronic = 0;
            
            for j = 1:nTimeSteps
                if ((timeVals(j) > taup) && (timeVals(j) <= (taup + tauc(i))))
                    if (nVal >= M) && (nVal <= min(nc, nparticlesConsidered))
                        binomialBit = 1;
                        if nc > nVal
                            logBinBit = nVal*log(p) + (nc - nVal)*log(1-p);
                            for z = (nc - nVal + 1):nc
                                logBinBit = logBinBit + log(z);
                            end
                            for z = 1:nVal
                                logBinBit = logBinBit - log(z);
                            end
                            binomialBit = exp(logBinBit);
                        end
                        integralChronic = integralChronic + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                    end
                end
            end
            
            integralPreAids = 0;
            
            for j = 1:nTimeSteps
                if ((timeVals(j) > (taup + tauc(i))) && (timeVals(j) <= (taup + tauc(i) + taua)))
                    if (nVal >= M) && (nVal <= min(na, nparticlesConsidered))
                        binomialBit = 1;
                        if na > nVal
                            logBinBit = nVal*log(p) + (na - nVal)*log(1-p);
                            for z = (na - nVal + 1):na
                                logBinBit = logBinBit + log(z);
                            end
                            for z = 1:nVal
                                logBinBit = logBinBit - log(z);
                            end
                            binomialBit = exp(logBinBit);
                        end
                        integralPreAids = integralPreAids + (timeVals(2) - timeVals(1))*probTransmitNvariantsGivenTimeTAndTransmitNparticles(M,j,nVal)*(1/(taup + tauc(i) + taua))*binomialBit;
                    end
                end
            end
            
            probTransmitMvariantsAndSPVLx(M,i) = probTransmitMvariantsAndSPVLx(M,i) + g(i)*(integralPrimary + integralChronic + integralPreAids);
            
        end
    end
end

probTransmitMvariantsAndSPVLx = probTransmitMvariantsAndSPVLx./(sum(sum(probTransmitMvariantsAndSPVLx)));

figure(6)
for i = 1:min(maxnvariants, nparticlesConsidered)
    for j = 1:length(g)
        if probTransmitMvariantsAndSPVLx(i,j) > 0
            hold on
            plot([j],[i],'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerSize', 55*sqrt(probTransmitMvariantsAndSPVLx(i,j)))
        end
    end
end
title({'The joint probability distribution for the numbers of'; 'particles founding new infections and donor SPVL'})
xlabel('Donor log(SPVL)')
xticklabels({'2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'})
ylabel('Number of variants')
figtitleRT('The relationship between the numbers of particles transmitted and donor SPVL (cf top left panel of Fig 5 of our publication)','fontweight','bold');

% Probability of transmission vs SPVL
probNoTransmission = zeros(length(g),1);
for i = 1:length(g)
    nc = ncVals(i);
    probNoTransmission(i) = (1 - f) + f*((taup/(taup + tauc(i) + taua))*(1 - p)^(np) + (tauc(i)/(taup + tauc(i) + taua))*(1 - p)^(nc) + (taua/(taup + tauc(i) + taua))*(1 - p)^(na));
end
figure(7)
subplot(2,2,1)
bar(1 - probNoTransmission, 'barWidth', 1)
xlim([0.5 11.5])
title('Probability of transmission vs SPVL')
xlabel('Donor log(SPVL)')
ylabel('Probability of transmission per sex act')
xticklabels({'2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'})

% Probability a randomly chosen infection comes from each SPVL
probInfectionFromEachSPVL = g.*(1 - probNoTransmission);
probInfectionFromEachSPVL = probInfectionFromEachSPVL./sum(probInfectionFromEachSPVL);
figure(7)
subplot(2,2,2)
bar(probInfectionFromEachSPVL, 'barWidth', 1)
xlim([0.5 11.5])
title({'Probability a randomly chosen'; 'infection comes from each SPVL'})
xlabel('Donor log(SPVL)')
ylabel({'Probability randomly chosen','transmission arose from SPVL x'})
xticklabels({'2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'})

% Probability a transmission from each SPVL is multi-variant
probTransmissionSPVLXMultivariant = zeros(length(g),1);
for i = 1:11
    probTransmissionSPVLXMultivariant(i) = 1 - (sum(probTransmitMvariantsAndSPVLx(1,i))/sum(probTransmitMvariantsAndSPVLx(:,i)));
end
figure(7)
subplot(2,2,3)
bar(probTransmissionSPVLXMultivariant, 'barWidth', 1)
hold on
plot([0 12], [0.3 0.3], 'k--')
xlim([0.5 11.5])
title({'Probability a transmission from a donor'; 'with each SPVL is multi-variant'})
xlabel('Donor log(SPVL)')
ylabel('Probability infection is multi-variant')
xticklabels({'2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'})

% Probability a randomly chosen multi variant transmission arose from SPVL x
probRandomInfectionCameFromEachSPVL = probInfectionFromEachSPVL.*probTransmissionSPVLXMultivariant;
figure(7)
subplot(2,2,4)
bar(probRandomInfectionCameFromEachSPVL./sum(probRandomInfectionCameFromEachSPVL), 'barWidth', 1)
xlim([0.5 11.5])
title({'Probability a randomly chosen multi-variant'; 'transmission arose from each SPVL'})
xlabel('Donor log(SPVL)')
ylabel({'Probability randomly chosen multi-'; 'variant transmission arose from SPVL x'})
xticklabels({'2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'})
figtitleRT('The relationship between multi-variant transmissions and donor SPVL (cf Fig 5 of our publication)','fontweight','bold');





























function [ fth ] = figtitleRT(titlestring,varargin)
% Adding titles to figure with several subplots,
% based on code by Chad A. Greene of the University of Texas at Austin
hca = gca; 
fontsize = get(hca,'fontsize'); 
h = axes('position',[0 0 1 1],'units','normalized');
axes('Units','normalized',...
    'Position',[0 0 1 1],...
    'Visible','off',...
    'XTick',[],...
    'YTick',[],...
    'Box','off');
fth = text(.5,1,titlestring,...
    'units','normalized',...
    'horizontalalignment','center',...
    'verticalalignment','top',...
    'fontsize',fontsize+2); 
if nargin>1
    set(fth,varargin{:});
end
delete(h)
set(gcf,'CurrentAxes',hca,'name',titlestring); 
if nargout==0
    clear fth; 
end
end
