%Script called by neuronalAnalysis to calculate various indices for number of channels active
%when we didn't necessarily use stimFreqs = [0 10 100 150 250 350 425]
hotVal = stimFreqs(7);
warmBounds = [stimFreqs(3) stimFreqs(6)];
coldUpperBound = stimFreqs(2);

ir1HotInds = find(allStimvecsFreq(:,1) == hotVal );
ir1WarmInds = find(allStimvecsFreq(:,1) >= warmBounds(1) & allStimvecsFreq(:,1) <= warmBounds(2));
ir1ColdInds = find(allStimvecsFreq(:,1) <= coldUpperBound);

ir2HotInds = find(allStimvecsFreq(:,2) == hotVal);
ir2WarmInds = find(allStimvecsFreq(:,2) >= warmBounds(1) & allStimvecsFreq(:,2) <= warmBounds(2));
ir2ColdInds = find(allStimvecsFreq(:,2) <= coldUpperBound);

ir3HotInds = find(allStimvecsFreq(:,3) == hotVal);
ir3WarmInds = find(allStimvecsFreq(:,3) >= warmBounds(1) & allStimvecsFreq(:,3) <= warmBounds(2));
ir3ColdInds = find(allStimvecsFreq(:,3) <= coldUpperBound);

ir4HotInds = find(allStimvecsFreq(:,4) == hotVal);
ir4WarmInds = find(allStimvecsFreq(:,4) >= warmBounds(1) & allStimvecsFreq(:,4) <= warmBounds(2));
ir4ColdInds = find(allStimvecsFreq(:,4) <= coldUpperBound);


%% Doublets 
ir12WarmInds = intersect(ir1WarmInds, ir2WarmInds); size(ir12WarmInds); %THIS WAS UNION
ir34ColdInds = intersect(ir3ColdInds, ir4ColdInds); size(ir34ColdInds); %THIS WAS UNION
irFrontInds = intersect(ir12WarmInds, ir34ColdInds);
size(irFrontInds);

ir34WarmInds = intersect(ir3WarmInds, ir4WarmInds);  %THIS WAS UNION
ir12ColdInds = intersect(ir1ColdInds, ir2ColdInds);  %THIS WAS UNION
irBackInds = intersect(ir34WarmInds, ir12ColdInds);

ir23WarmInds = intersect(ir2WarmInds, ir3WarmInds);  %THIS WAS UNION
ir14ColdInds = intersect(ir1ColdInds, ir4ColdInds);  %THIS WAS UNION
irLeftInds = intersect(ir23WarmInds, ir14ColdInds);

ir14WarmInds = intersect(ir1WarmInds, ir4WarmInds);  %THIS WAS UNION
ir23ColdInds = intersect(ir2ColdInds, ir3ColdInds);  %THIS WAS UNION
irRightInds = intersect(ir14WarmInds, ir23ColdInds);

%weird doublets  %THIS WAS not here
ir13WarmInds = intersect(ir1WarmInds, ir3WarmInds);
ir24ColdInds = intersect(ir2ColdInds, ir4ColdInds);
irDiagonalRightInds = intersect(ir13WarmInds, ir24ColdInds);

ir24WarmInds = intersect(ir2WarmInds, ir4WarmInds);
ir13ColdInds = intersect(ir1ColdInds, ir3ColdInds);
irDiagonalLeftInds = intersect(ir24WarmInds, ir13ColdInds);

numDoubletsInds = length(irFrontInds) + length(irBackInds) + length(irLeftInds) ...
                + length(irRightInds) + length(irDiagonalRightInds) + length(irDiagonalLeftInds);


%% Triplets (allStimvecsFreq)
ir123WarmInds = intersect(ir1WarmInds, intersect(ir2WarmInds, ir3WarmInds));
irFrontLeftInds = intersect(ir123WarmInds, ir4ColdInds);

ir124WarmInds = intersect(ir1WarmInds, intersect(ir2WarmInds, ir4WarmInds));
irFrontRightInds  = intersect(ir124WarmInds, ir3ColdInds);

ir234WarmInds = intersect(ir2WarmInds, intersect(ir3WarmInds, ir4WarmInds));
irBackLeftInds = intersect(ir234WarmInds, ir1ColdInds);  %usually empty

ir134WarmInds = intersect(ir1WarmInds, intersect(ir3WarmInds, ir4WarmInds));
irBackRightInds = intersect(ir134WarmInds, ir2ColdInds); %usually empty

numTripletsInds = length(irFrontLeftInds) + length(irFrontRightInds) ...
                + length(irBackLeftInds)  + length(irBackRightInds); 



%% All channels active
irAllInds = intersect(ir123WarmInds, ir4WarmInds);
numAllInds = length(irAllInds);

%% Individuals

ir234ColdInds = intersect(ir2ColdInds, intersect(ir3ColdInds, ir4ColdInds));
ir1WarmIndsHotInds = union(ir1WarmInds, ir1HotInds);
ir1OnlyInds = intersect(ir1WarmInds, ir234ColdInds);  %There are just not enough like this!
size(ir1OnlyInds) ; %28 substim :0

ir134ColdInds = intersect(ir1ColdInds, intersect(ir3ColdInds, ir4ColdInds));
ir2WarmHotInds = union(ir2WarmInds, ir2HotInds);
ir2OnlyInds = intersect(ir2WarmInds, ir134ColdInds);  %There are just not enough like this!
size(ir2OnlyInds);  %50 substim

ir124ColdInds = intersect(ir1ColdInds, intersect(ir2ColdInds, ir4ColdInds));
ir3WarmHotInds = union(ir3WarmInds, ir3HotInds);
ir3OnlyInds = intersect(ir3WarmInds, ir124ColdInds);  %There are just not enough like this!
size(ir3OnlyInds);  %88 substim

ir123ColdInds = intersect(ir1ColdInds, intersect(ir2ColdInds, ir3ColdInds));
ir4WarmHotInds = union(ir4WarmInds, ir4HotInds);
ir4OnlyInds = intersect(ir4WarmInds, ir123ColdInds);  %There are just not enough like this!
size(ir4OnlyInds);  %173 substim

numIndividualsInds = length(ir1OnlyInds) + length(ir2OnlyInds) + length(ir3OnlyInds) + length(ir4OnlyInds);


%% All off
irOffInds = intersect(ir123ColdInds, ir4ColdInds); 
numOffInds = length(irOffInds);



%% calculate proportion of each type
proportionChannelsActivated = [numOffInds numIndividualsInds numDoubletsInds numTripletsInds numAllInds]/totalNumSubstim;
if plotIndividualChannels
    figure;
    bar([0:4], proportionChannelsActivated);
    xlabel('Number of channels activated');
    ylabel('Proportion of substim');
    if printOn
        print
    end
end

%% Error check (recall total number substim: totalNumSubstim)
numWarmColdInds = numDoubletsInds + numTripletsInds + numAllInds + numIndividualsInds + numOffInds;
%get number with some hot (recall allStimvecsFreq is totalNumSubstim x 4 )
[hotRowInds, hotColInds] = find(allStimvecsFreq == hotVal);
numHotInds = length(unique(hotRowInds));
numSubstimEstimateInds = numWarmColdInds + numHotInds;
%Check total estimate against total actual
if numSubstimEstimateInds == totalNumSubstim
    display('Bullseye! Number of substim estimated in calcWarmHotColdInds_indexBased matches actual!');
else
    error('Substim estimate in calcWarmHotColdInds_indexBased is off!');
end





