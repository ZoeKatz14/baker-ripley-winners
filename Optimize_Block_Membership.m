%% Load data
clc, clear
data = readtable("Block Data v1.csv");
demoMat = table2array(data(:, 7 : 28));
blockX = (data.INTPTLON);
blockY = (data.INTPTLAT);

numBlock = length(blockX);

%% Run Voronoi Difference
clc
numDist = 11;
[optim, objVal] = optimizeVoronoiDifference(blockX, blockY, demoMat, numDist);

%% Visualize
clc
memberMat = centerVec2MemberMat(optim, blockX, blockY);
[fig, memberMat] = visualizeGeographic(memberMat, numDist, numBlock, blockX, blockY);
legend(string(1:11))
[largestDifference, difference, districtDemo] = computeDemographicDifference(memberMat, demoMat);

%% Edges
clc
centerX = optim(1 : numDist);
centerY = optim(numDist + 1 : 2 * numDist);
[vx, vy] = voronoi(centerX, centerY);
for ind =  1 : width(vx)
    geoplot(vy(:, ind), vx(:, ind))
    hold on
end
hold off

%% Export labels
label = memberMat2Labels(memberMat)
writematrix(label, "Block Labels v1.csv")

%% Functions

function label = memberMat2Labels(memberMat)
numBlocks = width(memberMat);
label = NaN(numBlocks, 1);
for ind = 1 : numBlocks
    label(ind) = find(memberMat(:, ind) == 1);
end
end

function ratio = computeCompactness(x, y)
[dBound, dArea] = boundary(x, y, 1);
[chBound, chArea] = convhull(x, y);
ratio = dArea / chArea;
end

function avgRatio = computeTotalCompactness(memberMat, numDist, blockX, blockY)
avgRatio = 0;
for distInd = 1 : numDist
    targetBlocks = memberMat(distInd, :) == 1;
    x = blockX(targetBlocks);
    y = blockY(targetBlocks);
    avgRatio = avgRatio + computeCompactness(x, y);
end
avgRatio = avgRatio / width(memberMat);
end

function [sumVariance, variance, districtDemo] = computeDemographicVariance(memberMat, demoMat)
districtDemo = memberMat * demoMat;
variance = var(districtDemo);
sumVariance = sum(variance);
end

function [largestDifference, difference, districtDemo] = computeDemographicDifference(memberMat, demoMat)
districtDemo = memberMat * demoMat;
highDemo = max(districtDemo);
lowDemo = min(districtDemo);
difference = (highDemo - lowDemo) ./ (highDemo + lowDemo) * 2;
largestDifference = max(abs(difference));
end

function memberMat = memberVec2Mat(memberVec, numDist, numBlock)
memberMat = NaN(numDist, numBlock);
for distInd = 1 : numDist
    (distInd - 1) * numBlock + 1
    distInd * numBlock
    memberMat(distInd, :) = memberVec((distInd - 1) * numBlock + 1 : distInd * numBlock); 
end
end

function objVal = multiObjFun(memberVec, numDist, numBlock, demoMat, blockX, blockY)
memberMat = memberVec2Mat(memberVec, numDist, numBlock);
[sumVariance, ~, ~] = computeDemographicVariance(memberMat, demoMat);
avgRatio = computeTotalCompactness(memberMat, numDist, blockX, blockY);
objVal = [sumVariance; avgRatio];
end

function objVal = objFun(memberVec, numDist, numBlock, demoMat, blockX, blockY, weight)
memberMat = memberVec2Mat(memberVec, numDist, numBlock);
[sumVariance, ~, ~] = computeDemographicVariance(memberMat, demoMat);
avgRatio = computeTotalCompactness(memberMat, numDist, blockX, blockY);
objVal = sumVariance * weight(1) + avgRatio * weight(2);
end


function pareto = runMultiObjGenAlg(numDist, numBlock, demoMat, blockX, blockY)
fun = @(memberVec) multiObjFun(memberVec, numDist, numBlock, demoMat, blockX, blockY);
nvars = numDist * numBlock;
A = [];
b = [];
beq = ones(numBlock, 1);
lb = zeros(nvar, 1);
ub = zeros(nvar, 1);
nonlcon = [];
intcon = 1 : nvars;

Aeq = [];
for distInd = 1 : numDist
    Aeq = horzcat(Aeq, eye(numBlock));
end

pareto = gamultiobj(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon);
end


function [optim, objVal] = runGenAlg(numDist, numBlock, demoMat, blockX, blockY, objWeight)
fun = @(memberVec) objFun(memberVec, numDist, numBlock, demoMat, blockX, blockY, objWeight);
nvars = numDist * numBlock;
A = [];
b = [];
beq = ones(numBlock, 1);
lb = zeros(nvar, 1);
ub = zeros(nvar, 1);
nonlcon = [];
intcon = 1 : nvars;

Aeq = [];
for distInd = 1 : numDist
    Aeq = horzcat(Aeq, eye(numBlock));
end

[optim, objVal] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon);
end

function [fig, memberMat] = visualizeGeographic(memberMat, numDist, numBlock, blockX, blockY)
labels = (1 : numDist) * memberMat;

fig = figure();
for distInd = 1 : numDist
    targetBlocks = find(labels == distInd);
    geoscatter(blockY(targetBlocks), blockX(targetBlocks));
    hold on
end
hold off
end

function [distancesSq, nearDist] = computeDistancesSq(numDist, numBlock, centerX, centerY, blockX, blockY)
distancesSq = NaN(numDist, numBlock);
for blockInd = 1 : numBlock
    distancesSq(:, blockInd) = (blockX(blockInd) - centerX') .^ 2 + (blockY(blockInd) - centerY') .^ 2;
end
[distancesSq, nearDist] = min(distancesSq);
end

function sumVariance = voronoiObjVariance(centerVec, blockX, blockY, demoMat)
memberMat = centerVec2MemberMat(centerVec, blockX, blockY);
[sumVariance, ~, ~] = computeDemographicVariance(memberMat, demoMat);
end

function largestDifference = voronoiObjDifference(centerVec, blockX, blockY, demoMat)
memberMat = centerVec2MemberMat(centerVec, blockX, blockY);
[largestDifference, ~, ~] = computeDemographicDifference(memberMat, demoMat);
end

function memberMat = centerVec2MemberMat(centerVec, blockX, blockY)
numDist = length(centerVec) / 2;
centerX = centerVec(1 : numDist);
centerY = centerVec(numDist + 1 : 2 * numDist);
numBlock = length(blockX);
[distancesSq, nearDist] = computeDistancesSq(numDist, numBlock, centerX, centerY, blockX, blockY);
memberMat = zeros(numDist, numBlock);
for blockInd = 1 : length(nearDist)
    memberMat(nearDist(blockInd), blockInd) = 1;
end
end

function [optim, objVal] = optimizeVoronoi(blockX, blockY, demoMat, numDist)
fun = @(centerVec) voronoiObjVariance(centerVec, blockX, blockY, demoMat);
nvars = numDist * 2;
A = [];
b = [];
Aeq = [];
beq = [];
lb = vertcat((-95.85 * ones(numDist, 1)), (29.53 * ones(numDist, 1)));
ub = vertcat((-94.95 * ones(numDist, 1)), (30.15 * ones(numDist, 1)));
nonlcon = [];
intcon = [];

options=optimoptions('GA','Display','iter');

[optim, objVal] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
end

function [optim, objVal] = optimizeVoronoiDifference(blockX, blockY, demoMat, numDist)
fun = @(centerVec) voronoiObjDifference(centerVec, blockX, blockY, demoMat);
nvars = numDist * 2;
A = [];
b = [];
Aeq = [];
beq = [];
lb = vertcat((-95.85 * ones(numDist, 1)), (29.53 * ones(numDist, 1)));
ub = vertcat((-94.95 * ones(numDist, 1)), (30.15 * ones(numDist, 1)));
nonlcon = [];
intcon = [];

options=optimoptions('GA','Display','iter');

[optim, objVal] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
end