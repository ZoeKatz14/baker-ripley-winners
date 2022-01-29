%% Load data
clc, clear
data = readtable("Block Group Data v0.csv");
demograph = table2array(data(:, 7 : 28));
blockX = (data.INTPTLON);
blockY = (data.INTPTLAT);

numBlock = length(blockX);

%% Run
clc
numDist = 11;
diffTol = 0.05;
[fun,nvar,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options] = initGA(numDist, numBlock, demograph, blockX, blockY, diffTol);
[optim, obj] = runGA(fun,nvar,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);

%% Functions
function memberMat = memberVec2Mat(memberVec, numDist, numBlock)
memberMat = NaN(numDist, numBlock);
for distInd = 1 : numDist
    memberMat(distInd, :) = memberVec((distInd - 1) * numBlock + 1 : distInd * numBlock); 
end
end

function compact = computeCompactness(x, y)
[~, dArea] = boundary(x, y, 1);
[~, chArea] = convhull(x, y);
compact = dArea / chArea;
end

function avgCompact = computeTotalCompactness(memberMat, numDist, blockX, blockY)
avgCompact = 0;
for distInd = 1 : numDist
    targetBlocks = memberMat(distInd, :) == 1;
    x = blockX(targetBlocks);
    y = blockY(targetBlocks);
    avgCompact = avgCompact + computeCompactness(x, y);
end
avgCompact = avgCompact / width(memberMat);
end

function [A, b] = buildConstraint(nDist, nBlock, demograph, tol)
%{
demograph(i,j) := block i, demography statistic j
%}

nDemo = width(demograph);
nvars = nDist * nBlock;

A = zeros(nDemo * nDist * 2, nvars);
b = zeros(nDemo * nDist * 2, 1);

counter = 1;
for demoId = 1 : nDemo
    for distId = 1 : nDist
        meanDemo = mean(demograph(:, demoId));
        for blockId = 1 : nBlock
            A(counter, (distId - 1 ) * nBlock + blockId) = demograph(blockId, demoId);
            A(counter + 1, (distId - 1 ) * nBlock + blockId) = -1 * demograph(blockId, demoId);
            b(counter) = meanDemo * (1 - tol);
            b(counter) = -1 * meanDemo * (tol - 1);
        end
        counter = counter + 2;
        counter / (nDemo * nDist * 2)
    end
end
end

function avgCompact = objFun(memberVec, numDist, numBlock, blockX, blockY)
memberMat = memberVec2Mat(memberVec, numDist, numBlock);
avgCompact = computeTotalCompactness(memberMat, numDist, blockX, blockY);
end

function [fun,nvar,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options] = initGA(numDist, numBlock, demograph, blockX, blockY, diffTol)
fun = @(memberVec) objFun(memberVec, numDist, numBlock, blockX, blockY);
nvar = numDist * numBlock;
[A, b] = buildConstraint(numDist, numBlock, demograph, diffTol);
beq = ones(numBlock, 1);
lb = zeros(nvar, 1);
ub = ones(nvar, 1);
nonlcon = [];
intcon = 1 : nvar;

Aeq = [];
for distInd = 1 : numDist
    Aeq = (horzcat(Aeq, eye(numBlock)));
    distInd / numDist
end

options=optimoptions('GA','Display','iter');
end

function [optim, obj] = runGA(fun,nvar,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options)
[optim, obj] = ga(fun,nvar,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
end

