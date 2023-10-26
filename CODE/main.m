close all; clear all;

% Main script implementation 

f = @(x,y) (1 + x.^2 + y.^2).*exp(x.^2 + y.^2);
g = @(x,y) exp((x.^2 + y.^2)/2);


f1 = @(x,y) (1 + x.^2).*(1 + y.^2).*exp(x.^2 + y.^2);
f2 = @(x,y) exp(x.^2 + y.^2).*(1+.5*(x+y).^2).*(1+.5*(y-x).^2);

% setting min and max values for solution matrix
minN = 4;
maxN = 7;
nValNum = length(minN:maxN);

stats = zeros(nValNum,maxN,2);
nVec = 2.^(minN+1:maxN+1)+1;
hVec = 2.^(-(minN+1:maxN+1));

errCell = cell(nValNum,maxN);
resCell = cell(nValNum,maxN);
frameCell = cell(nValNum,maxN);


exactSolCell = cell(nValNum,1);
numSolCell = cell(nValNum,1);

basesNum = 2; % Number of orthogonal bases - for 3x3 wide-stencil
iterVec = [5 50];

for i = minN:maxN
    n = 2^(i+1) + 1;
    h = 2/(n-1);
    xa = -1; xb = 1; ya = -1; yb = 1; tol = h^2/10;
    [X,Y] = meshgrid(xa:h:xb,ya:h:yb);

    F1 = f1(X,Y);
    F2 = f2(X,Y);
    F = min(F1,F2); %minimizing function
    

    G = g(X,Y);
    exactSolCell{i-minN+1} = G; %saving exact solution 
    A = zeros(n);
    u0 = initialize(F,g,n,h,X,Y);

    u0(:,1) = g(X(:,1),Y(:,1));
    u0(:,n) = g(X(:,n),Y(:,n));
    u0(1,:) = g(X(1,:),Y(1,:));
    u0(n,:) = g(X(n,:),Y(n,:)); 

    N = n;
    
    % four levels/depths of recursion
    for j = 1:4 
        [u,resMat,errMat,time,count] = iterate(F,g,n,N,j,2*iterVec,h,u0,xa,xb,ya,yb,tol);
        stats(i-minN+1,j,1) = norm(errMat(:,:,end),2); % storing error
        stats(i-minN+1,j,2) = count;    %storing iteration count
        errFrames = unique(ceil(linspace(1,count,4)));
        frameCell{i-minN+1,j} = errFrames;
        errCell{i-minN+1,j} = errMat(:,:,errFrames);
        resCell{i-minN+1,j} = resMat(:,:,errFrames);

    end
    
    numSolCell{i-minN+1} = u; %saving numerical solution
end
   
%% Error and iteration number plots

legendStrs = {'One level','Two levels','Three levels','Four levels'};
err_fig = figure;
plot(hVec,stats(:,1,1), 'o-')
xlabel('h')
ylabel('Error')
title('Abs. error vs. h (grid)','all depths of recursion')
axis tight
grid on
saveas(err_fig,'err_vs_h.fig')

iter_fig = figure;
semilogy(hVec,stats(:,:,2),'o-');
legend(legendStrs(1:4));
xlabel('h')
ylabel('Iterations')
title('Number of iterations vs. h', sprintf('for %d levels of recursion',4))
axis tight
grid on
saveas(iter_fig,'total_iters.fig')

%% Plotting - Exact solution, Numerical Solution, Error and Residuals

errorDir = 'error_surfs';
mkdir(errorDir);
resDir = 'res_surfs';
mkdir(resDir);
exactSolDir = 'exact_sol_surfs';
mkdir(exactSolDir);
numSolDir = 'numerical_sol_surfs';
mkdir(numSolDir);

for i = 1:nValNum
    
    nValDir1 = sprintf('%s/N_%d',errorDir,1/hVec(i)+1);
    mkdir(nValDir1)
    
    nValDir2 = sprintf('%s/N_%d',resDir,1/hVec(i)+1);
    mkdir(nValDir2)

    for j = 1:maxN
        if ~isempty(errCell{i,j}) 
            depthDir1 = sprintf('%s/depth_%d',nValDir1,j);
            mkdir(depthDir1)

            depthDir2 = sprintf('%s/depth_%d',nValDir2,j);
            mkdir(depthDir2)

            subplotNum = size(errCell{i,j},3);
            err = errCell{i,j};
            res = resCell{i,j};
            
            for k = 1:subplotNum

                fig1 = figure;
                surf(linspace(-1,1,nVec(i)),linspace(-1,1,nVec(i)),abs(err(:,:,k)),...
                    'linestyle','none');
                title(sprintf('h = %f, depth = %d levels, iteration = %d',...
                    hVec(i),j,frameCell{i,j}(k)), 'Error');
                xlabel('X'); ylabel('Y'); zlabel('error');
                saveas(fig1,sprintf('%s/count_%d.fig',depthDir1,frameCell{i,j}(k)));

                fig2 = figure;
                surf(linspace(-1,1,nVec(i)),linspace(-1,1,nVec(i)),res(:,:,k),...
                    'linestyle','none');
                title(sprintf('h = %f, depth = %d levels, iteration = %d',...
                    hVec(i),j,frameCell{i,j}(k)), 'Residual Error');
                xlabel('X'); ylabel('Y'); zlabel('residual');
                saveas(fig2,sprintf('%s/count_%d.fig',depthDir2,frameCell{i,j}(k))); 

            end
        end
    end
    
    fig = figure;
    surf(exactSolCell{i},'linestyle','none');
    xlabel('X'); ylabel('Y'); zlabel('u');
    title('Exact solution',sprintf('Evaluated for h = %f',hVec(i)));
    saveas(fig,sprintf('%s/N_%d.fig',exactSolDir,1/hVec(i)+1));
    
    fig2 = figure;
    surf(numSolCell{i},'linestyle','none');
    xlabel('X'); ylabel('Y'); zlabel('v');
    title('Numerical solution',sprintf('Evaluated for h = %f',hVec(i)));
    saveas(fig2,sprintf('%s/N_%d.fig',numSolDir,1/hVec(i)+1));
end