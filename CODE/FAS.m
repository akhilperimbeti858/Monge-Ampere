function [u,resMat,err] = FAS(F,g,n,N,levels,iterVec,h,u0,xa,xb,ya,yb,count)
%FAS_V2.m implements the Full Approximation Scheme. 

%Do the initial Gauss-Seidel iteration.
[u,resMat] = GS(F,g,iterVec,h,u0,xa,xb,ya,yb,0);
res = norm(resMat(:),inf);


%Go down one level finer mesh, setting a new n and h.
n = (n+1)/2;
h = 2*h;

%Set the coarse versions of u, F, and the residual matrix.
v = restrict(u);
resCoarse = restrict(resMat);

%Evaluate the MA equation on the coarse grid, pad it with zeros so that
%it's the same size as resCoarse, and add them together.
A = padarray(MA_operator(v,h),[1,1],0) + resCoarse;

%run the GS calculation if lowest level and set
%vNew as the output. Otherwise, set vNew as the recursive output of FAS.
if floor(log2(N)) - floor(log2(n)) ~= levels

    vNew = FAS(A,g,n,N,levels,iterVec,h,v,xa,xb,ya,yb,count);
    
else
    
    vNew = GS(A,g,iterVec,h,v,xa,xb,ya,yb,1);
    
end


%Calculate the coarse error matrix.
eCoarse = vNew - v;

%Interpolate the coarse error matrix to the fine level.
eFine = interpolate(eCoarse);

%Add the error matrix to the original matrix u.
u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + eFine(2:end-1,2:end-1);

%Come back up to the fine level.
n = 2*n-1;
h = h/2;

%Do a few more Gauss-Seidel iterations on u and calculate the residual 
%matrix.
[u,resMat] = GS(F,g,iterVec,h,u,xa,xb,ya,yb,0);
res = norm(resMat(:),inf);


%If we're on the finest level, calculate the error. Otherwise, set it to
%zero so as to avoid a 'not enough output arguments' error.
if n == N
    [X,Y] = meshgrid(xa:h:xb,ya:h:yb);
    G = g(X,Y);
    err = G-u;
else
    err = 0;
end

end