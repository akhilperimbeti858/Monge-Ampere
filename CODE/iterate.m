function [u,resRec,errMat,time,count] = iterate(F,g,n,N,levels,iterVec,h,u0,xa,xb,ya,yb,tol)
tic
count = 0;
u = u0;
res = 1;

while res > tol
    count = count + 1;
    [u,resMat,err] = FAS(F,g,n,N,levels,iterVec,h,u,xa,xb,ya,yb,count);
    res = norm(resMat(:),inf); 
    if count == 1
        errMat = err;
        resRec = resMat;
    else
        errMat = cat(3,errMat,err);
        resRec = cat(3,resRec,resMat);
    end
    
end

time = toc;

end