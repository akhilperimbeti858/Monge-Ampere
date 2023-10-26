function [u,resMat,err] = GS(F,g,iterVec,h,u0,xa,xb,ya,yb,coarse)
% GS performs Gauss Seidel iterations

%Lambda is a constant used for calculating a weighted combination of the
%existing u and the updated u.
lambda = 0.05;

%Set the grid on which the exact solution g(x,y) will be applied.
[X,Y] = meshgrid(xa:h:xb,ya:h:yb);
G = g(X,Y);

%Initialize the matrix u.
u = u0;

%Set the RHS of the MA equation and the number of iterations 
%in the loop depending on input matrix u is on the coarse or fine level.
if coarse == 1
    maxCount = iterVec(1);
else
    maxCount = iterVec(2);
end

indices = zeros(size(u));
indices(2:end-1,2:end-1) = 1;

indices = find(indices == 1);

for cnt = 1:maxCount
    
    for k = 1:length(indices)
        
        [i,j] = ind2sub(size(u),indices(k));
        
        A_xy = (1/h)^4*(u(i-1,j)+u(i+1,j)-2*u(i,j))...
                      *(u(i,j-1)+u(i,j+1)-2*u(i,j));
          
        A_vw = 1/(4*h^4)*(u(i-1,j-1)+u(i+1,j+1)-2*u(i,j))...
                        *(u(i+1,j-1)+u(i-1,j+1)-2*u(i,j));
       
        if A_xy <= A_vw

            u(i,j) = 0.25*(u(i+1,j)+u(i-1,j)+u(i,j-1)+u(i,j+1))...
                -0.5*sqrt(0.25*((u(i+1,j)+u(i-1,j)-u(i,j-1)-u(i,j+1))^2)...
                +h^4*F(i,j));

        else

            u(i,j) = 0.25*(u(i-1,j+1)+u(i+1,j-1)+u(i+1,j+1)+u(i-1,j-1))...
                -0.5*sqrt(0.25*((u(i-1,j+1)+u(i+1,j-1)-u(i+1,j+1)-u(i-1,j-1))^2)...
                +4*h^4*F(i,j));

        end
        
    end
    
    uNew = u(2:end-1,2:end-1);
    
    %uNew = notJacobi(u,F,h,2);
    
    %Update u with a weighted sum of points from the u_old and
    %the newly calculated uNew.
    u(2:end-1,2:end-1) = lambda*u(2:end-1,2:end-1) + (1-lambda)*uNew;
    
    resMat = padarray(F(2:end-1,2:end-1) - MA_operator(u,h),[1,1],0);
    
end

%Subtract u from the exact solution to find the error matrix.
err = G - u;

%Pad the residual matrix with zeros to make it the right size again.
resMat = padarray(F(2:end-1,2:end-1) - MA_operator(u,h),[1,1],0);

end