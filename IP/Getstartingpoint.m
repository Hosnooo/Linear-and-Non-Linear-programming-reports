function [x, lamda, s] = Getstartingpoint(A,b,c)
Xtelda = A'*inv(A*A')*b;
w = inv(A*A')*A*c;
Stelda = c - A'*w;
delX = max(-1.5*min(Xtelda),0);
delS = max(-1.5*min(Stelda),0);
delXdash =  delX + 0.5*(Xtelda + delX*ones(size(Stelda,1),1))'*...
    (Stelda + delS*ones(size(Xtelda,1),1))/(sum(Stelda + delS));

delSdash =  delS + 0.5*(Xtelda + delX*ones(size(Stelda,1),1))'*...
    (Stelda + delS*ones(size(Xtelda,1),1))/(sum(Xtelda + delX));
x = Xtelda + delXdash;
lamda = w;
s = Stelda + delSdash;
end