function [x,x_iter,Taw_iter] = IP_central(A,b,c,eps,sigma,x,lamda,s, option)
%%
if strcmp(option, 'adaptive')
    [x, lamda, s] = Getstartingpoint(A,b,c);
end
d_gap = s'*x;
n = length(x);
x_iter = [];
s_iter = [];
Taw_iter = [];
k = 0;
while (d_gap > eps)
    %% Mew
    mu = d_gap/n;
    %% Taw
    taw = sigma*mu;
    Taw = x.*s;
    %% Deltas
    D = diag(x)./s;
    Y = inv(A*D*A');
    y = x - taw./s;
    rp = b - A*x;
    rd = c - A'*lamda - s;
    delta_lamda = Y*(A*y + A*D*rd + rp);
    delta_s = -A'*delta_lamda + rd;
    delta_x = -y-D*delta_s;
    %% Alpha
    xratio = -x./delta_x;
    xratio_pos = xratio(delta_x < 0);
    alpha_p = min(xratio_pos);
    sratio = -s./delta_s;
    sratio_pos = sratio(delta_s < 0);
    alpha_d = min(sratio_pos);
    alpha_max = min([alpha_p alpha_d]);
    alpha = (1-10^(-6))*alpha_max;
    %% Sigma
    if strcmp(option, 'adaptive')
        mu_aff = ((x + alpha_p * delta_x)' * (s + alpha_d * delta_s)) / n;
        sigma = (mu_aff / mu) ^ 3;
    end
    %% Next iteration parameters
    x = x + alpha*delta_x;
    lamda = lamda + alpha*delta_lamda;
    s = s + alpha*delta_s;
    %% Gap
    d_gap  = s'*x;
    k = k+1;
    x_iter = [x_iter x];
    s_iter = [s_iter s];
    Taw_iter = [Taw_iter Taw];
end
end