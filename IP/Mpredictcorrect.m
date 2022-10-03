function [x,x_iter,s_iter] = Mpredictcorrect(A,b,c,eps)
%%
[x, w, s] = Getstartingpoint(A,b,c);
n = length(x);
x_iter = [];
s_iter = [];
k = 0;
rd = A'*w + s - c;
rp = A*x - b;
rc = diag(x)*s;
mu = x'*s/n;                                                       % duality measure

while (max([mu norm(rp) norm(rd)]) > eps)
    %% corrector
    % M & RHS
    M = A*diag(x)*inv(diag(s))*A';
    rhs = rp + A*inv(diag(s))*(-rc + diag(x)*rd);
    % Delta
    U = chol(M);
    d = U'\rhs;
    delta_w = chol(M)\d;
    delta_s = rd - A'*delta_w;
    delta_x = inv(diag(s))*(rc - diag(x)*delta_s);
    % Alpha
    xratio = x./delta_x;
    xratio_pos = xratio(delta_x > 0);
    alpha_p = min([1 xratio_pos']);
    sratio = s./delta_s;
    sratio_pos = sratio(delta_s > 0);
    alpha_d = min([1 sratio_pos']);
    % Sigma
    mu_aff = ((x - alpha_p * delta_x)' * (s - alpha_d * delta_s)) / n;
    sigma = (mu_aff / (mu)) ^ 3;
    %% predictor
    rc = rc - sigma*mu + delta_x.*delta_s;
    rhs = rp + A*inv(diag(s))*(-rc + diag(x)*rd);
    % Delta
    U = chol(M);
    d = U'\rhs;
    delta_w = chol(M)\d;
    delta_s = rd - A'*delta_w;
    delta_x = inv(diag(s))*(rc - diag(x)*delta_s);
    % Alpha
    xratio = x./delta_x;
    xratio_pos = xratio(delta_x > 0);
    eta = max([0.995 1-mu]);
    alpha_p = min([1 eta*xratio_pos']);
    sratio = s./delta_s;
    sratio_pos = sratio(delta_s > 0);
    alpha_d = min([1 eta*sratio_pos']);
    %% Next iteration parameters
    x = x - alpha_p*delta_x;
    w = w - alpha_d*delta_w;
    s = s - alpha_d*delta_s;
    mu = x'*s/n;
    rd = A'*w + s - c;
    rp = A*x - b;
    rc = diag(x)*s;
    %% Memory
    k = k+1;
    x_iter = [x_iter x];
    s_iter = [s_iter s];
end
end