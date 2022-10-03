function [f_opt,x,func_vals,x_vals,g_vals,iteration] = Fletcher_Reeves(f,x0,e)
n = length(x0);
vars = sym('x', [1 n]);
iteration = 0;
grad = gradient(f,vars);
A = hessian(f,vars);
x = x0;
s_old = 0;
grad_old = 1;
conv = 100;
x_vals = [x];
g_vals = [];
%s_vals = [];
%betas = [];
%lamdas = [];
func_vals = [double(subs(f,vars,x'))];
while conv > e
    iteration = iteration + 1;
    if iteration == n + 1
        s_old = 0;
    end
    grad_i = double(subs(grad,vars,x'));
    A_i = double(subs(A,vars,x'));
    beta = (norm(grad_i)/norm(grad_old))^2;
    s = - grad_i + beta*s_old;
    lamda = -s'*grad_i/(s'*A_i*s);   
    x_old = x;
    s_old = s;
    grad_old = grad_i;
    g_vals = [g_vals norm(grad_i)];
    %s_vals = [s_vals s];
    %betas = [betas beta];
    %lamdas = [lamdas lamda];
    x = x_old + lamda*s;
    x_vals = [x_vals x];
    conv = norm(double(subs(grad,vars,x')));
    func_vals = [func_vals double(subs(f,vars,x'))];
end
f_opt = double(subs(f,vars,x'));
g_vals = [g_vals conv];
