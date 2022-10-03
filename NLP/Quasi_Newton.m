function [f_opt,x,func_vals,x_vals,g_vals,iteration] = Quasi_Newton(f,x,e,t0,e_in,oneD)
n = length(x);
vars = sym('x', [1 n]);
iteration = 0;
grad = gradient(f,vars);
B = eye(n);
grad_old = double(subs(grad,vars,x'));
g_vals = [norm(grad_old)];
x_vals = [x];
func_vals = [double(subs(f,vars,x'))];
iteration = 0;
conv = 100;
while conv > e
    iteration = iteration + 1;
    s = -B*grad_old;
    f_lamda = f_as_lamda(f,s,x);
    grad_lamda = gradf_as_lamda(f,s,x);
    if strcmp(oneD,'cubic')
        [lamda,f_new] = cubic_interpolation(f_lamda,grad_lamda,s,t0,e_in);
    elseif strcmp(oneD,'quad')
        [lamda,f_new] = quadratic_interpolation(f_lamda,t0,e_in);
    elseif strcmp(oneD,'fibonacci')
        [lamda,f_new] = Fibonacci(f_lamda,[0 1],e_in);
    elseif strcmp(oneD,'golden')
        [lamda,f_new] = golden_section(f_lamda,[0 1],e_in);
    end
    func_vals = [func_vals f_new];
    x = x + lamda*s;
    x_vals = [x_vals x];
    grad_new = double(subs(grad,vars,x'));
    conv = norm(grad_new);
    g_vals = [g_vals conv];
    d = lamda*s;
    g = grad_new - grad_old;
    B = B + ((1 + g'*B*g/(d'*g))*(d*d') - d*g'*B - B*g*d')/(d'*g);
    grad_old = grad_new;
end
f_opt = f_new;
end