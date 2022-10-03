function [f_opt,x,func_vals,x_vals,g_vals,iteration,costs] = Marquardt(f,x0,alpha0,c1,c2,e)
n = length(x0);
vars = sym('x', [1 n]);
iteration = 0;
grad = gradient(f,vars);
A = hessian(f,vars);
x = x0;
alpha = alpha0;
conv = 100;
x_vals = [x];
g_vals = [];
costs = [];
func_vals = [double(subs(f,vars,x'))];
while conv > e
    iteration = iteration + 1;
    grad_i = double(subs(grad,vars,x'));
    g_vals = [g_vals norm(grad_i)];
    A_i = double(subs(A,vars,x'));
    f_old = double(subs(f,vars,x'));
    [x,alpha,s,f_old,cost] = get_x(f,f_old,A_i,alpha,c1,c2,grad_i,x,n,vars,0);
    costs = [costs cost];
    x_vals = [x_vals x];
    func_vals = [func_vals f_old];
    conv = norm(double(subs(grad,vars,x')));
end
f_opt = f_old;
g_vals = [g_vals conv];
end

function [x,alpha,s,f_old,cost] = get_x(f,f_old,A_i,alpha,c1,c2,grad_i,x,n,vars,cost)
    s = -inv(A_i + alpha*eye(n))*grad_i;
    x_new = x + s;
    f_new = double(subs(f,vars,x_new'));
    if f_new < f_old
        alpha = c1*alpha;
        x = x_new;
        f_old = f_new;
    else
        cost = cost + 1;
        alpha = c2*alpha;
        [x,alpha,s,f_old,cost] = get_x(f,f_old,A_i,alpha,c1,c2,grad_i,x,n,vars,cost);
    end
end
