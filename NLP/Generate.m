%%
clear all; clc;
%% Functions
syms x1  x2 x3 x4
f = 100*(x2-x1^2)^2+(1-x1)^2;
x0 = [-1.2 1]';
e = 10^-3;
f2 = (x1 + 10*x2)^2 + 5*(x3 - x4)^2 + (x2 - 2*x3)^4 + 10*(x1 - x4)^4;
x02 = [3,-1,0,1]';
%%
tic;
[f_F,x_F,func_vals,x_vals,g_vals,i_F] = Fletcher_Reeves(f2,x02,e);
t_F = toc;
%%
tic;
[f_M,x_M,func_vals,x_vals,g_vals,i_M,c_M] = Marquardt(f2,x02,10^4,0.25,2,e);
t_M = toc;
%%
tic;
[f_QN_q,x_QN_q,func_vals,x_vals,g_vals,i_QN_q] = Quasi_Newton(f2,x02,e,0.001,0.1,'quad');
t_QN_q = toc;
tic;
[f_QN_c,x_QN_c,func_vals,x_vals,g_vals,i_QN_c] = Quasi_Newton(f2,x02,e,0.1,10^-3,'cubic');
t_QN_c = toc;
tic;
[f_QN_f,x_QN_f,func_vals,x_vals,g_vals,i_QN_f] = Quasi_Newton(f2,x02,e,0.001,10^-4,'fibonacci');
t_QN_f = toc;
tic;
[f_QN_g,x_QN_g,func_vals,x_vals,g_vals,i_QN_g] = Quasi_Newton(f2,x02,e,0.001,10^-4,'golden');
t_QN_g = toc;
%%
Names = {'Fletcher-Reeves', 'Marquardt', 'Quasi-Newton (Quadratic-interpolation)' ...
    'Quasi-Newton (Cubic-interpolation)', 'Quasi-Newton (Fibonacci)'...
    'Quasi-Newton (Golden-section)'}';
vars = {'Method','Noi','x1','x2','x3','x4','f','Cpu'};
is = [i_F i_M i_QN_q i_QN_c i_QN_f i_QN_g]';
xs = [x_F, x_M , x_QN_q, x_QN_c, x_QN_f, x_QN_g]';
fs = [f_F f_M f_QN_q f_QN_c f_QN_f f_QN_g]';
cpus = [t_F t_M t_QN_q t_QN_c t_QN_f t_QN_g]';
%%
Powel = table(Names,is,xs(:,1),xs(:,2),xs(:,3),xs(:,4),fs,cpus,'VariableNames',vars,'RowNames',Names);
writetable(Powel)