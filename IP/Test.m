close; clc; clear all;
%% save
save  = 0; % please change this to zero when you are running the code, this saves the figures as pngs
%%
eps = 10^-8;
%% Case 1
%min z = x1 - x3 - 3*x4
%s.t. 2*x1 + 2*x3 + 3*x4 = 10
%  -2*x2 - 2*x3 - 6*x4 = -6
A = [2 0 2 3;0 -2 -2 -6];                                                   
c = [1 0 -1 -3]';
b = [10 -6]';  
% 
x0 = [1 1 1 1]';
lamda0 = [1 1]';
s0 = [1 1 1 1]';
[x,x_iter,taw_iter] = IP_central(A,b,c,eps,0.2,x0,lamda0,s0,'constant');
iteration = 1:size(x_iter,2);
f = x_iter(1,:) - x_iter(3,:) - 3*x_iter(4,:);
[f1,f2,f3] = plotter('Case study: 1 (constant)',x_iter,f,iteration,taw_iter,'Taw');
if save == 1
    saveas(f1,'path1constant.png')
    saveas(f2,'function1constant.png')
    saveas(f3,'condition1constant.png')
end
disp('Case 1 results: x =')
disp(x)
%
[x,x_iter,taw_iter] = IP_central(A,b,c,eps,0.2,x0,lamda0,s0,'adaptive');
iteration = 1:size(x_iter,2);
f = x_iter(1,:) - x_iter(3,:) - 3*x_iter(4,:);
[f1,f2,f3] = plotter('Case study: 1 (adpative)',x_iter,f,iteration,taw_iter,'Taw');
if save == 1
    saveas(f1,'path1adaptive.png')
    saveas(f2,'function1adaptive.png')
    saveas(f3,'condition1adaptive.png')
end
disp(x)
%
[x,x_iter,s_iter] = Mpredictcorrect(A,b,c,eps);
iteration = 1:size(x_iter,2);
f = x_iter(1,:) - x_iter(3,:) - 3*x_iter(4,:);
condition = [];
for i = 1:length(iteration)
    condition = [condition  x_iter(:,i).*s_iter(:,i)];
end
[f1,f2,f3] = plotter('Case study: 1 (Mehrotra)',x_iter,f,iteration,condition,'Mehrotra comp. condition');
if save == 1
    saveas(f1,'path1mehr.png')
    saveas(f2,'function1mehr.png')
    saveas(f3,'condition1mehr.png')
end
disp(x)
%% Case 2
%min z = 2*x1 + 9*x2 + 3*x3
%s.t. -2*x1 + 2*x2 + x3 >= 1
%  x1 + 4*x2 - x3 >= 1
A = [-2 2 1 -1 0;1 1 -1 0 -1];                                                   
c = [2 9 3 0 0]';
b = [1 1]';  
% 
x0 = [1 1 1 1 1]';
lamda0 = [1 1]';
s0 = [1 1 1 1 1]';
[x,x_iter,taw_iter] = IP_central(A,b,c,eps,0.2,x0,lamda0,s0,'constant');
iteration = 1:size(x_iter,2);
f = 2*x_iter(1,:) + 9*x_iter(2,:) + 3*x_iter(3,:);
[f1,f2,f3] = plotter('Case study: 2 (constant)',x_iter,f,iteration,taw_iter,'Taw');
if save == 1
    saveas(f1,'path2constant.png')
    saveas(f2,'function2constant.png')
    saveas(f3,'condition2constant.png')
end
disp('Case 2 results: x =')
disp(x)
%
[x,x_iter,taw_iter] = IP_central(A,b,c,eps,0.2,x0,lamda0,s0,'adaptive');
iteration = 1:size(x_iter,2);
f = 2*x_iter(1,:) + 9*x_iter(2,:) + 3*x_iter(3,:);
[f1,f2,f3] = plotter('Case study: 2 (adpative)',x_iter,f,iteration,taw_iter,'Taw');
if save == 1
    saveas(f1,'path2adaptive.png')
    saveas(f2,'function2adaptive.png')
    saveas(f3,'condition2adaptive.png')
end
disp(x)
%
[x,x_iter,s_iter] = Mpredictcorrect(A,b,c,eps);
iteration = 1:size(x_iter,2);
f = 2*x_iter(1,:) + 9*x_iter(2,:) + 3*x_iter(3,:);
condition = [];
for i = 1:length(iteration)
    condition = [condition  x_iter(:,i).*s_iter(:,i)];
end
[f1,f2,f3] = plotter('Case study: 2 (Mehrotra)',x_iter,f,iteration,condition,'Mehrotra comp. condition');
if save == 1
    saveas(f1,'path2mehr.png')
    saveas(f2,'function2mehr.png')
    saveas(f3,'condition2mehr.png')
end
disp(x)
%% Case 3
% min z = 5*x1 + 2*x2 - 4*x3
% st.     6*x1 +  x2  - 2*x3 >= 5
%           x1 +  x2      x3 <= 4
%         6*x1 + 4x2  - 2*x3 >= 10
A = [6 1 -2 -1 0 0;1 1 1 0 1 0;6 4 -2 0 0 -1];                              
c = [5 2 -4 0 0 0]';
b = [5 4 10]';          
% 
x0 = [1 1 1 1 1 1]';
lamda0 = [1 1 1]';
s0 = [1 1 1 1 1 1]';
[x,x_iter,taw_iter] = IP_central(A,b,c,eps,0.2,x0,lamda0,s0,'constant');
iteration = 1:size(x_iter,2);
f = x_iter(1,:) - x_iter(3,:) - 3*x_iter(4,:);
[f1,f2,f3] = plotter('Case study: 3 (constant)',x_iter,f,iteration,taw_iter,'Taw');
if save == 1
    saveas(f1,'path3constant.png')
    saveas(f2,'function3constant.png')
    saveas(f3,'condition3constant.png')
end
disp('Case 3 results: x =')
disp(x)
%
[x,x_iter,taw_iter] = IP_central(A,b,c,eps,0.2,x0,lamda0,s0,'adaptive');
iteration = 1:size(x_iter,2);
f = x_iter(1,:) - x_iter(3,:) - 3*x_iter(4,:);
[f1,f2,f3] = plotter('Case study: 3 (adpative)',x_iter,f,iteration,taw_iter,'Taw');
if save == 1
    saveas(f1,'path3adaptive.png')
    saveas(f2,'function3adaptive.png')
    saveas(f3,'condition3adaptive.png')
end
disp(x)
%
[x,x_iter,s_iter] = Mpredictcorrect(A,b,c,eps);
iteration = 1:size(x_iter,2);
f = x_iter(1,:) - x_iter(3,:) - 3*x_iter(4,:);
condition = [];
for i = 1:length(iteration)
    condition = [condition  x_iter(:,i).*s_iter(:,i)];
end
[f1,f2,f3] = plotter('Case study: 3 (Mehrotra)',x_iter,f,iteration,condition,'Mehrotra comp. condition');
if save == 1
    saveas(f1,'path3mehr.png')
    saveas(f2,'function3mehr.png')
    saveas(f3,'condition3mehr.png')
end
disp(x)
%% Built in results
builinresults()
%% Plotter
function [f1,f2,f3] = plotter(name,x_iter,f,iteration,condition,condition_name)
f1 = figure('Name',name);
plot3(x_iter(1,:),x_iter(2,:),x_iter(3,:),'-*','LineWidth',1.5)
xlabel('x1')
ylabel('x2')
zlabel('x3')
title('Central Path')
grid on
f2 = figure('Name',name);
plot(iteration,f,'-o','LineWidth',1.5)
xlabel('iteration')
ylabel('Objective function')
title('Objective function vs iteration')
grid on
f3 = figure('Name',name);
plot(iteration,condition,'-o','LineWidth',1.5)
xlabel('iteration')
ylabel(condition_name)
title('Condition vs iteration')
grid on
end
