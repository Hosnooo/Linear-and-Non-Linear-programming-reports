function builinresults()
disp('These are MATLAB built in function Linprog interior-point results')
%% Built in Linprog
%% Case 1
objfunc = [1 0 -1 -3];                                                                                                   
Aeq = [2 0 2 3;0 -2 -2 -6];                                                 
beq = [10 -6]';                                                              
A = [];                                                                     
b = [];                                                                     
lb = [0 0 0 0];
ub = [];
opt = optimoptions('linprog','Algorithm','interior-point');
X = linprog(objfunc,A,b,Aeq,beq,lb,ub,opt);
disp('case 1: x = ')
disp(X)
%% Case 2
objfunc = [2 9 3];                                                                                                  
Aeq = [];                                                 
beq = [];                                                              
A = -[-2 2 1 ;1 1 -1];                                                   
b = -[1 1];                                                                  
lb = [0 0 0];
ub = [];
opt = optimoptions('linprog','Algorithm','interior-point');
X = linprog(objfunc,A,b,Aeq,beq,lb,ub,opt);
disp('case 2: x = ')
disp(X)
%% Case 3
objfunc = [5 2 -4];                                                                                                
Aeq = [];                                                 
beq = [];                                                              
A = [-6 -1 2 ;1 1 1 ;-6 -4 2];                              
b = [-5 4 -10];                                                             
lb = [0 0 0];
ub = [];
opt = optimoptions('linprog','Algorithm','interior-point');
X = linprog(objfunc,A,b,Aeq,beq,lb,ub,opt);
disp('case 3: x = ')
disp(X)
end