function f_lamda = f_as_lamda(f,s,x)
syms lamda real
n = length(x);
vars = sym('x', [1 n]);
f_lamda(lamda) = subs(f,vars,(x + lamda*s)');
end