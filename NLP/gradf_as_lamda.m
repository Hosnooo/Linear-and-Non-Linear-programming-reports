function grad_lamda = gradf_as_lamda(f,s,x)
syms lamda real
n = length(x);
vars = sym('x', [1 n]);
grad_lamda(lamda) = subs(gradient(f),vars,(x + lamda*s)');
end