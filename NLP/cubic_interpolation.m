function [lamda_opt,f_opt] = cubic_interpolation(f_lamda,grad_lamda,s,t0,e)
[lamda_opt,f_opt] = cubic_iterate(f_lamda,0,t0,s,grad_lamda,e); 
end

function [A,B,fA,fB,fdash_A,fdash_B]= iterate(f,A,B)
dif = diff(f);
fdash_B = double(dif(B));
if fdash_B < 0
    A = B;
    B = 2*B;
    [A,B,fA,fB,fdash_A,fdash_B] = iterate(f,A,B);
else
    fdash_A = double(dif(A));
    fA = double(f(A));
    fB = double(f(B));
end
end

function [lamda_opt,f_opt] = cubic_iterate(f_lamda,A,B,s,grad_lamda,e)
[A,B,fA,fB,fdash_A,fdash_B]= iterate(f_lamda,A,B);
Z = 3*(fA - fB)/(B - A) + fdash_A + fdash_B;
Q = (Z^2 - fdash_A*fdash_B)^0.5;
lamda_1 = A + (fdash_A + Z + Q)/(fdash_A + fdash_B + 2*Z)*(B - A);
lamda_2 = A + (fdash_A + Z - Q)/(fdash_A + fdash_B + 2*Z)*(B - A);
lamdas = [lamda_1 lamda_2]; cond = lamdas > A & lamdas < B; lamda_opt = lamdas(cond); 
grad_opt = double(subs(grad_lamda,lamda_opt));
f_opt = double(f_lamda(lamda_opt));
conv = abs(s'*grad_opt);
if conv > e
    if f_opt > 0
        B = lamda_opt;
    else
        A = lamda_opt;
    end
    lamda_opt = cubic_iterate(f_lamda,A,B,s,grad_lamda,e);
end
end
