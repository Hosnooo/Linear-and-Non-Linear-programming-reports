function [lamda,f_opt] = Fibonacci(f_lamda,range,e)
a = range(1); b = range(2); lim = (b - a)/e; F = -inf; n = 0;
while F < lim
    n  = n + 1; F = fibonacci(n);
end
k = 0;
L2 = fibonacci(n - 1)/F*(b - a);
while k < n - 3
    L1 = (b - a);
    F_n1 = fibonacci(n - k - 1); F = fibonacci(n - k); r = F_n1 / F;
    if L2 > L1/2
        x1 = b - r*(b - a);
        x2 = a + r*(b - a);
    else
        x2 = b - r*(b - a);
        x1 = a + r*(b - a);
    end
    f1 = double(f_lamda(x1)); f2 = double(f_lamda(x2));
    if f1 > f2
        a = x1;
        L2 = r*L1;
    elseif f2 > f1 
        b = x2;
        L2 = r*L1;
    else
        a = x1; b = x2;
        L2 = r*(b - a);
        k = k + 1;
    end
    k = k + 1;
end
[f_opt,i] = min([f1 f2]);
lamdas = [x1 x2]; lamda = lamdas(i);