function [lamda,f_opt] = golden_section(f_lamda,range,e)
a = range(1); b = range(2);
d = 0.618*(b - a);
while d > e
    x1 = a + d;
    x2 = b - d;
    f1 = double(f_lamda(x1));
    f2 = double(f_lamda(x2));
    if f1 > f2
        b = x1;
    elseif f2 > f1
        a = x2;
    else
        b = x1;
        a = x2;
    end
    d = 0.618*(b - a);   
end
[f_opt,i] = min([f1 f2]);
lamdas = [x1 x2]; lamda = lamdas(i);
end