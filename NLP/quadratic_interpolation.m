function [lamda_opt,f_opt] = quadratic_interpolation(f_lamda,t0,e)
A = 0; B = t0; C = 2*B;
conv = 100;
i = 0;
while conv > e
    i = i + 1;
    fA = double(f_lamda(A));
    f1 = double(f_lamda(B));
    if f1 > fA
        fC = f1;
        C = B;
        B = B/2;
        fB = double(f_lamda(B));
    else
        t_1 = B;
        t_2 = C;
        f2 = double(f_lamda(t_2));
        while f2 < f1
            f1 = f2;
            t_1 = t_2;
            t_2 = 2*t_2;
            f2 = double(f_lamda(t_2));
        end
        B = t_1;
        C = t_2;
        fB = f1;
        fC = f2;
    end
    lamda_opt = (4*fB -3*fA - fC)/(4*fB - 2*fC - 2*fA)*B;
    if i ~=1 && (lamda_old - lamda_opt) < 10^-10
        break;
    end
    lamda_old = lamda_opt;
    a = fA;
    b = (4*fB - 3*fA - fC)/(2*B);
    c = (fC + fA - 2*fB)/(2^B^2);
    h = a + b*lamda_opt + c*lamda_opt^2;
    f_opt = double(f_lamda(lamda_opt));
    conv = abs((h - f_opt)/f_opt);
    if lamda_opt > B
        if f_opt < fB
            A = B;
            B = lamda_opt;
        else
            C = lamda_opt;
        end
    else
        if f_opt  < fB
            B = lamda_opt;
            C = B;
        else
            A = lamda_opt;
        end
    end
end