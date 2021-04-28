function eta = water_filling(t,P)
    n = length(t);
    eta_max = max(t);
    epsilon = 1e-3;
    if eta_max * n - sum(t) <= P
        eta = eta_max + (P- (eta_max * n - sum(t))) / n;
    else
        t = sort(t);
        l = t(1);
        r = eta_max;
        while 1
            eta = (l + r) / 2;
            for i = 1:n
                if eta <= t(i)
                    break;
                end
            end
            vol = (i-1) * eta - sum(t(1:i-1));
            if abs(vol- P) <= epsilon
                break;
            elseif vol < P
                l = eta;
            else
                r = eta;
            end
        end
    end
end
