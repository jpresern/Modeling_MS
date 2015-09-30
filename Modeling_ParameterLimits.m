%%% Written by Aleš Škorjanc at some point in 2011

function r2 = Modeling_ParameterLimits(param,lim1,lim2,tol)

if param > lim1 & param < lim2
    r2 = 0;
end;

if param <= lim1
    if isinf(lim2)
        dif = param-lim1;
        rel_dif = dif/lim1;
    else
        dif = lim1 - param;
        rel_dif = dif/(lim2 - lim1);
    end;
    over_tol = rel_dif/tol;
    r2 = over_tol^2;
end;

if param >= lim2
    if isinf(lim1)
        dif = param-lim2;
        rel_dif = dif/lim2;
    else
        dif = param - lim2;
        rel_dif = dif/(lim2 - lim1);
    end;
    over_tol = rel_dif/tol;
    r2 = over_tol^2;
end;
