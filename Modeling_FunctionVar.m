%%% Written by AleS Skorjanc at some point in 2011

function FunVar = Modeling_FunctionVar(x,y)

for a = 1 : size(x,1)
    d = max(x(a,:)) - min(x(a,:));
    dx = d/1000;
    xi = min(x(a,:)):dx:max(x(a,:));
    yi = interp1(x(a,:),y(a,:),xi,'spline');
    FunMean = (sum(yi)*dx)/d;
    FunVar(a) = sum((yi - FunMean).^2)/1000;
end;