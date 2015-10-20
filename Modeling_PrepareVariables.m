%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InputVariables = Modeling_PrepareVariables(fit_var_names,fit_var,const_var_names,const_var)

for a = 1:length(fit_var);
    FitVariable(a*2-1) = {horzcat('fit_',fit_var_names{a})};
    FitVariable(a*2) = {fit_var(a)};
end;

for b = 1:length(const_var);
    ConstVariable(b*2-1) = {horzcat('const_',const_var_names{b})};
    ConstVariable(b*2) = {const_var(b)};
end;

if isempty(fit_var)& isempty(const_var) == 0;
    InputVariables = ConstVariable;
elseif isempty(const_var) & isempty(fit_var) == 0;
    InputVariables = FitVariable;
elseif isempty(const_var) == 0 & isempty(fit_var) == 0;
    InputVariables = horzcat(FitVariable,ConstVariable);
else
    InputVariables = [];
end;


