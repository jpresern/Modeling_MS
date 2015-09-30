%%% Written by Aleš Škorjanc at some point in 2011

function [List] = Modeling_ListOfVariables(varargin)

List = [];

for a = 1 : 2 : length(varargin)
    names = varargin{a};
    values = varargin{a+1};
    for b = 1 : length(names)
        variable = {horzcat(names{b},' = ',num2str(values(b)))};
        List = vertcat(List,variable);
    end;
end;
