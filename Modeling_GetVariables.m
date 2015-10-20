%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output=Modeling_GetVariables(Project,ModelName,varargin)


inputDir = cd;
% in = load(horzcat(inputDir,'\',Project,'\',ModelName,'_Parameters.mat'));
% in_check = in;

% if an argument has been given and is not empty
if nargin > 2 && isempty(varargin{1}) == 0
    
    % go through all inputs and converts stringed values into numbers
    for st = 1: length(varargin);
        varargin(st) = {num2str(varargin{st})};
    end;
    
    % collect values and names of variables that will be fitted
    if strmatch('fit_',varargin);
        fit_var_ind = strmatch('fit_',varargin)';
        fit_var_count = 0;
        for r = 1 : length(fit_var_ind);
            fit_var_count = fit_var_count + 1;
            fit_var_name(fit_var_count) = {strrep(varargin{fit_var_ind(r)},'fit_','')};
            if fit_var_ind(r) < size(varargin,2);
                if isscalar(str2num(varargin{fit_var_ind(r)+1})); fit_var(fit_var_count) = str2num(varargin{fit_var_ind(r)+1});
                else fit_var(fit_var_count) = in.InitialValues.(strrep(varargin{fit_var_ind(r)},'fit_',''));
                end;
            else fit_var(fit_var_count) = in.InitialValues.(strrep(varargin{fit_var_ind(r)},'fit_',''));
            end;
        end;
    end;
    
    % collect tha values and the names of variables to fit
    if strmatch('const_',varargin);
        const_var_ind = strmatch('const_',varargin)';
        const_var_count = 0;
        for r = 1 : length(const_var_ind);
            const_var_count = const_var_count + 1;
            const_var_name(const_var_count) = {strrep(varargin{const_var_ind(r)},'const_','')};
            if const_var_ind(r) < size(varargin,2);
                if isscalar(str2num(varargin{const_var_ind(r)+1}));... % if you find a value for a constant variable change it in the InitialValue structure
                        in.InitialValues.(strrep(varargin{const_var_ind(r)},'const_','')) = str2num(varargin{const_var_ind(r)+1});
                end;
            end;
            const_var(const_var_count) = in.InitialValues.(strrep(varargin{const_var_ind(r)},'const_',''));
        end;
    end;
    
%     % if there are 
%     if isempty(strmatch('fit_',varargin));
%         names = fieldnames(in.InitialValues);
%         const_var_count = 0;
%         for r = 1 : length(names);
%             if isempty(strmatch(names{r},fit_var_name));
%                 const_var_count = const_var_count + 1;
%                 const_var(const_var_count) = in.InitialValues.(names{r});
%                 const_var_name(const_var_count) = {names{r}};
%             end;
%         end;
%     end;
%     
%     if isempty(strmatch('fit_',varargin));
%         names = fieldnames(in.InitialValues);
%         fit_var_count = 0;
%         for r = 1 : length(names);
%             if isempty(strmatch(names{r},const_var_name));
%                 fit_var_count = fit_var_count + 1;
%                 fit_var(fit_var_count) = in.InitialValues.(names{r});
%                 fit_var_name(fit_var_count) = {names{r}};
%             end;
%         end;
%     end;
    
else
    names = fieldnames(in.InitialValues);
    const_var_count = 0;
    
    for r = 1 : length(names)
        fit_var(r) = in.InitialValues.(names{r});
        fit_var_name(r) = {names{r}};
    end;
    
end;

% TrueVariables = fieldnames(in_check.InitialValues);
% 
% for vr = 1 : length(fit_var_name)
%     if isempty(strmatch(fit_var_name{vr},TrueVariables))
%         error(horzcat('Incorrect input variable ',fit_var_name{vr}));
%     end;
% end;
% 
% for vr = 1 : length(const_var_name)
%     if isempty(strmatch(const_var_name{vr},TrueVariables))
%         error(horzcat('Incorrect input variable ',const_var_name{vr}));
%     end;
% end;

output.fit_var = fit_var;
output.fit_var_name = fit_var_name;

if const_var_count ~= 0
    output.const_var = const_var;
    output.const_var_name = const_var_name;
else
    output.const_var = [];
    output.const_var_name = [];
end;
