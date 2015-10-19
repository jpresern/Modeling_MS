function Modeling_DRG_loop(Project_in)


global Project; 
Project = Project_in;
%Project = 'DRG_TCM';
global inputDir; inputDir = cd;

%   Loading input parameters and experimental data into the function 
%   Function executes loop in which it changes the starting parameters to
%   new settings. Suitable for overnight runs.

%   There are five adjustable parameters.
%   C_Weights      ...  which comparisons should be done and how weighted
%                       they are
%   FitVariables   ...  starting values of the fitted variables
%   FitVariables_limits ... limits whithin wich the FitVariables are
%                       rotated
%   FitVariables_tol... tolerance determining how much over the limit can
%                       the FitVariables be pushed

%%%%%%%%%%%%%% IMPORTANT

%   Manipulated variables should have the same number of columns!

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% START Load data files %%%%%%%%%%%%%%

%   opens variable but does not load it
InitFitParamFile = matfile(filename(inputDir,{inputDir,Project},[Project,'_InitialFitParameters.mat']),'Writable',true); 
%   loads components of the opened variable
variableFitParams = InitFitParamFile.variableFitParams;
variableFitParams_HiLims = InitFitParamFile.variableFitParams_HiLims;
variableFitParams_LoLims = InitFitParamFile.variableFitParams_LoLims;
variableFitParams_tol = InitFitParamFile.variableFitParams_tol;
variableFitC_Weights = InitFitParamFile.variableFitC_Weights;


for i = 1:size(variableFitParams,2)
   %    writes new values into appropriate variables
   ModelInput_FitVariables = variableFitParams (:,i)';
   ModelInput_FitVariables_limits(:,1) = variableFitParams_LoLims (:,i); 
   ModelInput_FitVariables_limits(:,2) = variableFitParams_HiLims (:,i);
   ModelInput_FitVariables_tol = variableFitParams_tol(:,i)';
   C_Weights = variableFitC_Weights(:,i)';
   
   %    writes the components of the opened variable back to the file
   InitFitParamFile.ModelInput_FitVariables = ModelInput_FitVariables;
   InitFitParamFile.ModelInput_FitVariables_limits = ModelInput_FitVariables_limits;
   InitFitParamFile.ModelInput_FitVariables_tol = ModelInput_FitVariables_tol;
   InitFitParamFile.C_Weights = C_Weights;
   
   %    starts the model with new parameters
   Modeling_DRG(Project);
   % pause;
end

end