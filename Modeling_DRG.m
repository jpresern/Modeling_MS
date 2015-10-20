%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function fits HH equations on mechanically activated currents in
%%% various neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.One starts the function with entering "Modeling_DRG('DRG')" at the command window
%   workspace MUST be set to folder containing Modeling_DRG.m
% 1.Folder containing script must contain experimental subfolder (in this 
%   example case 'DRG') in which the two files reside (see example):
%       -file with start-up parameters
%       -file with the initial fit parameters
% 2.The experimental subfolder must also contain folder 'Results' in which
%   the program dumps the results, figures etc....
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Modeling_DRG(Project_in)
%%

global Project; 
Project = Project_in;
global inputDir; inputDir = cd;

%   Loading input parameters and experimental data into the function 
%       the default parameters file is: DRG_TCM_InitialFitParametersTemp.mat
%       the default experimental data file is:
%       DRG_TCM_ExperimentalData_Hao2010.mat

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% START Load data files %%%%%%%%%%%%%%

%   Custom written function 'filename' is used to avoid problems using
%   scripts on different OS platforms
InitialFitParam = load(filename(inputDir,{inputDir,Project},[Project,'_InitialFitParameters.mat']));
ExpData = load(filename(inputDir,{inputDir,Project},[InitialFitParam.ExperimentalData,'.mat']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% END   Load data files %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% START prepare parameters loaded from file %%%%%%%%%%%%%%%%%%%%%
%%%% This section takes both parameter structure needed for fitting and
%%%% their names and:
%%%% 1) Merges together into a single *vector* the fitted variables and the
%%%% constants paired with their name (name, value, name2, value2, ....)
InputVariables = Modeling_PrepareVariables(InitialFitParam.ModelInput_FitVariables_names,...% Enter the names
                                           InitialFitParam.ModelInput_FitVariables,...      % Enter the values
                                           InitialFitParam.ModelInput_ConstVariables_names,...  % Enter the names
                                           InitialFitParam.ModelInput_ConstVariables);          % Enter the values


%%%% 2) Extracts the parameter structures needed for fitting and their names
%%%% and put them into entirely new parameter structure
if isempty(InputVariables);
    var = Modeling_GetVariables(Project,InitialFitParam.Model);
else
    var = Modeling_GetVariables(Project,InitialFitParam.Model,InputVariables{1:end});   % Extracts variables
end;
%%%% 3) Takes out only variable values, no constants no names. And puts
%%%% them into entirely new structure!
fit_var = var.fit_var;

%%%%%%%%%%%%%%%%%%%%% END prepare parameters loaded from file %%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%% START Prepare data for cost functions %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cost function C1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cost1_Stimulus = Modeling_CropStimulus(ExpData.Fig3A_Stimulus.t,ExpData.Fig3A_Stimulus.amp,InitialFitParam.C1.StimulusTimeSpan);
Cost1_Recording = Modeling_CropRecording(ExpData.Fig3A_Recording.ti,...
                    ExpData.Fig3A_Recording.ampi,...
                    InitialFitParam.dt,Cost1_Stimulus.stim_t,Cost1_Stimulus.RecShift);
Cost1_RecordingCostFun = Modeling_RecordingCostFun2(Cost1_Recording.rec_t,Cost1_Recording.rec_amp,...
                            InitialFitParam.C1.RecordingCostFun_Interval,...
                            InitialFitParam.C1.RecordingCostFun_Value,InitialFitParam.dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cost function C3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cost3_Stimulus = Modeling_CropStimulus(ExpData.Fig6D_10ms_Stimulus.t,ExpData.Fig6D_10ms_Stimulus.amp,InitialFitParam.C3.StimulusTimeSpan);
Cost3_Recording =  Modeling_CropRecording(ExpData.Fig6D_10ms_Recording.ti,ExpData.Fig6D_10ms_Recording.ampi,...
                    InitialFitParam.dt,Cost3_Stimulus.stim_t,Cost3_Stimulus.RecShift);
Cost3_RecordingCostFun = Modeling_RecordingCostFun2(Cost3_Recording.rec_t,Cost3_Recording.rec_amp,...
                            InitialFitParam.C3.RecordingCostFun_Interval,...
                            InitialFitParam.C3.RecordingCostFun_Value,InitialFitParam.dt);
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cost function C14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
Cost14_Stimulus = Modeling_CropStimulus(ExpData.Fig2A_Stimulus.t,ExpData.Fig2A_Stimulus.amp,InitialFitParam.C14.StimulusTimeSpan);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% START Run fitting %%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting done in parallel Matlab instances. "Parallel computing" toolbox required 
% Fitting process looks for global optimum
tic
% check if parallel pool is open                                  % Opens multiple matlab instances
% If no pool, do not create new one.
if isempty(gcp('nocreate'))
    parpool
end

opt = optimset('MaxFunEval',3000);
ff = @ChiSquare;                                % Assigns the handle of fitting function

%%%% various types of minimum search are written down here.
%%%% global search "simulannealbnd" requires options provided by "saoptimset"
%%%% provided by "saoptimset".

[xout,fval] = fminsearch(ff,fit_var,opt);
% [xout, fval] = fminunc(ff,fit_var,opt);
% [xout,fval] = patternsearch(ff,fit_var);
% [xout,fval] = ga(ff,length(fit_var));
% options = saoptimset('ObjectiveLimit',0.002);   % Sets objective limit (whatever that is)
% [xout,fval] = simulannealbnd(ff,fit_var,...     % Runs the fitting  - search for *global* maximum
%     zeros(length(fit_var),1),inf(length(fit_var),1),options);

tiktak = toc

% matlabpool close    % Closes multiple matlab instances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Run fitting %%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% START Save results & plot graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FittedVariables = xout;
FittedVariables_names = var.fit_var_name;

if isempty(var.const_var)
    ConstantVariables = [];
    ConstantVariables_names = [];
else
    ConstantVariables = var.const_var;
    ConstantVariables_names = var.const_var_name;
end;

%   Packs together variable & constant names
NamesOfVariables = vertcat({'Fitted variables:'},{[]},Modeling_ListOfVariables(FittedVariables_names,FittedVariables),{[]},...
    {'Constant variables:'},{[]},Modeling_ListOfVariables(ConstantVariables_names,ConstantVariables));
R2 = fval;


%   preparing the iteration number for the file names

%   Find the last existing number in the directory
%   Custom written function 'filename' is used to avoid problems using
%   scripts on different OS platforms
wh = what(filename(inputDir,{inputDir,Project},'Results')); %   looks whats in the dir


% fnd = find(strfind (wh.mat,[Project,'_Results_'])>0);

fnd = find(strncmp (wh.mat,[Project,'_Results_'],12)>0);

%   checks if there are (no) files in the dir and how many are there
if isempty(fnd)
    maxi = 0;
else
    maxi = NaN(1,length(fnd));
    for stri = 1 : length(fnd)
        StrNum = strrep(wh.mat{fnd(stri)},[Project,'_Results_'],'');
        StrNum = strrep(StrNum,'.mat','');
        maxi(stri) = str2double(StrNum);
    end;
end;
                          

MaxInd = max(maxi);   
%%% Following lines copy fitted variables, initial parameters, constants,
%%% its appropriate names and computing time consuption into a .mat file,
%%% with name Results_n+1.mat
save(filename(inputDir,{inputDir,Project,'Results'},[Project,'_Results_',num2str(MaxInd+1),'.mat']),...
    'InitialFitParam','FittedVariables','FittedVariables_names','ConstantVariables','ConstantVariables_names','R2',...
    'NamesOfVariables','tiktak');

%%% Following lines copy the used .m scripts into the same directory as
%%% the results. This ensures reproducibility of the results.

copyfile(filename(inputDir,{inputDir},'Modeling_DRG.m'),...
    filename(inputDir,{inputDir,Project,'Results'},[Project,'_Function_',num2str(MaxInd+1),'.m']),'f');
copyfile(filename(inputDir,{inputDir},'Modeling_DRG_TCM_Engine.m'),...
    filename(inputDir,{inputDir,Project,'Results'},[Project,'_Model_',num2str(MaxInd+1),'.m']),'f');

copyfile(filename(inputDir,{inputDir,Project},horzcat(Project,'_InitialFitParameters.mat')),...
    filename(inputDir,{inputDir,Project,'Results'},...
    [Project,'_InitialFitParameters_',num2str(MaxInd+1),'.mat']),'f');


clearvars -except Project Project_in MaxInd

%%% Following line outputs results, plots graphs and makes .pdf/.ps


% Modeling_DRG_Report(Project_in,horzcat('Results_',num2str(MaxInd+1)));
Modeling_DRG_Report(Project,horzcat(Project,'_Results_',num2str(MaxInd+1)));



%%%%%%%%%%%%%%%%%%%%%% END Save results & plot graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%% START Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y = ChiSquare(fit_var)
        
        model = InitialFitParam.Model;
        dt = InitialFitParam.dt;
        output = {};
        
        %%%%%%%%%%% Construct vector of variables %%%%%%%%%%%%%%%%%%
        if length(fieldnames(var)) > 2
            variables = horzcat(fit_var,var.const_var);
            variables_names = horzcat(var.fit_var_name,var.const_var_name);
        else
            variables = fit_var;
            variables_names = var.fit_var_name;
        end;
        
       %%%%%%%%%%% Penalize if variables are set outside limits %%%%%%%%%%%%
       FitVar_lim1 = InitialFitParam.ModelInput_FitVariables_limits(:,1);
       FitVar_lim2 = InitialFitParam.ModelInput_FitVariables_limits(:,2);
       FitVar_tol = InitialFitParam.ModelInput_FitVariables_tol(:);

%         for vr = 1:length(fit_var);
        parfor vr = 1:length(fit_var);
            vv(vr) = Modeling_ParameterLimits(fit_var(vr),...
                        FitVar_lim1(vr),FitVar_lim2(vr),FitVar_tol(vr));
        end;
        V = mean(vv(vv ~= 0));
        % prepares empty vector for paid penalties
        C = zeros(1,length(InitialFitParam.C_Weights)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%% Fit to cost functions  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% modeling the response to simulated stimulus
        %%%     peaksPredicted  ... simulated current peaks     

                    
%%% Fitting the Hao Figure 3A (control)
%%% Parameters
    % C1 compare current traces          
    % C7 restrict current traces from Fig3A to approach 0 mV in steady-state
    % C9 keep current traces before the stimulus in Fig3A above 0.02
    % C2 compare intensity response curves
    % C10 penalize if min peak goes below -1.1


        if InitialFitParam.C_Weights(1) ~= 0 || InitialFitParam.C_Weights(2) ~= 0 || InitialFitParam.C_Weights(7) ~= 0 ||...
            InitialFitParam.C_Weights(9) ~= 0 || InitialFitParam.C_Weights(10) ~= 0
        
            %   Preparing variables for parfor loop
            ProtocolIndex1 = InitialFitParam.C1.IndexOfStimuliToFit;
            tT = Cost1_Stimulus.stim_t(ProtocolIndex1,:);
            ampT = Cost1_Stimulus.stim_amp(ProtocolIndex1,:);
            cost1Rec_amp = Cost1_Recording.rec_amp(ProtocolIndex1,:);
            cost1Cost_amp = Cost1_RecordingCostFun.cost_amp(ProtocolIndex1,:);
            cost1Rec_t = Cost1_Recording.rec_t(ProtocolIndex1,:);
            peaksMeasured = -ExpData.Fig3A_Diagram.y'/100;  
            
            CW1 = InitialFitParam.C_Weights(1);
            CW2 = InitialFitParam.C_Weights(2);
            
            CW7 = InitialFitParam.C_Weights(7);
            CW9 = InitialFitParam.C_Weights(9);
            CW10 = InitialFitParam.C_Weights(10);
            
            wFig3 = InitialFitParam.C1.RecordingWeights;
            wFig3A = InitialFitParam.C2.PointWeights;
            
            [out,cost] = computeFig3 (model,tT,ampT,ProtocolIndex1,...
                               cost1Rec_amp,cost1Cost_amp, cost1Rec_t,...
                                variables,variables_names,dt,...
                                peaksMeasured,...
                                 CW1,CW7,CW9,CW2,CW10,...
                                 wFig3,wFig3A);
            C(1) = cost(1); C(2) = cost(2); C(7) = cost(3); C(9) = cost(4); C(10) = cost(5);
            output.Fig3A = out;
        end;
        
  
%%% Fitting the Hao Figure 6D & 6E
%%% Parameters
    % C3 compare current traces          
    % C12 compare recovery peaks of intensity response curves

        if InitialFitParam.C_Weights(3) ~= 0 || InitialFitParam.C_Weights(12) ~= 0
            
            %   Preparing variables for parfor loop (slicing?)
            ProtocolIndex3 = InitialFitParam.C3.IndexOfStimuliToFit;
            tT = Cost3_Stimulus.stim_t(ProtocolIndex3,:);
            ampT = Cost3_Stimulus.stim_amp(ProtocolIndex3,:);
            cost3Rec_amp = Cost3_Recording.rec_amp(ProtocolIndex3,:);
            cost3Cost_amp = Cost3_RecordingCostFun.cost_amp(ProtocolIndex3,:);
            peaksMeasured_recovery = ExpData.Fig6E_Diagram.y;
            
            
            CW3 = InitialFitParam.C_Weights(3);
            CW12 = InitialFitParam.C_Weights(12);
            wFig6E = InitialFitParam.C3.PointWeights;
            
            [out,cost] = computeFig6 (model,tT,ampT,ProtocolIndex3,...
                               cost3Rec_amp,cost3Cost_amp, ...
                                variables,variables_names,dt,...
                                 peaksMeasured_recovery,...
                                 CW3,CW12,...
                                 wFig6E);
            
            C(3) = cost(1); C(12) = cost(2);
            output.Fig6 = out;
        
        end;

% %%%%%%%%C8 & C13 Penalization for activation and inactivation
        
        %%%%%%%%%%% C8 penalize if h0 goes above 0.1 at x = 0 %%%%%%%%%%%%
        if InitialFitParam.C_Weights(8) ~= 0
            Hh = variables(strncmp('Hh', variables_names, length('Hh')));
            Xh = variables(strncmp('Xh', variables_names, length('Xh')));
            z = h0(Hh,Xh);
            if z > 0.1; C(8) = 1000; else C(8) = 0; end;
        end;
        
        %%%%%%%%%%% C13 penalize if m0 goes above 0.1 at x = 0 %%%%%%%%%%%%
        if InitialFitParam.C_Weights(13) ~= 0
            Hm = variables(strncmp('Hm', variables_names, length('Hm')));
            Xm = variables(strncmp('Xm', variables_names, length('Xm')));
            z = m0(Hm,Xm);
            if z > 0.1; C(13) = 1000; else C(13) = 0; end;
        end;
        

        
        
%%%%%%%%%%%%%%%%%%%%%% C14 compare adapted onset curves Fig2B %%%%%%%%%%%%%%%%%%%%%%%%
    %%% peaksPredicted_recovery     ... searches g(min) at certain
        %%%                                 part of response
        if InitialFitParam.C_Weights(14) ~= 0
            
            %   Preparing variables for parfor loop (slicing et al.)
            ProtocolIndex14 = InitialFitParam.C14.IndexOfStimuliToFit;
            tT = Cost14_Stimulus.stim_t(ProtocolIndex14,:);
            ampT = Cost14_Stimulus.stim_amp(ProtocolIndex14,:);
            peaksMeasured_rePoke = ExpData.Fig2B_Diagram.y'./max(ExpData.Fig2B_Diagram.y);
            wFig2 = InitialFitParam.C14.PointWeights;
            
            
            [out,cost] = computeFig2 (model,tT,ampT,ProtocolIndex14,...
                             variables,variables_names,dt,...
                              peaksMeasured_rePoke,...
                              wFig2);
            
            C(14) = cost(1);
            output.Fig2A = out;
        end;
        
%%%%%%%%%%%%%%%%%%%%%% C4, C5, C6 figure 3GHI Hao2010%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        %   C4 .. measures difference between modeled Tau Inact & Tau Adapt
        %           and measured Tau Inact & Tau Adapt
        %   C5 .. measures difference between modeled adaptive shift and
        %           measured adaptive shift
        %   C6 .. measures difference between modeled inactivated fraction
        %           and measured inactived fraction
        
        if InitialFitParam.C_Weights(4) ~= 0 || InitialFitParam.C_Weights(5) ~= 0 || InitialFitParam.C_Weights(6) ~= 0
        
            stimulus = ExpData.Fig3C_Stimulus;
            control_pks = output.Fig3A.Imax;
            control_amps = ExpData.Fig3A_Stimulus.amp(:,2)';
            
            Fig3I = ExpData.Fig3I.InactAdapt;
            Fig3I_PW = ExpData.Fig3I.PW_inact_adapt;
            tauFig3GH = ExpData.Fig3HG_InactAdapt;
            
            CW4 = InitialFitParam.C_Weights(4);
            CW5 = InitialFitParam.C_Weights(5);
            CW6 = InitialFitParam.C_Weights(6);
            
            
            [out, cost] = compute3GHI (model,stimulus,...
                                variables,variables_names,dt,...
                                control_pks,control_amps,...
                                tauFig3GH, Fig3I,...
                                CW4,CW5,CW6,Fig3I_PW);
        
            C(4) = cost(1); C(5) = cost(2); C(6) = cost(3);
            output.Fig3GHI = out;
        end;
 
%%%%%%%%%%Final remarks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        weightC = InitialFitParam.C_Weights/sum(InitialFitParam.C_Weights);
        
        C;
        CC = sum(C.*weightC);
        if isnan(CC); CC = Inf; end;
        y = nanmean([CC,V]);
 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END Fitting results %%%%%%%%%%%%%%%%%%%%%%%%%

%% %%% Additional functions for evaluation COST
    function z = h0(Hh,Xh)
        
        z = 1./(1 + exp(-Hh*(0 - Xh)));
        
    end

    function z = m0(Hm,Xm)
        
        z = 1./(1 + exp(-Hm*(0 - Xm)));
        
    end

end
