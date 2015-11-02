%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function draws figures and produces .ps (& .pdf) files
%   During normal modeling it is called by function Modeling_DRG.m. It is
%   possible to execute this function "solo" for which it is required to
%   uncomment/add the line beneath the "function [output]..." line,
%   containing project name ('Project=' which is working folder, 
%   containing Results subfolder, which in turn contains results you are 
%   interested in). The modeling results are specified in 'fn' where the
%   correct number should be added to '_Results_'.
%   See examples below.

%   Input:
%       Project ... project folder, containing results subfolder.
%       fn      ... file name, containing the name of result file

%   Output
%       a .mat file, called Output_DRG_Results_X.mat, X matching the number
%       in input parameters. Output_DRG_Results_X contains most of the
%       computations of this single run, also useful for a crosscheck.

function [output] = Modeling_DRG_Report(Project,fn)
%% % Loads results from Results_X.mat file

inputDir = cd;
load(filename(inputDir,{inputDir,Project,'Results'},[fn,'.mat']));
ExpData = load(filename(inputDir,{inputDir,Project},[InitialFitParam.ExperimentalData,'.mat']));
SavePath = filename(inputDir,{inputDir,Project,'Results'},[]);

var = horzcat(FittedVariables,ConstantVariables);
var_names = horzcat(FittedVariables_names,ConstantVariables_names);


%%% Set new color map
cmap = [0 0 1;          % blue
        0.15 0.5 0.7;  % light blue/cyan 
        1.0 0.6 0.2;    % orange
        0.5 0.0 0.5;    % violet
        1 0 0;          % red
        0 0.7 0.3;      % green
        0.7 0.7 0.7;    % gray
        0 0 1;          % blue
        0 0.7 0.7;      % light blue/cyan
        1.0 0.6 0.2;    % orange
        0.5 0.0 0.5;    % violet
        1 0 0;          % red
        0 0.7 0;        % green
        0.5 0.5 0.5];   % darker gray
    
cmapBlue = [0 0 1;          % blue
        0.15 0.5 0.8;  % light blue/cyan 
        0.15 0.5 0.7;
        0.1 0.5 0.6;
        0.1 0.6 0.7;    
        0.05 0.6 0.6;   
        0.05 0.7 0.6;
        0.1 0.7 0.6;    
        0.15 0.75 0.6;      
        0.1 0.7 0.6;
        0.05 0.7 0.6;
        0.15 0.5 0.7;    
        0.15 0.5 0.8;    
        1 0 0;          
        0 0.7 0;        
        0.5 0.5 0.5];   
    
%%% Set a marker map
marmap = ['o','+','v','^','p','x'];

%%
f = [];
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Current response to simple stimulus and I-R curve%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fig, output.BasicIR] = Draw_BasicIR(InitialFitParam.dt, ExpData.BasicIR_Stimulus.amp, ExpData.BasicIR_Stimulus.t,...
                        ExpData.BasicIR_Recording.ampi, ExpData.BasicIR_Recording.ti, var,...
                        [InitialFitParam.ModelInput_FitVariables,ConstantVariables],...
                        InitialFitParam.ModelInput_FitVariables_limits, var_names,...
                        InitialFitParam.C_Weights,InitialFitParam.C1.RecordingWeights,...
                        cmap, cmapBlue, fn);
    
f = fig;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% m0 and h0 curves   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = Draw_MH(ExpData.BasicIR_Stimulus.amp,...
        output.BasicIR.results{1,size(output.BasicIR.results,2)}.Predicted_m0.x,...
        output.BasicIR.results{1,size(output.BasicIR.results,2)}.Predicted_m0.y,...
        output.BasicIR.results{1,size(output.BasicIR.results,2)}.Predicted_h0.x,...
        output.BasicIR.results{1,size(output.BasicIR.results,2)}.Predicted_h0.y);
    

f = [f,fig];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Drawing of inactivated fraction via double-pulse %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fig, output.InactFract] = Draw_InactFract(InitialFitParam.dt,...
                        ExpData.InactFract_Stimulus.amp,ExpData.InactFract_Stimulus.t,...
                        var, var_names,...
                        cmap, cmapBlue, fn);
    
f = [f,fig];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Recovery from inactivations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fig, output.Recovery] = Draw_Recovery (InitialFitParam.dt,...
                        ExpData.Recovery_Stimulus.amp,...
                        ExpData.Recovery_Stimulus.t,...
                        var, var_names,...
                        cmap, [cmapBlue(1:13,:);cmapBlue(1:13,:)], fn);

f = [f, fig]; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Two-step conditioning, drawing selected cond. amplit. %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
[fig] = Draw_Conditioning (InitialFitParam.dt,...
                    ExpData.Conditioning_Stimulus,...
                    output.BasicIR.stimulus.ampMax,...
                    output.BasicIR.model.Imax,...
                    output.BasicIR.model.x50k50,...
                    var, var_names,...
                    cmap, cmapBlue, marmap, fn);

f = [f, fig];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Two-step conditioning, computation only %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output.Conditioning = computeConditioning (horzcat(InitialFitParam.Model,'_Report'),...
                    ExpData.Conditioning_Stimulus,...,
                    var,var_names,InitialFitParam.dt);               



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Two-step conditioning, adaptive shift and inact. frac %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fig, output.InactAdapt_1] = Draw_InactAdapt_1(output.Conditioning.model.stimAmp,...
                                output.Conditioning.model.peakRecovery,...
                                output.BasicIR.model.x50k50,...
                                cmap, cmapBlue, marmap);

f = [f, fig]; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Two-step conditioning, draw tau inac, adapt, separation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fig, output.InactAdapt_2] = Draw_InactAdapt_2 (output.InactAdapt_1.model.tauAct,...
                                    output.InactAdapt_1.model.tauInact,...
                                    output.InactAdapt_1.model.adaptShift,...
                                    output.Conditioning.model.peakRecovery,...
                                    cmap);
                                
f = [f, fig];

%%
%%%%%%%%%  Create PDF   - print me off, before you go-go...%%%%%%%%%%%%%%%%%%%%%%%
CreatePdf(f,fn,SavePath);


%% Output to variable "output."
fnn = horzcat('Output_',fn);
% SavePath = horzcat(inputDir,'\',Project,'\Results\');
save(horzcat(SavePath,fnn),'output');

%% Closes created figures

close (f);
end
