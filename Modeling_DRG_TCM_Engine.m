%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%................DRG_TransducerCurrentModel.............................
%
% This function simulates transduction currents of DRG neurons.
% Models are based on the paper Hao and Delmas 2010.
%
%
%
%
%
% The model is a phenomenological description of the DRG receptor
% current recorded in a voltage-clamp mode at holding potential -60mV.
% Current measured at a constant holding potential is linearly proportional
% to the conductance of the transduction channels at that potential, therefore
% we use conductance instead of currents in the model (g=I/E).
%




% Input parameters:
% 
% ModelName... name of the model (use ModelName_Report for using the
% extended version, which produces data for report)
% time... ramp stimulus time intervals (length of all dynamic and static
% phases in ordered as they come)
% amplitude... amplitudes of individual stimulus phases(amplitude of 
% static phases is 0; amplitude of dynamic phases determines the relative
% shift of stimulus strength during the dynamic phase)
% variables... values of variables used by the model
% variables_names... names of variables used by the model
% dt... integration step




function [g,varargout]=Modeling_DRG_TCM_Engine(ModelName,time,amplitude,variables,variables_names,dt)



% List of variables:

% time constant of activation
tau_m = variables(strcmp(variables_names,'tau_m'));
% time constant of inactivation
tau_h = variables(strcmp(variables_names,'tau_h'));
% scaling
A = variables(strcmp(variables_names,'A'));
% factor describing the linear relationship between adaptation and inactivation
B = variables(strcmp(variables_names,'B'));
% Boltzmann for activation: Am..scaling, Hm..slope, Xm..midpoint
Am = variables(strcmp(variables_names,'Am'));   %   constant
Hm = variables(strcmp(variables_names,'Hm'));   %   fitted variables
Xm = variables(strcmp(variables_names,'Xm'));   %   fitted variables
% Boltzmann for inactivation: Ah..scaling, Hh..slope, Xh..midpoint
Ah = variables(strcmp(variables_names,'Ah'));   %   constant
Hh = variables(strcmp(variables_names,'Hh'));   %   fitted variable
Xh = variables(strcmp(variables_names,'Xh'));   %   fitted variable
%   Unused at the moment
N = variables(strcmp(variables_names,'N'));     %   constant
M = variables(strcmp(variables_names,'M'));     %   constant
% P = variables(strmatch('P',variables_names,'exact'));
% R = variables(strmatch('R',variables_names,'exact'));



varargout = {[]};           %%%% Prepare empty structure

%%%%%%%%%%%%% Create a list of variables for output %%%%%%%%%%%%%%%

% if findstr('_Report',ModelName)     
if strfind(ModelName,'_Report')    %%% If "_Report" is found in model name, do write values in out.parameters
    out.parameters = [];
    if isempty(tau_m) == 0; out.parameters = vertcat(out.parameters,{horzcat('tau m = ',num2str(tau_m))}); end;
    if isempty(tau_h) == 0; out.parameters = vertcat(out.parameters,{horzcat('tau h = ',num2str(tau_h))}); end;
    if isempty(A) == 0; out.parameters = vertcat(out.parameters,{horzcat('A = ',num2str(A))}); end;
    if isempty(B) == 0; out.parameters = vertcat(out.parameters,{horzcat('B = ',num2str(B))}); end;
    if isempty(Hm) == 0; out.parameters = vertcat(out.parameters,{horzcat('Hm = ',num2str(Hm))}); end;
    if isempty(Xm) == 0; out.parameters = vertcat(out.parameters,{horzcat('Xm = ',num2str(Xm))}); end;
    if isempty(N) == 0; out.parameters = vertcat(out.parameters,{horzcat('N = ',num2str(N))}); end;
    if isempty(M) == 0; out.parameters = vertcat(out.parameters,{horzcat('M = ',num2str(M))}); end;
    if isempty(Am) == 0; out.parameters = vertcat(out.parameters,{horzcat('Am = ',num2str(Am))}); end;
    if isempty(Ah) == 0; out.parameters = vertcat(out.parameters,{horzcat('Ah = ',num2str(Ah))}); end;
end;

% All variables and their initial values are stored in workspace file in
% folder DRG_TCM
% coeff...coefficients of Boltzman functions fitted to data in Hao2010 for m_inf, a_inf and h_inf
% curves...data for m_inf, a_inf and h_inf obtained from Hao2010 with digitize2

t = 0:dt:sum(time);     % integration time interval computed from sum of the phase durations
g = NaN(1,length(t));   % preparing the empty variable for computed current
m = NaN(1,length(t));   % preparing the empty variable for computed activation
h = NaN(1,length(t));   % preparing the empty variable for computed inactivation

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL mh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Model equations:
%
% 1) g = g_tot * m^N * (1 - h)^M;
% 2) tau_m * (dm/dt) = - m + m0(x - B * h);
% 3) tau_h * (dh/dt) = - h + h0(x);
%
% 1) conductancy of mechanosensitive channels
% g...relative conductance of transduction channels
% g_tot...total conductance of transduction channels (=1)
% m...activation parameter
% N...power reflecting the number of channel subunits
% h...inactivation parameter (1-h ... fraction of un-inactivated channels)
% M...power reflecting the number of channel subunits

% 2) dynamics of activation
% tau_m...time constant of activation
% dm/dt...time derivative of activation parameter
% m0...onset activation curve (I/I_max(x)) (Hao and Delmas 2010, Fig3E)
% x...stimulus [microns](Hao and Delmas 2010)
% B...linear coefficient inactivation effect on adaptive shift 

% 3) dynamics of inactivation
% tau_h...time constant of inactivation
% dh/dt...time derivative of inactivation parameter
% h0...un-inactivated fraction (Hao and Delmas 2010, Fig3E)
% x...stimulus [microns](Hao and Delmas 2010)



% Calculate the stimulus in whole
x0 = 0;
y0 = 0;
x = x0;
y = y0;
for val1 = 1:length (time)
    if time(val1) ~= 0
        x1 = x0 + time(val1);
        y1 = y0 + amplitude(val1);
        x0 = x1;
        y0 = y1;
        %   x & y cannot be prealocated because the size of the stimulus segment
        %   is not the same in every loop
        x = [x,x1];
        y = [y,y1];
     end;
end;
stim = interp1(x,y,t,'linear');


% Calculate intial conditions: 
h_temp = mh_h0(stim(1));
m_temp = mh_m0(stim(1) - B*(h_temp));

%   calculate the response
if strcmp(ModelName,'DRG_TCM_Model_mh') % used for fitting
    for c = 1:length(t);
        m_temp = dt/tau_m * (-m_temp + mh_m0(stim(c) - B*(h_temp))) + m_temp;
        h_temp = dt/tau_h*(-h_temp + mh_h0(stim(c))) + h_temp;
        g(c) = -A*(m_temp^N)*(1 - h_temp)^M;
    end;
end;


if strcmp(ModelName,'DRG_TCM_Model_mh_Report') % used for report
    % Calculate intial conditions:
    % further abuse of computing power
    Measured_m0 = @(x)1./(1+exp(-1.315*(x-4.934)));
    [results.Measured_m0.x,results.Measured_m0.y] = fplot(Measured_m0,[0,10]);
    Measured_h0 = @(x)1./(1+exp(-0.94*(x-4.837)));
    [results.Measured_h0.x,results.Measured_h0.y] = fplot(Measured_h0,[0,10]);
    Predicted_m0 = @mh_m0; 
    [x,y] = fplot(Predicted_m0,[0,10]);
    results.Predicted_m0.x = x; results.Predicted_m0.y = y.^N;
    Predicted_h0 = @mh_h0;
    [x,y] = fplot(Predicted_h0,[0,10]);
    results.Predicted_h0.x = x; results.Predicted_h0.y = y.^M;
    
    results.t = t;
        
    for c = 1:length(t);
        m_temp = dt/tau_m * (-m_temp + mh_m0(stim(c) - B*(h_temp))) + m_temp;
        h_temp = dt/tau_h*(-h_temp + mh_h0(stim(c))) + h_temp;
        m(c) = m_temp^N;
        h(c) = (1 - h_temp)^M;
        g(c) = -A*(m_temp^N)*(1 - h_temp)^M;
    end;
    results.m = m;
    results.h = h;
    results.g = g;
    results.stim = stim;

    out.equations = {'1) g = g_tot * m^N * (1 - h)^M';...
        '2) tau_m * (dm/dt) = - m + m0(x - B * (1 - h))';...
        '3) tau_h * (dh/dt) = - h + h0(x)';...
        '4) m0 = Am./(1 + exp(-Hm*(w - Xm)))';...
        '5) h0 = y = Ah./(1 + exp(-Hh*(w - Xh)))'};
    
    varargout(1) = {out};
    varargout(2) = {results};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%................Subfunctions of Model_mh, Model_mh2..............
%
%
% Functions describe the onset transduction current m0 and inactivated fraction
% of transduction channels h0. All functions are Boltzman functions
% with formula: A/(1 + exp(-H*(x - X))). x is stimulus strength in
% microns.

    function y = mh_m0(w)
    
    y = Am./(1 + exp(-Hm*(w - Xm)));
    
    end
    
    
    function y = mh_h0(w)
    
    y = Ah./(1 + exp(-Hh*(w - Xh)));
    
    end
   
%......................................................................
end