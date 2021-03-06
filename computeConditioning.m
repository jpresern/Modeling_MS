%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function computes the responses to tests with different conditioning
%   amplitudes and durations. 
%   Function requires:
%       model   ..         model type
%       stimulus      ..   stimulus - structure containing 3D matrix with 
%                           stim times & stimulus amplitudes
%       varr        ..     variable values as inserted by fminsearch
%       varr_names ..      names of variables

%   Function outputs:
%   out.model.peakInitial ..    a matrix of peaks elicited by conditionig
%                               stimuli
%   out.model.peakRecovery..    a matrix of peaks elicited by test stimuli
%   out.model.stimAmp ..        a matrix of all stimuli amplitudes for all
%                               stimuli used


function out = computeConditioning (model,stimulus,...
                                varr,varr_names,dt)
                            
%   prepare empty array for the stimulus
stimAmp = nan(size(stimulus(1,1).t,1),size(stimulus,2));

%%%%%%%%%%%% calculates the stimulus amplitudes only, not for stimulation
for val0 = 1:size (stimulus,2)
    for val1 = 1:size (stimulus(1,val0).t,1)
        stimAmp(val1, val0) = sum (stimulus(1,val0).amp(val1,1:4,1));
    end;
end;
                            

%%%%%%%%%%%%% calculating the model responses

%   prepares the variables for the both the conditioning stimulus peak and
%   test stimulus    
peakInitial = NaN(size(stimulus(1,1).t,1),size(stimulus(1,1).t,3),size(stimulus,2));
peakRecovery = peakInitial;

%%%%%% Conditioning stimulus loop
% for l = 1 : size(stimulus,2)
parfor l = 1 : size(stimulus,2)
%   It turned out to be very efficient to have the parfor loop as the
%   outernmost for loop. The calculation time was cut down from 300-500 s
%   to ~ 150 s including opening the matlab pool (10 s) and closing it (5
%   s)
    %   parfor: prepare dummy variables for peakInitial and peakRecovery
    dummy1=nan(size(stimulus(1,l).t,1),size(stimulus(1,l).t,3));
    dummy2=nan(size(stimulus(1,l).t,1),size(stimulus(1,l).t,3));
    %   parfor: prepare dummy variables for index ix1 & ix2
    dummy3=nan(size(stimulus(1,l).t,1),size(stimulus(1,l).t,3));
    condResp = [];
    
    %%%%%% Delay loop
    for n = 1 : size (stimulus(1,l).t,3) 
%         n
        %%%%%% Amplitude loop, starts at value l to skip the rows with NaN    
        for m = l : size(stimulus(1,l).t,1)  
            %   clears the variables, makes sure no artefacts are drawn (or
            %   extracted)            
            ix2 = [];
            pks = [];
            %    checks if stimulus is a number or NaN
            if ~isnan(stimulus(1,l).t(m,1:3,n)) 
               g = Modeling_DRG_TCM_Engine(model,stimulus(1,l).t(m,:,n),...
                            stimulus(1,l).amp(m,:,n),varr,varr_names,dt);
                %   calculates index to grab peaks in the correct time window
                %   stores the current computed in the response to the
                %   conditioning stimuli, storing only max. delay current to ensure the
                %   proper length of the vector for further computation;
                %   However it stores peak values and index
               
                if m == l
                    condResp = g;
                    [dummy1(m,n), dummy3(m,n)] = min(condResp);

                    %   problem with no peaks at the 0 stimuli
                    %   if dummy3(m,n) == length(condResp)
                    if dummy3(m,n) >= int64(sum(stimulus(1,l).t(m,1:3,n))/dt);
                        dummy3(m,n) = int64(sum(stimulus(1,l).t(m,1:3,n))/dt);
                    end;
                else
                    %   if no value is stored, copy the value along the
                    %   dimensions of the matrix for simpler use
                    dummy1(m,:) = dummy1(l,:);
                    dummy3(m,:) = dummy3(l,:);
                end;
   
                %   recovery peak extraction
                [pks, ix2] = findpeaks (-g(dummy3(m,n):int64(sum(stimulus(1,l).t(m,1:end-2,n))/dt)),...
                                'npeaks',1);
                ix2 = dummy3(m,n) + ix2;
            end;
            dummy2(m,n) = 0;
            %   subtracting the value of the conditioning response (offset)
            if ~isempty(pks)
                dummy2(m,n) = -pks-condResp(ix2);
            end;
        end;
    end;
    peakInitial(:,:,l) = dummy1;
    peakRecovery(:,:,l) = dummy2;
end;

out.model.peakInitial = peakInitial;
out.model.peakRecovery = peakRecovery;
out.model.stimAmp = stimAmp;

end