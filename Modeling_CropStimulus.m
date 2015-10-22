%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Utility cropping stimulus at the desired length (at the end of
%   relevance) as set in InitFitParam

function output = Modeling_CropStimulus(stim_t,stim_amp,TimeSpan)

% Function crops the ramp stimulus. Stimulus is defined by the vector of 
% time intervals of each ramp phase and their corresponding amplitudes.
% Amplitudes of static phases have value 0. Amplitudes of dynamic phases 
% have the value of the relative change in the stimulus strength.

for dr = 1:size(stim_t,1)
    if isnan(sum(stim_t(dr,:))) == 0
        if sum(TimeSpan) ~= 0  
            IntervTemp = stim_t(dr,:);
            AmpTemp = stim_amp(dr,:);
            VelTemp = stim_amp(dr,:)./stim_t(dr,:);
            IntervCumsum = cumsum(IntervTemp);
            
            %%%%%%% Determine the intervals within the TimeSpan limits %%%%%%%%
            
            drind1 = find(TimeSpan(1) <= IntervCumsum, 1 );
            if isempty(drind1); drind1 = 1; end; % if the left limit is higher then the length of the stimulus take the first interval
            drind2 = find(TimeSpan(2) <= IntervCumsum, 1 );
            if isempty(drind2); drind2 = length(IntervTemp); end; % if the right limit is higher then the length of the stimulus take the last interval
            interval(dr,[drind1:drind2]) = IntervTemp(drind1:drind2);
            amplitude(dr,[drind1:drind2]) = AmpTemp(drind1:drind2);
                        
            %%%%%%% Determine the new borders of the stimulus %%%%%%%%%%%%%%%
            
            if drind1 ~=1; shift1 = TimeSpan(1)-IntervCumsum(drind1-1);
            else 
                shift1 = TimeSpan(1);
            end;
            
            shift2 = IntervCumsum(drind2)-TimeSpan(2);
            
            interval(dr,1) = interval(dr,1) - shift1;
            interval(dr,drind2) = interval(dr,drind2) - shift2;
            
            amplitude(dr,drind2) = VelTemp(drind2)*interval(dr,drind2);
            
        else
            interval(dr,:) = stim_t(dr,:);
            amplitude(dr,:) = stim_amp(dr,:);
        end;
    end;
end;

output.stim_t = interval; 
output.stim_amp = amplitude;
output.RecShift = TimeSpan(1);