%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Utiliy used in generating stimulus time course from stimulus segment
%   changes and durations.

function [stim_t,stim_amp] = Modeling_GenerateStimulus(t,StimInterval,StimAmplitude)

stim_amp = NaN(size(StimInterval,1),length(t));
stim_t = repmat(t,size(StimInterval,1),1);
for val2 = 1:size (StimInterval,1)
        x0 = 0;
        y0 = 0;
        x = x0;
        y = y0;
        
        for val1 = 1:size (StimInterval,2)
            if StimInterval(val2,val1) ~= 0
                x1 = x0 + StimInterval(val2,val1);
                y1 = y0 + StimAmplitude(val2,val1);
                x0 = x1;
                y0 = y1;
                x = [x,x1];
                y = [y,y1];
            end;
        end;
%   interp does not work with NaN values, therefore the NaN should be removed beforehand       
        stim_amp(val2,:) = interp1(x,y,t,'linear');

end;
    


