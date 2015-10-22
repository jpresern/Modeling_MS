%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Utiliy, currently not in use

function [tT,ampT] = Modeling_DRG_DeadSpace(t,amp,x0)

%%%
% This function finds the ramps in the stimuli. The stimulus itself is
% written as the duration/max.amplitude pairs. This means that ascending
% ramps have positive max.amplitude/tmax-tmin ratio, while descending
% ramps have negative ratio. In tflat parts of the stimuli this ratio
% equals 0.

%We are interested only in those ramps which start at the bottom, e.g.
%amplitude = 0. As stimuli can be different (simple ramps, conditioning
%ramps - ramp on the top of the ramp, ramp after ramp for recovery) it is
%necessary to distinguish the ramps which start at the initial amplitude
%(0) from those in conditioning stimuli which start at the amplitudes >0.
%To provide this information we will calculate the cumulative sum of
%mentioned ratio for all the parts of the stimuli. We are interested in
%parts where cumulative sum equals zero and search for the their index
%values. We are interested in the last 0 value just before the ramp and
%first zero value after each ramp. In this case we will take index value+1.


amp_sum = cumsum(amp);
c1 = 0;
c2 = 0;
tT =  t;
ampT = amp;

for a = 1:length(amp_sum)-1
    if amp_sum(a) == 0 && amp_sum(a+1) > 0
        c1 = c1 + 1;
        indL(c1) = a + 1;
        v = amp(a+1)/t(a+1);
        tT(a) = tT(a) + x0/v;
        tT(a+1) = tT(a+1) - x0/v;
        ampT(a+1) = ampT(a+1) - x0; 
    end;
    if amp_sum(a) ~= 0 && amp_sum(a+1) == 0
        c2 = c2 + 1;
        indR(c2) = a + 1;
        v = amp(a+1)/t(a+1);
        tT(a+1) = tT(a+1) + x0/v;
        ampT(a+1) = ampT(a+1) + x0; 
        if a + 2 <= length(amp_sum)
            tT(a+2) = tT(a+2) - x0/v;
        end;
    end;
end;


        
    