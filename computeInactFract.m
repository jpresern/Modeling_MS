%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function computes the desensitization with a prepulse and compares the
%   I-R curve with the one published by Hao2010, fig 2.
%   Function requires:
%       model   ..        model type
%       tT      ..        stimulus time
%       ampT    ..        stimulus amplitude at times tT
%       ProtocolIndex ..  which of the traces/points shall we fit
%       variables  ..     variable values as inserted by fminsearch
%       variables_names ..names of variables
%       peaksMeasured_rePoke .. experimentally measured recovery peaks
%       wFig2   ..        point weights

%   Function outputs:
%       outputs.peakRePoke..  maximum current g at various stimuli amplitudes
%       cost              ...       computed costs

function [out,varargout] = computeInactFract (model,tT,ampT,ProtocolIndex,...
                                variables,variables_names,dt,...
                                peaksMeasured_rePoke,...
                                wFig2)


%   Preparing variables for parfor loop (slicing et al.)
peaksPredicted = nan(size(ProtocolIndex));
peaksPredicted_rePoke = nan(size(ProtocolIndex));
            
%%% Fitting the figure 2A - trapezoid stimuli in poke - repoke fashion

% for w14 = 1:length(ProtocolIndex)    
parfor w14 = 1:length(ProtocolIndex)

    %%% computes the g elicited by stimulus
    g = Modeling_DRG_TCM_Engine(model,tT(w14,:),ampT(w14,:),...
        variables,variables_names,dt);
    %%% extracts the peaks
    peaksPredicted(w14) = min(g(1:825));  
    peaksPredicted_rePoke(w14) = min(g(826:end)); % selecting the part of the modeled 
                                                  % response we are interested in                
end;

%%%%%%%%%%% C14 compare recovery curves %%%%%%%%%%%%%%%%%
%%%         compares the measured peaks against the model

peaksRepoke = peaksPredicted_rePoke./min(peaksPredicted_rePoke);
c1 = sum((wFig2(ProtocolIndex)/length(wFig2(ProtocolIndex))).*...
    (peaksRepoke - peaksMeasured_rePoke(ProtocolIndex)).^2)./...                       
    abs((mean(peaksMeasured_rePoke(ProtocolIndex)))*length(peaksPredicted)); 

%%%%%  output section
cost = c1;
out.peaksRepoke = peaksPredicted_rePoke;
varargout = {cost};
end