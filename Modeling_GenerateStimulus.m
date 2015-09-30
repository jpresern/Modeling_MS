%%% Written by Aleš Škorjanc at some point in 2011

%%% Adapted and accelerated by Janez Prešern, 2014

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
    


% for a = 1:size(StimInterval,1);
%     
%     if sum(StimInterval(a,:)) ~= 0
%         
%         Stim_Interv_temp(a,1) = 0;
%         Stim_Interv_temp(a,[2:length(StimInterval(a,:))+1]) = cumsum(StimInterval(a,:));
%         
%         for b = 1:length(StimInterval(a,:));
%             
%             ind = find(t >= Stim_Interv_temp(a,b) & t <= Stim_Interv_temp(a,b+1));
%             
%             if isempty(ind) == 0
%                 
%                 if StimAmplitude(a,b) == 0;
%                     stim_amp(a,ind) = sum(StimAmplitude(a,1:b-1));
%                 else
%                     stim_amp(a,ind) = sum(StimAmplitude(a,1:b-1)) +...
%                         (t(ind) - sum(StimInterval(a,1:b-1)))*(StimAmplitude(a,b)/StimInterval(a,b));
%                 end;
%             end;
%         end;
%         
%         ind2 = find(t > sum(StimInterval(a,:)));
%         stim_amp(a,ind2) = sum(StimAmplitude(a,:));
%         
%     end;
%     
% end;


