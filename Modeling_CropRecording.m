%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function output = Modeling_CropRecording(rec_ti,rec_ampi,dt,stim_t,RecShift)

for dr = 1:length(rec_ti);
    if isnan(rec_ti{dr}) == 0 % ignore, if only NaNs
        RecordingTime(dr,:) = 0:dt:sum(stim_t(dr,:));
        Recording_t = rec_ti{dr} - RecShift;
        Recording_amp = rec_ampi{dr};
        RecordingAmp(dr,:) = interp1(Recording_t,Recording_amp,RecordingTime(dr,:));
        
        %%%%%%%%%% Replace NaNs with the closest non-NaN number %%%%%%%%%
        
        if isnan(sum(RecordingAmp(dr,:)))
            Ind_NaN = find(isnan(RecordingAmp(dr,:)));
            if isempty(find(diff(Ind_NaN) > 1, 1));
                if Ind_NaN(1) == 1
                    RecordingAmp(dr,Ind_NaN) = RecordingAmp(dr,Ind_NaN(end) + 1);
                else
                    RecordingAmp(dr,Ind_NaN) = RecordingAmp(dr,Ind_NaN(1) - 1);
                end;
            else
                ind_break = find(diff(Ind_NaN) ~= 1);
                RecordingAmp(dr,Ind_NaN(1:ind_break)) = RecordingAmp(dr,Ind_NaN(ind_break) + 1);
                RecordingAmp(dr,Ind_NaN(ind_break:end)) = RecordingAmp(dr,Ind_NaN(ind_break+1) - 1);
            end;
        end;
    end;
end;

output.rec_t = RecordingTime;
output.rec_amp = RecordingAmp;