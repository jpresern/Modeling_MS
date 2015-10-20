%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = Modeling_RecordingCostFun2(rec_t,rec_amp,CostFunInterv,CostFun_Value,dt)


CostFun_Time = rec_t;
CostFun_Amp = NaN(size(rec_t));
CostFun_Value = CostFun_Value./sum(CostFun_Value);

% checks if all the parameters are real or any contains string (peak)
if iscell(CostFunInterv)
    for st = 1: length(CostFunInterv);
        CostFunIntervChar(st) = {num2str(CostFunInterv{st})};
    end;
end;

for csf = 1 : size(CostFun_Time,1)
    
    if sum(rec_amp(csf,:)) ~= 0
        
        if iscell(CostFunInterv)
            
            if strmatch('peak',CostFunIntervChar);
                
                CostFunInterv_temp = CostFunIntervChar;
                
                for p = 1:length(strmatch('peak',CostFunIntervChar))

                    indStrPeak = find(strcmp('peak',CostFunInterv_temp));
                    
                    indPeak1 = dsearchn(rec_t(csf,:)',str2num(CostFunInterv_temp{indStrPeak(1) - 1}));
                    indPeak2 = dsearchn(rec_t(csf,:)',str2num(CostFunInterv_temp{indStrPeak(1) + 1}));
                    
                    %computing of peak index (zakajtakokomplicirano?)
                    [~,I]=max(abs(rec_amp(csf,indPeak1:indPeak2)));
                    peak_ind = I + indPeak1;
                    %computing of peak time
                    PeakTime = rec_t(csf,peak_ind);
                    
                    CostFunInterv_temp(1,indStrPeak(1)+2:end+1) = CostFunInterv_temp(1,indStrPeak(1) + 1:end);
                    CostFunInterv_temp(1,indStrPeak(1):indStrPeak(1)+1) = {num2str(PeakTime)};
                    
                end;
            
            else
                CostFunInterv_temp = CostFunIntervChar;
            end;
            
            for st = 1: length(CostFunInterv_temp);
                CostFunIntervNum(st) = str2num(CostFunInterv_temp{st});
            end;
            
            CostFunIntervNumTemp(csf,:)=CostFunIntervNum;
            CostFunInterv_temp = [];
            
            
            for csfn = 1 : 2 : size(CostFunIntervNum,2);
                indd1 = dsearchn(CostFun_Time(csf,:)',CostFunIntervNum(csfn));
                indd2 = dsearchn(CostFun_Time(csf,:)',CostFunIntervNum(csfn+1));
                varr_temp(csfn) = mean((rec_amp(csf,indd1:indd2)-mean(rec_amp(csf,indd1:indd2))).^2);
                mean_temp(csfn) = mean(rec_amp(csf,indd1:indd2));
            end;
            

            
            for csf2 = 1 : 2 : size(CostFunIntervNum,2);
                NumOfPoints = (CostFunIntervNum(csf2+1)-CostFunIntervNum(csf2))/dt;
                indd1 = dsearchn(CostFun_Time(csf,:)',CostFunIntervNum(csf2));
                indd2 = dsearchn(CostFun_Time(csf,:)',CostFunIntervNum(csf2+1));
%                 CostFun_Amp(csf,indd1:indd2) = CostFun_Value((csf2+1)/2)/(NumOfPoints*varr);
                CostFun_Amp(csf,indd1:indd2) = abs(CostFun_Value((csf2+1)/2)/mean_temp(csf2));
            end;
            
        else
            
            varr = Modeling_FunctionVar(rec_t(csf,:),rec_amp(csf,:));
%             CostFun_Amp(csf,1:end) = 1./(length(CostFun_Time)*varr);
            CostFun_Amp(csf,1:end) = abs(1/mean(rec_amp(csf,:)));
            
        end;
    end;
end;






output.cost_t = CostFun_Time;
output.cost_amp = CostFun_Amp;
