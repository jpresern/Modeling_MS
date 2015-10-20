%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, output] = Draw_8A(dt,stimAmp, stimTime,...
                                Diagram_x, Diagram_y,...
                                var, var_names,...
                                cmap1, cmap2, fn)
                                
f = figure;

s(2) = axes('OuterPosition', [0 0.4 1 0.4]); 

hold on;
b = 10;
for a = 1:10
    [~,~,results]=Modeling_DRG_TCM_Engine(horzcat('DRG_TCM_Model_mh','_Report'),...
        stimTime(b,:),stimAmp(b,:),...
        var,var_names,dt);
    p1 = plot(results.t,results.g,'-','Color',cmap2(b,:),'LineWidth',2);
    output.Fig8A.model.Imax(b) = min(results.g);
    output.Fig8A.model.tIMax(b) = results.t(min(results.g)==results.g);
    line ([output.Fig8A.model.tIMax(b) output.Fig8A.model.tIMax(b)],[output.Fig8A.model.Imax(b) 0],'Color','r','LineWidth',1.5);
    b = b-1;
end;
plot (output.Fig8A.model.tIMax, output.Fig8A.model.Imax, '*r','MarkerSize',5);
hold off;
set(gca,'Ylim',[-1.1 0]);
xlim([0 220]);
xlabel ('time [ms]');
ylabel ('I/I_{max}');
grid on;

s(1) = axes('OuterPosition', [0 0.8 1 0.2]); 
set(gca,'XTickLabel',[]);
title(horzcat(fn,': ','Fig8A: Stimulus speed encoding in DRG neuron'));
for val2 = 1:10
    x0 = 0;
    y0 = 0;
    for val1 = 1:5
        x1 = x0 + stimTime(val2,val1);
        y1 = y0 + stimAmp(val2,val1);
        line([x0 x1],[y0 y1],'Color',cmap2(val2,:),'LineWidth',1.5);
        x0 = x1;
        y0 = y1;
    end;
    line ([output.Fig8A.model.tIMax(val2) output.Fig8A.model.tIMax(val2)],[8.5 0],'Color','r','LineWidth',1.5);
end;
xlim ([0 220]);
ylabel ('[stim amplitude]');
set (gca,'YTick',[8.5]);
set (gca,'YTickLabel',[8.5]);


s(3) = axes('OuterPosition', [0 0.0 1 0.3]); 

velocity = (stimAmp(:,2)./stimTime(:,2))*1000;
hold on;
plot (velocity,output.Fig8A.model.Imax./min(output.Fig8A.model.Imax), '*r','LineWidth',1.5,'MarkerSize',10);
plot (Diagram_x, Diagram_y, 's', 'LineWidth',1.5, 'MarkerSize',5,'Color',cmap1(14,:));
set(s(3), 'Xdir','reverse')
xlabel ('stimulus up-ramp velocity [\mum/s]');
ylabel ('I/I_{max}');
grid on;
title(horzcat(fn,': ','Fig8C: Stimulus speed encoding in DRG neuron'));

