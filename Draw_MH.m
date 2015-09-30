function fig  = Draw_MH(stimAmp, pred_m0_x, pred_m0_y,...
                            pred_h0_x, pred_h0_y)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;

s(1) = axes('OuterPosition', [0 0.0 1 1]); 

hold on;

p2 = plot(pred_m0_x, pred_m0_y,'-b','LineWidth',2);
p4 = plot(pred_h0_x, pred_h0_y,'-r','LineWidth',2);
set(gca,'XTick',stimAmp(:,2));

hold off;

xlabel('x [\mum]');
ylabel('rel. value');
title('Activation and inactivation parameters');
legend([p2,p4],'predicted m0','predicted h0','Location','southeast');

