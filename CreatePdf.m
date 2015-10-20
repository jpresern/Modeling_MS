%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   prepares .ps and .pdf files (currently commented out)

function CreatePdf(ff,filenamee,pathh)

for a = 1:length(ff)
    set(ff(a),'PaperType','A4','PaperUnits','normalized','PaperPosition',[0 0 1 1]);
    print('-dps2c',ff(a),horzcat(pathh,filenamee),'-append');
end;


% if strcmp (computer,'PCWIN64')
%     ps2pdf('psfile',horzcat(pathh,filenamee,'.ps'),'pdffile',...
%         horzcat(pathh,filenamee,'.pdf'),'gspapersize', 'a4','deletepsfile',1);
% else
% %     !ps2pdf(horzcat(pathh,filenamee,'.ps'),horzcat(pathh,filenamee,'.pdf')
% end;


% for a = 1:length(ff)
% %     set(ff(a),'PaperType','A4','PaperUnits','normalized','PaperPosition',[0 0 1 1]);
%     print('-dpdf',ff(a),horzcat(pathh,filenamee),'-append');
% end;
