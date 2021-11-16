function save_pdf_min_size(filename)
%SAVE_PDF Summary of this function goes here
%   Detailed explanation goes here
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,filename,'-dpdf','-r0')
end

