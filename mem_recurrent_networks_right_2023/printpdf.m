function printpdf(h, outfilename, pos)
%PRINTPDF Print figure object 'h' nicely.

if nargin == 2
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos = get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    print('-dpdf',outfilename);
elseif nargin == 3
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[pos(1) pos(2) pos(3) pos(4)]);
    print('-dpdf',outfilename);
end

end

