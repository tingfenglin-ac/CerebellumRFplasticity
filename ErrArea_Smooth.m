function ErrArea_Smooth(X,Y,ERR,LineColor)

err=fill([X;flipud(X)],[Y-ERR;flipud(Y+ERR)],LineColor(1:3),'linestyle','none');
% set(gco,'linestyle','none')
% Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
if length(LineColor)==4
    coloralpha=LineColor(4);
else
    coloralpha=.2;
end
set(err,'facealpha',coloralpha,'marker','none')
set(get(get(err,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

end