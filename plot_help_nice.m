function [] = plot_help_nice(tit,miin,maax)
    title(tit);
    ax2 = gca;
    ax2.XTick = [miin,maax];
    ax2.XTickLabel = datestr(todatenum([cdfepoch(miin),cdfepoch(maax)]));
    xlim([miin,maax]);
end