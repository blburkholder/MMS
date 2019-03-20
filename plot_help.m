function [] = plot_help(tit,miin,maax)
    title(tit);
    ax1 = gca;
    ax1.XTick = [];
    xlim([miin,maax]);
end