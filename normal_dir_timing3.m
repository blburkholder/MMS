function [n,x1] = normal_dir_timing3(fig_num,nin,jt)

    figure(fig_num);
    subplot(6,1,3)
    [x,~] = ginput(1);

    x1 = find(abs(jt - x) == min(abs(jt - x)));
    n = nin(x1,:,1);
end