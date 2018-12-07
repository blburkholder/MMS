function [v] = boundary_normal(i)
    [~,~,~,~,bx_comp,by_comp,bz_comp,~,~] = MMS_fgm(i);
    for s = 1:4;
        bx = bx_comp(:,s);
        by = by_comp(:,s);
        bz = bz_comp(:,s);

    %     avg = 100;
    %     ind = floor(length(bx)/avg)*avg;
    %     bxm = mean(reshape(bx(1:ind),[avg,length(bx(1:ind))/avg])',2);
    %     bym = mean(reshape(by(1:ind),[avg,length(by(1:ind))/avg])',2);
    %     bzm = mean(reshape(bz(1:ind),[avg,length(bz(1:ind))/avg])',2);
    %     tm = mean(reshape(t(1:ind),[avg,length(t(1:ind))/avg])',2);

        figure(1)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        plot(bx)
        hold on
        plot(by)
        plot(bz)

        [x,~] = ginput(1);
        x = round(x);
        ymax = max([max(bx) max(by) max(bz)]);
        ymin = min([min(bx) min(by) min(bz)]); 
        plot([x,x],[ymin,ymax])
        if s == 1
            [l,~] = ginput(1);
            l = round(l);
        end
        dif = abs(x-l);
        %plot([x-dif,x-dif],[ymin,ymax])
        %plot([x+dif,x+dif],[ymin,ymax])
        close(1)

        [v,~] = normal_dir_var(bx(x-dif:x+dif),...
            by(x-dif:x+dif),bz(x-dif:x+dif),1,1,1)
end