function [n,x1,x2] = normal_dir_timing4(fig_num,n,jt,b_x,b_y,b_z,x,y,z)

    figure(fig_num);
    legend('B_x','B_y','B_z')
    [xx,~] = ginput(2);

    ccc = input('x (1) or y (2) or z (3) or MVA (4)');
    if ccc == 4 
        [n,x1] = normal_dir_timing3(300,n,jt);
        x2 = x1;
    else
        ddd = input('max (1) or min (2) or steepest gradient (3) ');

        clf
        px = zeros(4,1);
        py = zeros(4,1);
        pz = zeros(4,1);
        t = zeros(4,1);
        for i = 1:4;
            Bt = jt;
            pos_time = jt;

            x_comp = x(:,i);    
            y_comp = y(:,i);    
            z_comp = z(:,i);

            bx = b_x(:,i);
            by = b_y(:,i);
            bz = b_z(:,i);

            x1 = find(abs(Bt - xx(1)) == min(abs(Bt - xx(1))));
            x2 = find(abs(Bt - xx(2)) == min(abs(Bt - xx(2))));

            if ccc == 1
                bbb = bx(x1:x2);
                title('B_x')
            elseif ccc == 2
                bbb = by(x1:x2);
                title('B_y')
            else
                bbb = bz(x1:x2);
                title('B_z')
            end
            Btt = Bt(x1:x2);
            plot(Btt,bbb)
            hold on
            if ddd == 1
                bbbb = max(bbb);
            elseif ddd == 2
                bbbb = min(bbb);
            else
                gradi = abs(bbb(1:end-1)-bbb(2:end));
                bbbb = bbb((gradi == max(gradi)));
            end
            scatter(Btt(bbb == bbbb),bbbb)
            px(i) = interp1(pos_time,x_comp,Btt(bbb == bbbb));
            py(i) = interp1(pos_time,y_comp,Btt(bbb == bbbb));
            pz(i) = interp1(pos_time,z_comp,Btt(bbb == bbbb));
            t(i) = Btt(bbb == bbbb);
        end
        ax = gca;
        ax.XTick = [Btt(1),Btt(end)];
        ax.XTickLabel = datestr(todatenum([cdfepoch(Btt(1)),cdfepoch(Btt(end))]));

        r12 = [px(1) - px(2),py(1) - py(2),pz(1) - pz(2)];
        r13 = [px(1) - px(3),py(1) - py(3),pz(1) - pz(3)];
        r14 = [px(1) - px(4),py(1) - py(4),pz(1) - pz(4)];
        t12 = t(1) - t(2);
        t13 = t(1) - t(3);
        t14 = t(1) - t(4);

        nt = [r12;r13;r14]\[t12;t13;t14];
        n = nt/sqrt(nt(1)^2+nt(2)^2+nt(3)^2);   
    end