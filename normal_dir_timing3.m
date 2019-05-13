function n = normal_dir_timing3(jt,b_x,b_y,b_z,x,y,z)

        px = zeros(3,1);
        py = zeros(3,1);
        pz = zeros(3,1);
        t = zeros(3,1);
        
        figure; hold on;

        for i = 1:3;
            x_comp = x(:,i);    
            y_comp = y(:,i);    
            z_comp = z(:,i);

            bx = b_x(:,i);
            by = b_y(:,i);
            bz = b_z(:,i);

            gradix = abs(bx(1:end-1)-bx(2:end));
            gradiy = abs(by(1:end-1)-by(2:end));
            gradiz = abs(bz(1:end-1)-bz(2:end));

            if max(gradix) >= max(gradiy) && max(gradix) >= max(gradiz)
                bbb = bx;
                gradi = gradix;
            elseif max(gradiy) >= max(gradix) && max(gradiy) >= max(gradiy)
                bbb = by;
                gradi = gradiy;
            else
                bbb = bz;
                gradi = gradiz;
            end

            plot(jt,bx,'r')
            plot(jt,by,'g')
            plot(jt,bz,'b')

            bbbb = bbb((gradi == max(gradi)));

            px(i) = interp1(jt,x_comp,(jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2);
            py(i) = interp1(jt,y_comp,(jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2);
            pz(i) = interp1(jt,z_comp,(jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2);
            t(i) = (jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2;
            scatter(t(i),(bbb(bbb == bbbb)+bbb(find(bbb == bbbb)+1))/2);
        end

        meep = input('do you fucking accept?');
        if meep 
            r12 = [px(1) - px(2),py(1) - py(2),pz(1) - pz(2)];
            r13 = [px(1) - px(3),py(1) - py(3),pz(1) - pz(3)];
            t12 = t(1) - t(2);
            t13 = t(1) - t(3);
        else
            splooge = 0;
            while splooge == 0
                peak = ginput(2);
                t11 = find(abs(jt-peak(1,1)) == min(abs(jt-peak(1,1))));
                t22 = find(abs(jt-peak(2,1)) == min(abs(jt-peak(2,1))));
                letter = input('x or y or z');
                if letter == 1
                    b1 = b_x(t11:t22,1); b2 = b_x(t11:t22,2); b3 = b_x(t11:t22,3);
                elseif letter == 2
                    b1 = b_y(t11:t22,1); b2 = b_y(t11:t22,2); b3 = b_y(t11:t22,3);
                else
                    b1 = b_z(t11:t22,1); b2 = b_z(t11:t22,2); b3 = b_z(t11:t22,3); 
                end
                extrema = input('min or max');
                t = jt(t11:t22);
                if extrema == 1
                    t1 = t(b1 == min(b1)); t2 = t(b2 == min(b2)); t3 = t(b3 == min(b3));
                    px1 = interp1(jt,x(:,1),t1); py1 = interp1(jt,y(:,1),t1); pz1 = interp1(jt,z(:,1),t1);
                    px2 = interp1(jt,x(:,2),t2); py2 = interp1(jt,y(:,2),t2); pz2 = interp1(jt,z(:,2),t2);
                    px3 = interp1(jt,x(:,3),t3); py3 = interp1(jt,y(:,3),t3); pz3 = interp1(jt,z(:,3),t3);
                    scatter(t1,min(b1)); scatter(t2,min(b2)); scatter(t3, min(b3));
                else
                    t1 = t(b1 == max(b1)); t2 = t(b2 == max(b2)); t3 = t(b3 == max(b3));
                    px1 = interp1(jt,x(:,1),t1); py1 = interp1(jt,y(:,1),t1); pz1 = interp1(jt,z(:,1),t1);
                    px2 = interp1(jt,x(:,2),t2); py2 = interp1(jt,y(:,2),t2); pz2 = interp1(jt,z(:,2),t2);
                    px3 = interp1(jt,x(:,3),t3); py3 = interp1(jt,y(:,3),t3); pz3 = interp1(jt,z(:,3),t3);
                    scatter(t1,max(b1)); scatter(t2,max(b2)); scatter(t3, max(b3));
                end

                r12 = [px1 - px2,py1 - py2,pz1 - pz2];
                r13 = [px1 - px3,py1 - py3,pz1 - pz3];
                t12 = t1 - t2;
                t13 = t1 - t3;
                splooge = input('do you fucking accept?');
            end
        end
        nt = [r12;r13]\[t12;t13];
        n = nt/sqrt(nt(1)^2+nt(2)^2+nt(3)^2);   
end