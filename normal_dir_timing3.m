function [n,coond] = normal_dir_timing3(jt,b_x,b_y,b_z,x,y,z)

        px = zeros(3,1);
        py = zeros(3,1);
        pz = zeros(3,1);
        t = zeros(3,1);
        
        %figure; hold on;

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

            %plot(jt,bbb)

            bbbb = bbb((gradi == max(gradi)));

            px(i) = interp1(jt,x_comp,(jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2);
            py(i) = interp1(jt,y_comp,(jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2);
            pz(i) = interp1(jt,z_comp,(jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2);
            t(i) = (jt(bbb == bbbb)+jt(find(bbb == bbbb)+1))/2;
            %scatter(t(i),(bbb(bbb == bbbb)+bbb(find(bbb == bbbb)+1))/2);
        end

        r12 = [px(1) - px(2),py(1) - py(2),pz(1) - pz(2)];
        r13 = [px(1) - px(3),py(1) - py(3),pz(1) - pz(3)];
        t12 = t(1) - t(2);
        t13 = t(1) - t(3);

        nt = [r12;r13]\[t12;t13];
        n = nt/sqrt(nt(1)^2+nt(2)^2+nt(3)^2);   
        coond = cond([r12;r13]);
end