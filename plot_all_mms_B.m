function [] = plot_all_mms_B(jt,bx_int,by_int,bz_int,sc)
    for i = 1:sc
         plot(jt,bx_int(:,i),'r')
         plot(jt,by_int(:,i),'g')
         plot(jt,bz_int(:,i),'b')
%        plot(jt,by_int(:,i),'k')
%         plot(jt,bx_int(:,i)-mean(bx_int(:,i)),'r')
%         plot(jt,by_int(:,i)-mean(by_int(:,i)),'g')
%         plot(jt,bz_int(:,i)-mean(bz_int(:,i)),'b')
    end
end