function [hi] = hamperl(dfft,wind)
    h = zeros(size(dfft));
    std_mult = 1;
    %figure; hold on
    %loglog(abs(dfft),'k')
    for i = 1:length(dfft)
        h(i) = dfft(i); 
        if i > floor(wind/2)
            if i > (length(dfft)-floor(wind/2))
                rr = dfft(end-wind:end);
            else
                rr = dfft(i-floor(wind/2):i+floor(wind/2));                
            end
            mr = median(rr);
            stdr = std(rr);

            if abs(abs(dfft(i)) - abs(mr)) > std_mult*stdr
                h(i) = mr;
            end
        end
    end
    %loglog(abs(h),'r')
    hi = ifft(h,'symmetric');
    %loglog(abs(h(1:131)),'k')
end