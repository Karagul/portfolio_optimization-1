function ret = price2ret(price)
    [N, M] = size(price);
    ret = zeros(N, M);
    for i=1:N-1
        ret(i, :) = price(i, :) ./ price(i+1, :) - 1;
%         ret(i, :) = price(i, :) - price(i+1, :);
    end
end