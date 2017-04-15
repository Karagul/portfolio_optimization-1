function r_size = getMin(A, B, C)
    r_size = length(A);
    if(length(B) < r_size)
        r_size = length(B);
    end
     if(length(C) < r_size)
        r_size = length(C);
     end
end