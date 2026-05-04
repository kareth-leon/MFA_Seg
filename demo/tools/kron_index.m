function out = kron_index(in)
    in(isnan(in))=true; % this solves the problem of inexisting neighbours at the field's boundaries
    out = ~in;
end    