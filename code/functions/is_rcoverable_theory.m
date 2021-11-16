function is_it = is_rcoverable_theory(w,k,d,range_sols)
    % determine if the dynamics is recoverable based on theory of
    % intersection between functions
    % d is the value of exciting, range_sols is the range in which all
    % states should be (x_high, x_low)
    dc = find_dc_by_wk_theory_intersection(w,k,range_sols);
    is_it = d>=dc;       
end
