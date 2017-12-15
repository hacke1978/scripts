function status = check_to_cat(toCat)
% check if it is cattible
status = cellfun(@(x) (length(x)), {toCat});