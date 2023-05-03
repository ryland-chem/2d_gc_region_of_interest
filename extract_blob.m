function [subtensor] = extract_blob(tensor)

    % Find linear indices of non-zero values
    linear_indices = find(tensor(:));
    
    % Convert linear indices to 3D subscripts
    [i,j,k] = ind2sub(size(tensor), linear_indices);
    
    % Find boundaries of non-rectangular blob
    i_min = min(i);
    i_max = max(i);
    j_min = min(j);
    j_max = max(j);
    k_min = min(k);
    k_max = max(k);
    
    % Extract subregion of tensor containing non-zero values
    subtensor = tensor(i_min:i_max, j_min:j_max, k_min:k_max);

end