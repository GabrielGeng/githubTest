function [ ROW_idx, COL_idx ] = find_2d_peaks( Spec, max_peak_no )

% Spec = [1 2 3 2 1; 2 3 8 2 2; 2 9 3 1 2; 1 2 5 2 1];

no_row = size(Spec,1);
no_col = size(Spec,2);

row_diff_1 = Spec(2:no_row-1,:) - Spec(1:no_row-2,:);
row_diff_2 = Spec(2:no_row-1,:) - Spec(3:no_row,:);
[row_dif_idx_r, row_dif_idx_c] = find(row_diff_1>0 & row_diff_2>0);
row_ID = [row_dif_idx_r+1, row_dif_idx_c];

col_diff_1 = Spec(:,2:no_col-1) - Spec(:,1:no_col-2);
col_diff_2 = Spec(:,2:no_col-1) - Spec(:,3:no_col);

[col_dif_idx_r, col_dif_idx_c] = find(col_diff_1>0 & col_diff_2>0);
col_ID = [col_dif_idx_r, col_dif_idx_c+1];

ID = [];
for i=1:size(row_ID,1)
    for j = 1:size(col_ID,1)
        if row_ID(i,:) ==  col_ID(j,:);
            ID = [ID;row_ID(i,:)];
        end
    end
end

[~,idx] = sort(diag(Spec( ID(:,1),ID(:,2) )),'descend');
if length(idx)>max_peak_no
    idx = idx(1:max_peak_no);
end
ID = ID(idx,:);

ROW_idx = ID(:,1);
COL_idx = ID(:,2);

end

