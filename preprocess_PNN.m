function S=preprocess_PNN(S,p)
% preprocess_PNN sparsifies the similarity matrix S by keeping, for each
% miRNA/Disease, the p nearest neighbors and discarding the rest.

    NN_mat = zeros(size(S));

    % for each miRNA/Disease...
    for j=1:length(NN_mat)
        row = S(j,:);                           % get row corresponding to current miRNA/Disease
        row(j) = 0;                             % ignore self-similarity
        [~,indx] = sort(row,'descend');         % sort similarities descendingly
        indx = indx(1:p);                       % keep p NNs
        NN_mat(j,indx) = S(j,indx);             % keep similarities to p NNs
        NN_mat(j,j) = S(j,j);                   % also keep the self-similarity (typically 1)
    end

    % symmetrize the modified similarity matrix
    S = (NN_mat+NN_mat')/2;

end