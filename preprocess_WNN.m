function Y=preprocess_WNN(Y,Sd,St,cv_setting,eta)
%preprocess_WNN preprocesses the interaction matrix Y using WNN (weighted
%nearest neighbors)
%
% Y = preprocess_WNN(Y,Sd,St,cv_setting,eta)
%
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise row similarities matrix
%  St:          pairwise column similarities matrix
%  cv_setting:  cross validation setting ( 'cv_p')
%  eta:         decay rate (used by WNN)
%
% OUTPUT:
%  Y:           preprocessed interaction matrix
%
%

    y = Y;
    if strcmp(cv_setting,'cv_d')
        empty_rows = find(any(Y,2) == 0);   % get indices of empty rows
        w = eta .^ (0:length(Sd)-1);
        w(w < 10^-4) = [];
        k = length(w);
        for r=1:length(empty_rows)
            i = empty_rows(r);  % i = current empty row
            miRNA_sim = Sd(i,:); % get similarities of miRNA i to other miRNAs
            miRNA_sim(i) = 0;    % set self-similarity to ZERO
            [~,indx]=sort(miRNA_sim,'descend');  % sort descendingly
            indx = indx(1:k);
            y(i,:) = w * Y(indx,:);     % multiply sorted similarities by decreasing decay values
        end
    elseif strcmp(cv_setting,'cv_t')
        empty_cols = find(any(Y) == 0);   % get indices of empty columns
        w = eta .^ (0:length(St)-1);
        w(w < 10^-4) = [];
        k = length(w);
        for c=1:length(empty_cols)
            j = empty_cols(c);      % j = current empty column
            disease_sim = St(j,:);   % get similarities of disease j to other diseases
            disease_sim(j) = 0;      % set self-similarity to ZERO
            [~,indx]=sort(disease_sim,'descend'); % sort descendingly
            indx = indx(1:k);
            y(:,j) = Y(:,indx) * w';    % multiply sorted similarities by decreasing decay values
        end
    end
    Y = y;

end