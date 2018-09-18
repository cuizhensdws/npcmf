function folds=get_folds(Y,cv_setting,nr_fold,left_out)

    % get training set
    if strcmp(cv_setting,'cv_p')        % 'left_out' is left-out pairs
        trainingSet = 1:numel(Y);
        trainingSet(left_out) = [];
    end
    len = length(trainingSet);
    rand_ind = randperm(len);
    rand_ind = trainingSet(rand_ind);

    folds = cell(nr_fold);
    for i=1:nr_fold
        if strcmp(cv_setting,'cv_p')
            folds{i} = rand_ind((floor((i-1)*len/nr_fold)+1:floor(i*len/nr_fold))');
        end
    end
    
end