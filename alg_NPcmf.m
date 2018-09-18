function y3=alg_NPcmf(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix)

    % get best parameters
    [k,lambda_l,lambda_d,lambda_t,num_iter] = alg_NPcmf_parest(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix);
    fprintf('k%g\t\t%g\t%g\t%g\t\t',k,lambda_l,lambda_d,lambda_t);

    % preprocessing Y
    if use_WKNKN
        Y = preprocess_WKNKN(Y,Sd,St,K,eta);
    end

    % initialize A & B
    [A,B] = initializer(Y,k);

    % predict
    test_ind = get_test_indices(Y,cv_setting,left_out);
    W = ones(size(Y));
    W(test_ind) = 0;
    [A,B] = alg_NPcmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W);

    % compute prediction matrix
    y3 = A*B';

    %--------------------------------------------------------------------

end