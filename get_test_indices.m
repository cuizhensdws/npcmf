function test_ind=get_test_indices(Y,cv_setting,left_out)

    % 'left_out' is left-out interactions
    if strcmp(cv_setting,'cv_p')
        test_ind = left_out;
end