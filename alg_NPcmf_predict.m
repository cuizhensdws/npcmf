function [A,B]=alg_NPcmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_NPcmf_predict is a helper function of NPCMF that executes the alternating
%least squares method to obtain the solution. 
    
    K = size(A,2);
    alpha = 0.5;
    Sd = alpha*Sd + (1-alpha)*getGipKernel(Y);
    St = alpha*St + (1-alpha)*getGipKernel(Y');
 
	lambda_d_Sd = lambda_d*Sd;
    lambda_t_St = lambda_t*St;
    lambda_l_eye_K = lambda_l*eye(K);

  Sd(logical(eye(length(Sd)))) = 0;   % remove self-similarities
  [maxx, indx] = max(Sd);             % get nearest neighbor for each miRNA
   for i=1:length(Sd)
     Sd(i, :) = 0;                   % reset all similarities to 0...
     Sd(i, indx(i)) = maxx(i);       % except that of the nearest neighbor
   end
  St(logical(eye(length(St)))) = 0;   % remove self-similarities
  [maxx, indx] = max(St);             % get nearest neighbor for each disease
   for j=1:length(St)
      St(j, :) = 0;                   % reset all similarities to 0...
      St(j, indx(j)) = maxx(j);       % except that of the nearest neighbor
 end
 
    if nargin < 10
        AtA = A'*A;
        BtB = B'*B;
        for z=1:num_iter
            A = (Y*B + lambda_d_Sd*A)  / (BtB + lambda_l_eye_K + lambda_d*(AtA));
            AtA = A'*A;
            B = (Y'*A + lambda_t_St*B) / (AtA + lambda_l_eye_K + lambda_t*(BtB));
            BtB = B'*B;
        end
        
    else
        H = W .* Y;
        for z=1:num_iter
            A_old = A;
            HB_plus_lambda_d_Sd_A_old = H*B + lambda_d_Sd*A_old;
            lambda_l_eye_k_plus_lambda_d_A_oldt_A_old = lambda_l_eye_K + lambda_d*(A_old'*A_old);
            for a=1:size(A,1)
                A(a,:) = HB_plus_lambda_d_Sd_A_old(a,:) / (B'*diag(W(a,:))*B + lambda_l_eye_k_plus_lambda_d_A_oldt_A_old);
            end
            B_old = B;
            HtA_plus_lambda_t_St_B_old = H'*A + lambda_t_St*B_old;
            lambda_l_eye_k_plus_lambda_t_B_oldt_B_old = lambda_l_eye_K + lambda_t*(B_old'*B_old);
            for b=1:size(B,1)
                B(b,:) = HtA_plus_lambda_t_St_B_old(b,:) / (A'*diag(W(:,b))*A + lambda_l_eye_k_plus_lambda_t_B_oldt_B_old);
            end

        end
    end
    function Dij=L21f(Q)
 
   QQ=sqrt(sum(Q.*Q,1) + 1e-6);
  
  QQA=1./QQ;
  
   Dij = diag (QQA);
 
end
end