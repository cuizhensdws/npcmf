function krnl=getGipKernel(y)
%getGipKernel computes the GIP (Gaussian Interaction Profile) kernel from
%the interaction matrix y
%
% krnl = getGipKernel(y)

	krnl = y*y';
	krnl = krnl / mean(diag(krnl));
	krnl = exp(-kernel_to_distance(krnl));

end


function d=kernel_to_distance(k)
	di = diag(k);
	d = repmat(di,1,length(k)) + repmat(di',length(k),1) - 2 * k;
end