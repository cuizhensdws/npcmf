function AUPR=calculate_aupr(diseases,predicts)

	if nargin > 1
		[~,i] = sort(predicts,'descend');
		diseases = diseases(i);
	end
	
	cumsums = cumsum(diseases)./reshape(1:numel(diseases),size(diseases));
	AUPR = sum(cumsums(~~diseases));
	pos = sum(diseases);
	AUPR = AUPR / pos;
end