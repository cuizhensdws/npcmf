function AUC = calculate_auc(diseases,predicts)
[~,i] = sort(predicts,'descend');
roc_y = diseases(i);
stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
AUC=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));
plot(stack_x,stack_y);
xlabel('1-Sp');
ylabel('Sn');
end