function Rsq = sqres(Y,yCalc)
if size(Y,1)>1
    for i = 1:size(Y,1)
        Y(i,isnan(Y(i,:))) = mean(Y(:,isnan(Y(i,:))));
        Rsq(i) = sum((Y(i,:) - yCalc).^2);
    end
else
    Rsq = (Y - yCalc).^2;
end
