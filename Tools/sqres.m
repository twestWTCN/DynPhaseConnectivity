function Rsq = sqres(Y,yCalc)
if size(Y,1)>1
    for i = 1:size(Y,1)
        Y1 = Y(i,:);
        Rsq(i) = sum((Y1(~isnan(Y1)) - yCalc(~isnan(Y1))).^2);
    end
else
    Rsq = (Y - yCalc).^2;
end
