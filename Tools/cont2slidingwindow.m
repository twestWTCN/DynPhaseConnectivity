function X = cont2slidingwindow(Y,window,overlap)
for i = 1:size(Y,2);
    X(:,i)= mean(slideWindow(Y(:,i),window,overlap),1);
end