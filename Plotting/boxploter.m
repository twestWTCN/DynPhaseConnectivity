function a = boxploter(x,y,cmap,R,ylab)
pv = [];
% subplot(1,2,2)
x = x(:); x(x==0) = [];
y = y(:); y(y==0) = [];
G = [repmat(1,size(x)); repmat(2,size(y))];
[pv h] = ranksum(x,y);
boxplot([x;y],G,'labels',{R.condname},...
    'BoxStyle','filled','Widths',0.8)
%%
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
boxes(1).Color = cmap;
boxes(2).Color = cmap.*0.7;
set(boxes,'linewidth',35); % Set width
idx=strcmpi(t,'Whisker');  % Find Box objects
whisker=a(idx);          % Get the children you need
set(whisker,'LineWidth',2); % Set width
idx=strcmpi(t,'Median');  % Find Box objects
whisker=a(idx);          % Get the children you need
set(whisker,'LineWidth',2); % Set width
grid on
ylabel(ylab); ylim([-0.4 0.4]);
set(gcf,'Position',[678         668        1164         310])
