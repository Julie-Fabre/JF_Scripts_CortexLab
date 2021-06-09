
function makeprettyLarge()
% set some graphical attributes of the current axis

set(get(gca, 'XLabel'), 'FontSize', 25);
set(get(gca, 'YLabel'), 'FontSize', 25);
set(gca, 'FontSize', 17);

set(get(gca, 'Title'), 'FontSize', 30);

ch = get(gca, 'Children');
set(gcf,'color','w');
for c = 1:length(ch)
    thisChild = ch(c);
    if strcmp('line', get(thisChild, 'Type')) 
        if strcmp('.', get(thisChild, 'Marker'))
            set(thisChild, 'MarkerSize', 15);
        end
        if strcmp('-', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2.0);
        end
    end
end
