
function makeprettyLarge()
% set some graphical attributes of the current axis

set(get(gca, 'XLabel'), 'FontSize', 30);
set(get(gca, 'YLabel'), 'FontSize', 30);
set(gca, 'FontSize', 25);

set(get(gca, 'Title'), 'FontSize', 40);

ch = get(gca, 'Children');
set(gcf,'color','w');
for c = 1:length(ch)
    thisChild = ch(c);

    
    if strcmp('line', get(thisChild, 'Type')) 
        if strcmp('.', get(thisChild, 'Marker'))
            set(thisChild, 'MarkerSize', 40);
        end
        if strcmp('-', get(thisChild, 'LineStyle'))
            set(thisChild, 'LineWidth', 2.0);
        end
    end
end
