function makeprettyNoText()
% set some graphical attributes of the current axis
set(gcf,'color','w')

set(get(gca, 'XLabel'), 'FontSize', 6);
set(get(gca, 'YLabel'), 'FontSize', 6);
set(gca, 'FontSize', 7);

set(get(gca, 'Title'), 'FontSize', 9);

ch = get(gca, 'Children');
axis tight;
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
