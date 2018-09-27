function [labels] = create_labels(storage_names,attr)
tmp = [storage_names storage_names];
for m = 1:length(tmp)
    if m <= length(storage_names)
        tmp{m} = [attr '_{' tmp{m} '}'];
    else
        tmp{m} = [attr '_{' tmp{m} ', SWMM}'];
    end
end

labels = tmp;