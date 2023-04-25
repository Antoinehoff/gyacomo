function [] = top_title(title)
% Home made function to use suptitle or sgtitle
try % valid for matlab 2016
    suptitle(title)
catch % valid for matlab 2023
    sgtitle(title)
end

end