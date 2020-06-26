function [ output ] = escape( input )
% ESCAPE Escape special characters in string (i.e., for figure titles)

    output = strrep(input, '_', '\_');
end

