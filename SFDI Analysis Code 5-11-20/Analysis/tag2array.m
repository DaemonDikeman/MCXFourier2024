function [ values ] = tag2array( xmlString, tag )
%tag2array Takes a string representation of an xml file and returns the
%numerical values within the specified tag.

%Find the tag in the file
tagStart = strfind(xmlString,['<',tag]);
tagEnd = strfind(xmlString, ['</',tag]);
%Warn if there are multiple matches
if length(tagStart)>1
    warning('String entered matches more than one tag in the file. The first will be used');
%Error if there are no matches
elseif isempty(tagStart)
    error('Tag not found');
end
%Get characters within the tags
tagString = xmlString(tagStart(1):tagEnd(1));
%Split along the square brackets
splitString= strsplit(tagString,{'<','>'});
%Ignore all non numeric characters
values=str2double(splitString(~isnan(str2double(splitString))));

end

