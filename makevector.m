function vec = makevector( inpt, makeLen )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


if isscalar(inpt)
    vec = repmat(inpt,[1,makeLen]);
elseif length(inpt) ~= makeLen
    warning('A parameter input did not match the length of the number of input files.\nOnly using first element.');
    vec = repmat(inpt(1),[1,makeLen]);
else
    vec = inpt;
end


end

