function o = optionmerge(varargin)
%
% o = optionmerge(o1,o2,...)
%
%    Merges option sets into one. Each option set is a structure with
%    named fields. The merge is again such a structure. If the same
%    parameter is set multiple times the last occurence counts. Since
%    each function must recognize its own parameters and disregard the
%    others it is a good practice to prepend a function identifier to
%    each parameter separated by an underscore (_).
%
%  o = optionmerge(o1,'o2',...)
%
%    The 'o2' is the name of the option variable. If it exist, it is
%    dereferenced, otherwise it is ignored. This help avoid the
%    exists('opt','var') test prior to calling optionmerge. Do not use this
%    option if you are chasing CPU time.
%
% See also: Options

% (c) Radim Sara (sara@cmp.felk.cvut.cz) FEE CTU Prague, 24 Jan 03

% Revisions:
%  30 Dec 12 : works with variable names.

 o = varargin{1};
 if ~isstruct(o); error 'Parameter set is not a structure'; end
 
 for i = 2:length(varargin)
  oi = varargin{i};
  if ~isempty(oi) && ~isstruct(oi) && ~ischar(oi)
   error([mfilename, ':Arguments'], 'Argument is neither a structure nor a string')
  end
   
  % if only name is passed, get the value, if the var exists
  if ischar(oi)
   ex = evalin('caller', sprintf('exist(''%s'',''var'')',oi));
   if ex
    oi = evalin('caller','opt');
   else
    continue
   end
  end
   
  if ~isempty(oi)
   for fld = fieldnames(oi)'
    %o = setfield(o,fld{1},getfield(oi,fld{1}));
    o.(fld{1}) = oi.(fld{1}); % this is about 3x faster
   end
  end
 end

end