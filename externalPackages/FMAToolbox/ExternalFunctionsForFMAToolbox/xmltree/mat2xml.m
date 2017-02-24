function varargout = mat2xml(matfile, xmlfile)
%MAT2XML Convert a MAT-file into an XML file
%  MAT2XML(MATFILE, XMLFILE) loads MAT-file MATFILE and converts it 
%  into an XML file saved in XMLFILE.
%  T = MAT2XML(MATFILE) also converts the MAT-file MATFILE to XML
%  but returns it instead as an XMLTree object T.
%
%  See also XML2MAT, XMLTREE, SAVE.

%  Copyright 2003 Guillaume Flandin. 
%  Revision: 1.0 $  $Date: 2003/07/01 21:43 $

error(nargchk(1,2,nargin));

s = load(matfile);
names = fieldnames(s);

t = xmltree(['<matfile name="' matfile '"/>']);

for i=1:length(names)
	[t, uid] = add(t,root(t),'element',names{i});
	t = attributes(t,'add',uid,'type',class(getfield(s,names{i})));
	t = attributes(t,'add',uid,'size',sub_num2str(size(getfield(s,names{i}))));
	t = sub_var2xml(t,getfield(s,names{i}),uid);
end

if nargin == 2
	save(t,xmlfile);
else
	varargout{1} = t;
end

%=======================================================================
function t = sub_var2xml(t,v,uid)

	switch class(v)
		case 'double'
			t = add(t,uid,'chardata',sub_num2str(v));
		case 'sparse'
			% TODO % better names for sparse elements
			[i,j,s] = find(v);
			[t, uid2] = add(t,uid,'element','row');
			t = attributes(t,'add',uid2,'size',sub_num2str(size(i)));
			t = add(t,uid2,'chardata',sub_num2str(i));
			[t, uid2] = add(t,uid,'element','col');
			t = attributes(t,'add',uid2,'size',sub_num2str(size(j)));
			t = add(t,uid2,'chardata',sub_num2str(j));
			[t, uid2] = add(t,uid,'element','val');
			t = attributes(t,'add',uid2,'size',sub_num2str(size(s)));
			t = add(t,uid2,'chardata',sub_num2str(s));
		case 'struct'
			names = fieldnames(v);
			for j=1:prod(size(v))
				for i=1:length(names)
					[t, uid2] = add(t,uid,'element',names{i});
					t = attributes(t,'add',uid2,'index',num2str(j));
					t = attributes(t,'add',uid2,'type',...
						class(getfield(v(j),names{i})));
					t = attributes(t,'add',uid2,'size', ...
						sub_num2str(size(getfield(v(j),names{i}))));
					t = sub_var2xml(t,getfield(v(j),names{i}),uid2);
				end
			end
		case 'cell'
			for i=1:prod(size(v))
				[t, uid2] = add(t,uid,'element','cell'); 
				% TODO % special handling of cellstr ?
				t = attributes(t,'add',uid2,'index',num2str(i));
				t = attributes(t,'add',uid2,'type',class(v{i}));
				t = attributes(t,'add',uid2,'size',sub_num2str(size(v{i})));
				t = sub_var2xml(t,v{i},uid2);
			end
		case 'char'
			% TODO % char values should be in CData
			t = add(t,uid,'chardata',v);
		case {'int8','uint8','int16','uint16','int32','uint32'}
			[t, uid] = add(t,uid,'element',class(v));
			% TODO % Handle integer formats (cannot use sprintf or num2str)
		otherwise
			if ismember('serialize',methods(class(v)))
				% TODO % is CData necessary for class output ?
				t = add(t,uid,'cdata',serialize(v));
			else
				warning(sprintf(...
				'[MAT2XML] Cannot convert from %s to XML.',class(v)));
			end
	end
	
%=======================================================================
function s = sub_num2str(n)
	% TODO % use format ?
	if isempty(n)
		s = '[]';
	else
		s = ['[' sprintf('%g ',n(1:end-1))];
		s = [s num2str(n(end)) ']'];
	end
