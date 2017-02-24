

A = dictArray({ { 'abc', 42 }, 
		{ 'ciccio', pi }, 
		{ 'never', 0 }
		{ 'qqqq', 'romeo' } } );

A

A{'ciccio'}

A{'lollo'} = 'pinco';

B = dictArray({ { 'minchia', 'sonoio' }, 
		{ 'tusorella', 69 } 
		{ 'tonino', 'carino' } } );
B

A = update(A, B);

A

disp('del');
A = del(A, 'ciccio');

A

has_key(A, 'minchia')

c = get(A, 'minchia', 'maddeche')
iscell(c)
get(A, 'minchione', 'maddeche')

len(A)
len(B)

keys(A)

values(A)

items(A)