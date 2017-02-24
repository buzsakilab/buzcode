To compile mym, type:

>> mex -I/usr/include/mysql -lmysqlclient -lz mym.cpp

If matlab complains that the resulting mex-file is not valid, you may need
to force mex to use a specific (usually older) version of gcc, e.g.

>> mex -I/usr/include -I/usr/include/mysql -lmysqlclient -lz CXX=g++-3.3 mym.cpp

