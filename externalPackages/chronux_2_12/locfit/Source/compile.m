if ispc()
	mex -v mexlf.c liblfev.c liblocf.c libmut_win.c libtube.c mlfut.c 
	mex -v mexpp.c liblfev.c liblocf.c libmut_win.c libtube.c mlfut.c 
else
	mex -v mexlf.c liblfev.c liblocf.c libmut.c libtube.c mlfut.c 
	mex -v mexpp.c liblfev.c liblocf.c libmut.c libtube.c mlfut.c 
end