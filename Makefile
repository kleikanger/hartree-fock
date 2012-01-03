HFmatrix: HFmatrix.cpp HFmatrix.h lib/lib.cpp lib/lib.h
	c++ -o HFmatrix.out HFmatrix.cpp lib/lib.cpp #-Wall -O3
		
run: HFmatrix 
	make HFmatrix	
	./HFmatrix.out 

clear:
	rm *.out
	make HFmatrix
	echo 'Har fjaena alle .out filer foer make HFmatrix'		
