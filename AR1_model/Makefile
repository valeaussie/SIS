all:
	g++ -o AR1 -pedantic -std=c++11 -Wall -Werror AR1.cpp
	g++ -c -pedantic -std=c++11 -Wall AR1_model.cpp
	g++ -c -pedantic -std=c++11 -Wall SIS_corrections.cpp
	g++ AR1_model.o SIS_corrections.o -o AR1_simulation
	g++ -c -pedantic -std=c++11 -Wall AR1_interp.cpp
	g++ AR1_model.o AR1_interp.o -o AR1_interp
