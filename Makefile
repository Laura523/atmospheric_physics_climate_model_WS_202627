convection: convection_main.c convection_functions.c convection_functions.h
	clang convection_main.c convection_functions.c -o convection_main -lm