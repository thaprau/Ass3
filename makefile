jtest: johantest.c
	gcc -O3 johantest.c graphics/graphics.c graphics/graphics.h -o johantest -lpthread -lm -lX11
test: johantest.c
	gcc -O3 test.c graphics/graphics.c graphics/graphics.h -o test -lpthread -lm -lX11
run: 
	./test 3000 input_data/ellipse_N_03000.gal 100 0.00001 0
jrun: 
	./johantest 3000 input_data/ellipse_N_03000.gal 100 0.00001 0
clean:
	rm test
runbir:
	./biran 2 input_data/sun_and_planet_N_2.gal 1000 0.00001 1
run_test:
	./test 2 result.gal 1000 0.00001 1e