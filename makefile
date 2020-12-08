OBJS = D_update.o E_update.o H_update.o calc_fdtd.o main.o memory_allocate3d.o
HEADERS = main.h fdtd3d.h
OPTS = -O3 -std=c++1z

main: $(OBJS)
	g++ -o $@ $(OBJS)
%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)