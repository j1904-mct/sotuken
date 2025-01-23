HEADERS = vector3D.h particle.h bhtree.h
CCC = g++ -O3 -Wall $(PGCFLAGS)

BHtree: BHtree.cpp $(HEADERS)
	$(CCC) -o BHtree -DTEST BHtree.cpp
