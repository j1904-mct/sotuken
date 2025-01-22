#PGPLOT_DIR = /usr/local/pgplot
#PGINCLUDE =  -I$(PGPLOT_DIR)
#PGCFLAGS = -O2  $(PGINCLUDE)  
#PGCLIB     =  -L$(PGPLOT_DIR) -lcpgplot -lpgplot -lX11 -lUfor -lfor -lFutil -lm -lprof1 -lpdf
#PGCLIB     =  -L$(PGPLOT_DIR) -lcpgplot -lpgplot -lm  -L/usr/X11/lib -lX11 -lg2c
#PGLIB     =  -L$(PGPLOT_DIR)   -lpgplot -L/usr/X11/lib  -lX11
#PGCLIB     =  -L$(PGPLOT_DIR) -lcpgplot -lpgplot -lm  -L/usr/X11R6/lib -lX11 -lg2c
#PGLIB     =  -L$(PGPLOT_DIR)   -lpgplot -L/usr/X11R6/lib  -lX11
HEADERS = vector.h particle.h BHtree.h
CCC = g++ -Wall $(PGCFLAGS)

BHtreetest: BHtree.C $(HEADERS)
	$(CCC) -o BHtreetest -DTREETEST BHtree.C
BHtree: BHtree.C $(HEADERS)
	$(CCC) -o BHtree -DTEST BHtree.C
BHtree.o: BHtree.C $(HEADERS)
	$(CCC) -c  BHtree.C
treecode.o: treecode.C $(HEADERS)
	$(CCC) -c  treecode.C
treecode: treecode.o BHtree.o
	$(CCC) -o  treecode  treecode.o BHtree.o $(PGCLIB)
mk_plummer: mk_plummer.C
	$(CCC) -o  mk_plummer   mk_plummer.C
