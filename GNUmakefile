# ===========================================================================
#  Makefile                                        Chun Shen Apr. 16, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -g -Wall 

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) -L/sw/lib -lgsl -lgslcblas
INCLUDE         =       -I/sw/include
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	calEmissionrates.e
endif

SRC		=	main.cpp Arsenal.cpp ParameterReader.cpp \
                  Collinear_Kernel.cpp gauss_quadrature.cpp \
                  Physicalconstants.cpp

INC		= 	Arsenal.h ParameterReader.h \
                  Collinear_Kernel.h gauss_quadrature.h \
                  Physicalconstants.h Stopwatch.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(LDFLAGS) $(INCLUDE) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
./main.cpp : ParameterReader.h Collinear_Kernel.h Stopwatch.h Arsenal.h
./Collinear_Kernel.cpp : Collinear_Kernel.h ParameterReader.h Arsenal.h gauss_quadrature.h 
./ParameterReader.cpp : ParameterReader.h Arsenal.h
./Arsenal.cpp : Arsenal.h
./Physicalconstants.cpp : Physicalconstants.h
./gauss_quadrature.cpp : gauss_quadrature.h
