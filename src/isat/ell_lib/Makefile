# Author: Varun Hiremath <vh63@cornell.edu>
# Last Modified:  Thu, 12 Nov 2009 16:27:00 -0500

include ../build-files/vars.mk

LIB_NAME = libell

# ------------- OBJECTS --------------------

include SOURCES.mk

OBJS=$(FILESf90:.f90=.$(OBJ)) $(FILESf:.f=.$(OBJ))

# ------------- TARGETS --------------------

all: $(FILES) $(CREATE_LIBS)

include ../build-files/build-targets.mk
