CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
    main.cpp \
    zhalninrv.cpp \
	  parshinad.cpp \
	  malovki.cpp \
    landyshevav.cpp \
    garinma.cpp \
    simatovvv.cpp \
    guskovas.cpp \
    kozinasa.cpp \
    sayfetdinovsf.cpp \
    borisovayu.cpp \
    lab.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vvm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o vvm
