CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
    maslovaes.cpp \
    main.cpp \
    zhalninrv.cpp \
    isokovaa.cpp \
    golovatyukam.cpp \
	  kirdyushkindv.cpp \
	  puzinva.cpp \
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
