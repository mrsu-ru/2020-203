CC=g++
CFLAGS=-c -Wall -g
LDFLAGS= -g
SOURCES= \
    maslovaes.cpp \
    main.cpp \
    zhalninrv.cpp \
    edelevaup.cpp \
    ashryatovarr.cpp \
    kotkovsn.cpp \
    kvashninka.cpp \
    bochkarevda.cpp \
    kazakovais.cpp \
    isokovaa.cpp \
    golovatyukam.cpp \
	gorbunovaa.cpp \
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
	kozlovaes.cpp \
	loginovvv.cpp \
    manindi.cpp \
    zevaykinae.cpp \
    dvoryaninovada.cpp \
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
