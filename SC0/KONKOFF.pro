QT -= gui

CONFIG += c++17
QMAKE_CXXFLAGS += -std=c++11

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    bcr.cpp \
    cell.cpp \
    chemokines3D.cpp \
    GC3D.cpp \
    lattice.cpp \
    mafalda.cpp \
    main.cpp \
    output.cpp \
    parameters.cpp \
    random.cpp \
    setparam.cpp \
    trackball.cpp \
    vector3d.cpp

HEADERS += \
	bcr.h \
	cell.h \
	chemokines3d.h \
	dynarray.h \
	GC3D.h \
	lattice.h \
	mafalda.h \
	main.h \
	output.h \
	parameters.h \
	performance.h \
	random.h \
	setparam.h \
	trackball.h \
	vector3d.h
