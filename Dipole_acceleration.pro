TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ground_state.cpp \
        main.cpp

HEADERS += \
    constants.h \
    ground_state.h

LIBS += /usr/local/lib/libfftw3.a
