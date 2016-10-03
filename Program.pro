#-------------------------------------------------
#
# Project created by QtCreator 2016-07-12T11:57:58
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Program
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    vypocetni.cpp

HEADERS  += mainwindow.h \
    vypocetni.hpp

FORMS    += mainwindow.ui


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/Qwt-6.1.3/lib/ -lqwt
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/Qwt-6.1.3/lib/ -lqwtd

INCLUDEPATH += $$PWD/Qwt-6.1.3/include
DEPENDPATH += $$PWD/Qwt-6.1.3/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/Qwt-6.1.3/lib/libqwt.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/Qwt-6.1.3/lib/libqwtd.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/Qwt-6.1.3/lib/qwt.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/Qwt-6.1.3/lib/qwtd.lib
