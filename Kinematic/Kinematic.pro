TEMPLATE	= app
LANGUAGE	= C++

CONFIG	+= qt warn_on debug 
DEFINES += BUILD_DOUBLE USE_OSG

SOURCES	+= main.cpp ObjectFormImpl.cpp KinematicFormImpl.cpp 
HEADERS += ObjectFormImpl.h KinematicFormImpl.h KinematicForm.h ObjectForm.h
INCLUDEPATH += ../include /usr/include/libxml2 /usr/include/qt3 
LIBS += -L.. -losg -losgGA -losgDB -losgViewer -losgManipulator -losgUtil -lOpenThreads -lMoby -lqwt 

unix {
  MOC_DIR = .moc
  OBJECTS_DIR = .obj
}


QT +=  qt3support 

