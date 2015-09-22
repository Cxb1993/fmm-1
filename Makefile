CXXFLAGS =	-std=c++0x -O0 -g -pg -Wall -fmessage-length=0

SRCS=	matrix.cpp key.cpp point.cpp grid.cpp testfmm.cpp

OBJS =	$(subst .cpp,.o,$(SRCS))

LIBS =

TARGET =	testfmm

$(TARGET):	$(OBJS)
	$(CXX) -pg -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
