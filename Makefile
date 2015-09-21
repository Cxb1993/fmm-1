CXXFLAGS =	-O0 -g -pg -Wall -fmessage-length=0

OBJS =		testfmm.o 

LIBS =

TARGET =	testfmm

$(TARGET):	$(OBJS)
	$(CXX) -pg -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
