CXX := g++
CXXFLAGS := -c -O3 -fopenmp
INCLUDE := 
LDFLAGS := -fopenmp 
SRCDIR := src
OBJDIR := obj
BINDIR := bin

EXEFILE := SimIPS

cppfiles := $(wildcard $(SRCDIR)/*.cpp)
objects := $(cppfiles:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
deps := $(objects:.o=.d)

.PHONY: all clean

all: $(BINDIR)/$(EXEFILE)

-include $(deps)

$(OBJDIR)/%.d: $(SRCDIR)/%.cpp
	mkdir -p $(@D)
	$(CXX) -MM -MT "$@ $(patsubst %.d, %.o, $@)" -MF $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $<

$(BINDIR)/$(EXEFILE): $(objects)
	mkdir -p $(@D)
	$(CXX) -o $@ $(LDFLAGS) $^

clean:
	$(RM) -r -f $(OBJDIR)
	$(RM) -r -f $(BINDIR)


