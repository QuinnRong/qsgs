
CXX = g++

CXXFLAGS = -std=c++11 -Iinclude

TARGET = qsgs_3d

SRCDIR = ./src/
SRCS = $(wildcard $(SRCDIR)*.cpp)

OBJDIR = ./bin/
OBJS = $(addprefix $(OBJDIR), $(notdir $(patsubst %.cpp, %.o, $(SRCS))))

ifeq ($(wildcard $(OBJDIR)), )
$(shell mkdir -p $(OBJDIR))
endif

$(TARGET):$(OBJS)
	$(CXX) $^ -o $@

$(OBJDIR)%.o:$(SRCDIR)%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

OUTDIR = ./output/
ifeq ($(wildcard $(OUTDIR)), )
$(shell mkdir -p $(OUTDIR))
endif

run_windows:
	./qsgs_3d.exe

run_hpc:
	./qsgs_3d

run: $(TARGET)
	sbatch qsgs_3d.slurm

.PHONY:clean
clean:
	-rm -rf $(TARGET) bin $(OUTDIR)* *.err *.out
