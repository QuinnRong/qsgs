
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

run: $(TARGET)
	sbatch qsgs_3d.slurm

run_mac: $(TARGET)
	./$(TARGET) run_1_valid 2 20
	./$(TARGET) run_2_train 8 20
	./$(TARGET) run_3_test 2 20

.PHONY:clean
clean:
	-rm -rf $(TARGET) bin $(OUTDIR) *.err *.out
