CC         = g++
LD         = g++ 
CFLAGS      = -Wall -g -std=c++11
PROG_NAME  = sconce
BOOSTFLAGS = -lboost_system -lboost_filesystem -lboost_program_options
GSLFLAGS   = `gsl-config --cflags --libs` 

SRC_DIR    = ./src
BUILD_DIR  = ./build
SRC_LIST   = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_LIST   = $(addprefix $(BUILD_DIR)/, $(notdir $(SRC_LIST:.cpp=.o)))

.PHONY: all clean

all: $(PROG_NAME)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(PROG_NAME): $(OBJ_LIST)
	$(LD) $(GSLFLAGS) $(BOOSTFLAGS) $(OBJ_LIST) -o $@

clean:
	rm -f $(PROG_NAME) $(BUILD_DIR)/*.o

