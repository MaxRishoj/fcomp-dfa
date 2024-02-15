TARGET=comfa
CXX=g++

SRC_DIR=src
OBJ_DIR=obj

LDFLAGS=-lstdc++
CXXFLAGS=-g -O4 -Wall

SOURCES=$(wildcard $(SRC_DIR)/*.cpp)
OBJECTS=$(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))
HEADERS=$(wildcard $(SRC_DIR)/*.h)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@ $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf $(OBJ_DIR)/*.o $(TARGET)

submit_tcp: all
	bsub < tools/hpc_submit_tcp.sh

submit_stat_k: all
	bsub < tools/hpc_submit_stat-k.sh

submit_test: all
	bsub < tools/hpc_submit_test.sh

extract_rules:
	python3 tools/rule_extract.py        data/raw_rules/snort.rules    data/full_rules_by_type snort
	python3 tools/rule_extract.py        data/raw_rules/suricata.rules data/full_rules_by_type suricata
	python3 tools/rule_extract.py --zeek data/raw_rules/zeek.rules     data/full_rules_by_type zeek
