CC = mpic++ 
IFLAGS = -I ../include/
CFLAGS = -lpthread -lrt -lm -fopenmp
TARGETS = test_pagerank test_sssp test_cc test_bfs test_bc
.PHONY: all
all: $(TARGETS)	

test_pagerank: test_pagerank.cpp
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $^
test_sssp: test_sssp.cpp
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $^
test_cc: test_cc.cpp
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $^
test_bfs: test_bfs.cpp
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $^
test_bc: test_bc.cpp
	$(CC) $(IFLAGS) $(CFLAGS) -o $@ $^


clean:
	rm -rf $(TARGETS)
