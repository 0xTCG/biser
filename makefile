all: sedef unionfind biser 

sedef:
	cd src/sedef && make -j release

unionfind:
	g++ -O3 -w -o exec/uf src/uf.cpp
	

SEQ := $(shell command -v seqc)
CLANG := $(shell command -v clang)

SEQ_PATH := $(shell dirname $(SEQ))

export LD_LIBRARY_PATH := $(shell cd $(SEQ_PATH) && cd .. && cd lib/seq && pwd)
export SEQ_LIBRARY_PATH := $(LD_LIBRARY_PATH)

biser:
ifndef SEQ
	$(error "No Seq found. Please install Seq to proceed.")
endif
ifndef CLANG
	$(error "No CLANG found. Please install CLANG to proceed.")
endif
	seqc build -release src/biser_search.seq -o exec/biser_search
	seqc build -release src/big_cluster.seq -o exec/big_cluster
	seqc build -release src/decompose.seq -o exec/decompose


