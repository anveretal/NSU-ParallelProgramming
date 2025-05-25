cCC = gcc-14
pCC = mpicxx

T1_DIR = Task1
T2_DIR = Task2


CT1_SRC   = $(T1_DIR)/consistent_task_1.cpp
PT1V1_SRC = $(T1_DIR)/parallel_task_1_version_1.cpp
PT1V2_SRC = $(T1_DIR)/parallel_task_1_version_2.cpp

CT2_SRC   = $(T2_DIR)/consistent_task_2.cpp
PT2V1_SRC = $(T2_DIR)/parallel_task_2_version_1.cpp
PT2V2_SRC = $(T2_DIR)/parallel_task_2_version_2.cpp


CT1_BIN   = ct1
PT1V1_BIN = pt1v1
PT1V2_BIN = pt1v2

CT2_BIN   = ct2
PT2V1_BIN = pt2v1
PT2V2_BIN = pt2v2


ct1: $(CT1_SRC) | $(T1_DIR)
	$(cCC) $(CT1_SRC) -o $(CT1_BIN)

pt1v1: $(PT1V1_SRC) | $(T1_DIR)
	$(pCC) $(PT1V1_SRC) -o $(PT1V1_BIN)

pt1v2: $(PT1V2_SRC) | $(T1_DIR)
	$(pCC) $(PT1V2_SRC) -o $(PT1V2_BIN)

ct2: $(CT2_SRC) | $(T2_DIR)
	$(cCC) $(CT2_SRC) -o $(CT2_BIN)

pt2v1: $(PT2V1_SRC) | $(T2_DIR)
	$(pCC) $(PT2V1_SRC) -o $(PT2V1_BIN)

pt2v2: $(PT2V2_SRC) | $(T2_DIR)
	$(pCC) $(PT2V2_SRC) -o $(PT2V2_BIN)