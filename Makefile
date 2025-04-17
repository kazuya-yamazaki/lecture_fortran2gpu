subdirs=DoConcurrent OpenACC OpenMP_for_CPU OpenMP_for_GPU

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)

