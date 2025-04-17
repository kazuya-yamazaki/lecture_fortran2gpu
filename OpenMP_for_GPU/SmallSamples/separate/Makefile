subdirs=01_single_loop 02_nested_loop 03_nested_loop_collapse 04_reduction_loop 05_reduction_nested_loop 06_atomic_loop 07_workshare_loop 08_function_call

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)

