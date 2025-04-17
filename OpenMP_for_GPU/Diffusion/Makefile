subdirs=ompgpu ompgpu_unified ompgpu_nowait_overlap

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)

