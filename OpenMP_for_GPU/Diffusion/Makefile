subdirs=ompgpu ompgpu_unified

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)

