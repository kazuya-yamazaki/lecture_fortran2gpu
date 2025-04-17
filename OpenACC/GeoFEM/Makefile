subdirs=ACC_S_a ACC_S_c ACC_U_a ACC_U_c

.PHONY: all $(subdirs)

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(subdirs)

