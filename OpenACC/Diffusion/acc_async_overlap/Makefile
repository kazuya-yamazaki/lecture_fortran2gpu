F90  = mpif90
RM  = rm -f

FFLAGS    = -fast -acc -Minfo=accel -gpu=cc90,mem:separate -cpp
LDFLAGS   = 

SRCS   = misc.f90 diffusion.f90 main.f90
TARGET = ../run/acc_async_overlap
DISTTARGET = $(TARGET)_1.0.0

OBJS += $(filter %.o,$(SRCS:%.f90=%.o))


.PHONY: all
all : $(TARGET)

$(TARGET) : $(OBJS)
	$(F90) $(FFLAGS) $(TARGET_ARCH) $(OBJS) -o $@ $(LDFLAGS)

%.o : %.f90
	$(F90) $(FFLAGS) $(TARGET_ARCH) -c $<


.PHONY: clean
clean :
	$(RM) $(TARGET)
	$(RM) $(OBJS)
	$(RM) *.mod
	$(RM) *~


