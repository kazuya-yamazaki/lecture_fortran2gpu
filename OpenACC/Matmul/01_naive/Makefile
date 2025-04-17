F90  = nvfortran
RM  = rm -f

FFLAGS    = -fast -mp -acc -gpu=cc90,mem:separate -Minfo=accel -Mpreprocess
LDFLAGS   = 

SRCS   = main.f90
#TARGET += $(filter %.exe,$(SRCS:%.f90=%.exe))
TARGET = run
DISTTARGET = $(TARGET)_1.0.0

OBJS += $(filter %.o,$(SRCS:%.f90=%.o))


.PHONY: all
all : $(TARGET)

$(TARGET) : $(OBJS)
	$(F90) $(FFLAGS) $(OBJS) -o $@ $(LDFLAGS)

# %.exe : %.o
# 	$(F90) $(FFLAGS) $< -o $@ $(LDFLAGS)

%.o : %.f90
	$(F90) $(FFLAGS) -c $<


.PHONY: clean
clean :
	$(RM) $(TARGET)
	$(RM) $(OBJS)
	$(RM) *.mod
	$(RM) *~


