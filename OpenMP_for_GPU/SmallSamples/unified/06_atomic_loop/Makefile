OMPVER=5 #4 or 5

F90  = nvfortran
RM  = rm -f

FFLAGS    = -fast -mp=gpu -gpu=cc90,mem:unified -Minfo=mp -Mpreprocess -DOMP=${OMPVER}
# Unified Memory非対応機では、mem:unified を mem:managed に変える
# mem:unified:nomanagedallocだとSegmentation faultになる
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


