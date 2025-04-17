OMPVER=5 #4 or 5

F90  = nvfortran
RM  = rm -f

FFLAGS    = -O3 -mp=gpu -gpu=cc80 -Minfo=mp -Mpreprocess -DOMP=${OMPVER}
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


