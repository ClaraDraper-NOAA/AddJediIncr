SHELL=	/bin/ksh
FC      =$(FCMP)
LIBS    =$(LIBSM)
OBJS	=apply_incr_noahmp_snow.o

CMD     =apply_incr
$(CMD):	$(OBJS)
	${FC} ${DEBUG} $(FFLAGS) $(LDFLG) $(OBJS) $(LIBS) -o $(CMD)

apply_incr_noahmp_snow.o: apply_incr_noahmp_snow.f90
	$(FC)  $(FFLAGS) $(OMPFLAG) -c apply_incr_noahmp_snow.f90

install:
	-cp -p $(CMD) ../exec/.
clean:
	-rm -f $(OBJS) *.mod $(CMD)
