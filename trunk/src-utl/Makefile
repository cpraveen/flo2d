include ../Makefile.in

FFLAGS = $(CFLAGS)
LDFLAGS= -I../src-flo

TARGETS = deform_ellipse

deform_ellipse = deform_ellipse.o \
                 set_bdpts.o \
                 rbf_train.o \
                 rbf_eval.o \
                 divergence.o \
                 ../src-flo/read_input.o \
                 ../src-flo/read_grid.o

all: $(TARGETS)

deform_ellipse: $(deform_ellipse)
	$(FC) -o deform_ellipse $(deform_ellipse) -llinpack -lblas

.f.o:	; $(FC) $(FFLAGS) $(LDFLAGS) -c -o $*.o $*.f

clean:
	rm -f *.o $(TARGETS)
