#FC=ifort
FC=gfortran
#FLAGS_IFORT= -mkl -CB -traceback -g 
FLAGS_IFORT= -mkl -O2
FLAGS_GCC= -O2 -llapack -lblas 
# コンパイルのために順番が大事。下層ほど先に書く。 
SRCS=\
     global_parameters.f90 \
     subroutines.f90
OBJS=$(SRCS:.f90=.o)
MAIN_SRCS=1d_chain.f90
MAIN_OBJ=1d_chain.o
PROG=1d_chain.exe

#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(MAIN_OBJ)
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -o $@ $(OBJS) $(MAIN_OBJ) $(LIBS)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(MAIN_OBJ) 
endif

# moduleをコンパイルするときの依存性を解消
%.o: %.f90
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -c $<
else
	$(FC) $(FLAGS_IFORT) -c $<
endif
%.mod: %.f90 %.o
	@true

# moduleの依存性
subroutines.o: \
  global_parameters.o 

.PHONY: clean
clean:
	mv $(PROG) $(PROG).bak; rm -f *.o *.mod core 
