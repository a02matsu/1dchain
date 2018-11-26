module subroutines
use global_parameters
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_parameters
use global_parameters
implicit none


open(INPUT_FILE, file=INPUT_FILE_NAME, status='old', action='READ')
  read(INPUT_FILE,*) job_number
  read(INPUT_FILE,*) NumBall
  read(INPUT_FILE,*) MASS
  read(INPUT_FILE,*) k_spring
  read(INPUT_FILE,*) r_NL
  read(INPUT_FILE,*) new_config
  read(INPUT_FILE,*) write_output
  read(INPUT_FILE,*) T_tot
  read(INPUT_FILE,*) T_av
  read(INPUT_FILE,*) T_diff
  read(INPUT_FILE,*) deltaT
close(INPUT_FILE)

k_nonlinear = k_spring * r_NL
N_tot = nint(T_tot/deltaT)
N_av = nint(T_av / deltaT)
matrix_size = NumBall

!!!
!allocate( pos1(1:NumBall) )
!allocate( mom1(1:NumBall) )
!allocate( force1(1:NumBall) )
!allocate( int_pos1(1:NumBall) )
!allocate( int_mom1(1:NumBall) )
!allocate( int_force1(1:NumBall) )
!!!
!allocate( pos2(1:NumBall) )
!allocate( mom2(1:NumBall) )
!allocate( force2(1:NumBall) )
!allocate( int_pos2(1:NumBall) )
!allocate( int_mom2(1:NumBall) )
!allocate( int_force2(1:NumBall) )
!!!
!allocate( int_XX(1:NumBall,1:NumBall) )
!allocate( int_PP(1:NumBall,1:NumBall) )
!allocate( int_FF(1:NumBall,1:NumBall) )


end subroutine set_parameters


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initial values
subroutine set_initial_values(pos1,mom1,force1,pos2,mom2,force2,time)
use global_parameters
implicit none

double precision, intent(out) :: pos1(1:NumBall) ! pos_X1(1:NumBall)
double precision, intent(out) :: mom1(1:NumBall) ! mom_X1(1:NumBall)
double precision, intent(out) :: force1(1:NumBall) ! force1(1:NumBall)
double precision, intent(out) :: pos2(1:NumBall) ! pos_X1(1:NumBall)
double precision, intent(out) :: mom2(1:NumBall) ! mom_X1(1:NumBall)
double precision, intent(out) :: force2(1:NumBall) ! force1(1:NumBall)
double precision, intent(out) :: time


integer :: n, i,j


! ディレクトリ生成
!call system('./make_directory.sh')
call system('FILE="OUTPUT"; if [ ! -d $FILE ]; then mkdir -p $FILE; fi')
call system('FILE="CONFIG"; if [ ! -d $FILE ]; then mkdir -p $FILE; fi')
!call system('FILE="SV"; if [ ! -d $FILE ]; then mkdir -p $FILE; fi')

! LASTCONF_FILE_NAME（シンボリックリンク）の設定
write(LASTCONF_FILE_NAME,'("lastconfig_N",i3.3,"NL",f5.3)') NumBall,r_NL

!! set initial Xmat1, Vmat1 and Fmat1
if( new_config == 1 ) then 
  job_number=0
  time=0d0

  !! position is zero
  do i=1,Numball
    pos1=dble(i)
  enddo
  mom1=0d0
  !! momentum is random
  !call set_random_seed
  !call BoxMuller(mom1)
  !call make_total_momentum_zero(mom1)
  call set_random_seed
  call BoxMuller(k_spring2)
  k_spring2 = dabs(k_spring2)


else
  if( job_number== 0 ) then 
    INCONF_FILE_NAME=trim("CONFIG/" // LASTCONF_FILE_NAME)
  else
    write(Inconf_FILE_NAME,'("CONFIG/config_N",i3.3,"NL",f5.3,"_",i4.4)') NumBall,r_NL,job_number-1
  endif

  open(unit=Inconf_FILE, file=Inconf_FILE_NAME, status='old', action='read',form='unformatted')

  read(Inconf_File) job_number
  job_number=job_number+1

  read(Inconf_File) time
  read(Inconf_File) pos1
  read(Inconf_File) mom1

  close(Inconf_FILE)
endif

call calc_force2(force1,pos1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! set Xmat2, Vmat2, Fmat2
pos2=pos1
mom2=mom1
force2=force1
call time_evolution(pos2,mom2,force2,T_diff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Names
write(OUTCONF_FILE_NAME,'("config_N",i3.3,"NL",f5.3,"_",i4.4)') NumBall,r_NL,job_number
write(MXX_FILE_NAME,'("OUTPUT/SVXX_N",i3.3,"NL",f5.3,"Tdiff",E9.3,"Tav",E9.3,"_",i4.4)') NumBall,r_NL,(T_diff),(T_av),job_number
write(MPP_FILE_NAME,'("OUTPUT/SVVV_N",i3.3,"NL",f5.3,"Tdiff",E9.3,"Tav",E9.3,"_",i4.4)') NumBall,r_NL,(T_diff),(T_av),job_number
write(MFF_FILE_NAME,'("OUTPUT/SVFF_N",i3.3,"NL",f5.3,"Tdiff",E9.3,"Tav",E9.3,"_",i4.4)') NumBall,r_NL,(T_diff),(T_av),job_number

end subroutine set_initial_values


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! make the total momentum to be 0
subroutine make_total_momentum_zero(mom)
implicit none

double precision, intent(inout) :: mom(:)
integer i,N
double precision :: total_mom

N=size(mom,1)
total_mom=0d0

do i=1,N
  total_mom=total_mom+mom(i)
enddo
total_mom = total_mom / dble(N)
do i=1,N
  mom(i)=mom(i) - total_mom
enddo
end subroutine make_total_momentum_zero


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Hamiltonian 
subroutine calc_hamiltonian(Ham,pos,mom)
use global_parameters
implicit none

double precision, intent(out) :: Ham
double precision, intent(in) :: pos(1:NumBall), mom(1:NumBall)
integer :: i,next

Ham=0d0
do i=1,NumBall
  if( i==NumBall ) then 
    next=1
  else
    next=i+1
  endif

  Ham = Ham &
    + mom(i)**2 / (2d0*MASS) &
    + 0.5d0*k_spring * (pos(next)-pos(i))**2 &
    + k_nonlinear * ( pos(next)-pos(i) )**4
enddo

end subroutine calc_hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate forces
subroutine calc_force(force,pos)
use global_parameters
implicit none

double precision, intent(out) :: force(1:NumBall)
double precision, intent(in) :: pos(1:NumBall)
integer i,prev,next

do i=1,NumBall
  prev=i-1
  next=i+1
  if( prev==0 ) prev=NumBall
  if( next==NumBall+1 ) next=1

  force(i) = &
      k_spring * (pos(next)-pos(i)) &
    - k_spring * (pos(i)-pos(prev)) &
    + 4d0*k_nonlinear * (pos(next)-pos(i))**3 &
    - 4d0*k_nonlinear * (pos(i)-pos(prev))**3 
enddo

end subroutine calc_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Hamiltonian 
subroutine calc_hamiltonian2(Ham,pos,mom)
use global_parameters
implicit none

double precision, intent(out) :: Ham
double precision, intent(in) :: pos(1:NumBall), mom(1:NumBall)
integer :: i,next

Ham=0d0
do i=1,NumBall
  if( i==NumBall ) then 
    next=1
  else
    next=i+1
  endif

  Ham = Ham &
    + mom(i)**2 / (2d0*MASS) &
    + 0.5d0*k_spring2(i) * pos(i)**2 & !(pos(next)-pos(i))**2 &
    + k_nonlinear * ( pos(next)-pos(i) )**4
enddo

end subroutine calc_hamiltonian2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate forces
subroutine calc_force2(force,pos)
use global_parameters
implicit none

double precision, intent(out) :: force(1:NumBall)
double precision, intent(in) :: pos(1:NumBall)
integer i,prev,next

do i=1,NumBall
  prev=i-1
  next=i+1
  if( prev==0 ) prev=NumBall
  if( next==NumBall+1 ) next=1

  force(i) = &
    - k_spring2(i) * pos(i) &
    + 4d0*k_nonlinear * (pos(next)-pos(i))**3 &
    - 4d0*k_nonlinear * (pos(i)-pos(prev))**3 
enddo

end subroutine calc_force2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! time evolution for a given time
subroutine time_evolution(pos,mom,force,t)
use global_parameters
implicit none

double precision, intent(inout) :: pos(1:NumBall)
double precision, intent(inout) :: mom(1:NumBall)
double precision, intent(inout) :: force(1:NumBall)
double precision, intent(in) :: t
integer :: num_t
integer :: i

num_t = nint( t / deltaT)

do i = 1, num_t
  call LeapFrog(pos,mom,force)
enddo
end subroutine time_evolution


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2nd order time evolution
subroutine LeapFrog(pos,mom,force)
use global_parameters
implicit none

double precision, intent(inout) :: pos(1:NumBall)
double precision, intent(inout) :: mom(1:NumBall)
double precision, intent(inout) :: force(1:NumBall)
double precision :: force_2(1:NumBall)
integer n

pos = pos + deltaT*mom/Mass + 0.5d0*deltaT*deltaT*force/Mass

call calc_force2(force_2,pos)
mom = mom + 0.5d0*deltaT*(force+force_2)
force=force_2

end subroutine LeapFrog

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate 
!!  1. \int X1(s) ds
!!  2. \int X2(s+t) ds
!!  3. \int X1(s)X2(s+t) ds
subroutine integration(int_val1, int_val2, int_VV, val1, val2, counter)
use global_parameters
implicit none

double precision, intent(inout) :: int_val1(1:NumBall)
double precision, intent(inout) :: int_val2(1:NumBall)
double precision, intent(inout) :: int_VV(1:NumBall,1:NumBall)
double precision, intent(in) :: val1(1:NumBall)
double precision, intent(in) :: val2(1:NumBall)
integer, intent(in) :: counter
double precision :: rate 
integer :: i,j

if( counter == 0 .or. counter == N_av-1 ) then
  rate=0.5d0
else
  rate=1d0
endif

int_val1=int_val1 + deltaT*rate*val1
int_val2=int_val2 + deltaT*rate*val2

do j=1,NumBall
  do i=1,NumBall
    int_VV(i,j)=int_VV(i,j) + deltaT*rate*val1(i)*val2(j)
  enddo
enddo

if( counter == N_av-1 ) then
  int_val1 = int_val1/T_av
  int_val2 = int_val2/T_av
  int_VV = int_VV/T_av
endif
end subroutine integration




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write singular value of M^{XX} to file
subroutine write_SV(vec1, vec2, mat, time, FILE_NUM)
use global_parameters
implicit none

double precision, intent(in) :: vec1(1:NumBall)
double precision, intent(in) :: vec2(1:NumBall)
double precision, intent(in) :: MAT(1:NumBall,1:NumBall)
double precision, intent(in) :: time
integer, intent(in) :: FILE_NUM

double precision :: MXX(1:NumBall,1:NumBall)
double precision :: SV(1:NumBall)
integer j,i

character(50) FMT1, FMT2

FMT1="(" // trim(FMT_time) // ",2X)"
FMT2="(" // trim(FMT_vals) // ",2X)"

write(FILE_NUM,FMT1,advance='no') time 

do j=1,NumBall
  do i=1,NumBall
    MXX(i,j) = Mat(i,j) - vec1(i)*vec2(j)
    !write(*,*) i,j,MXX(i,j)
  enddo
enddo
!write(*,*) "#####"


call Matrix_Singularvalues(SV,MXX)

do i=1,NumBall
  write(FILE_NUM,FMT2,advance='no') dlog(SV(i))
enddo
write(FILE_NUM,*)

end subroutine write_SV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Matrix_Singularvalues
subroutine Matrix_Singularvalues(sval,MAT)
implicit none

double precision, intent(in) :: MAT(:,:)
double precision, intent(out) :: sval(:) !


!lapack
character JOBU, JOBVT
integer M,N,LDA,LDU,LDVT,LWORK,INFO
integer, allocatable :: IWORK(:)
!double precision :: VL,VU
double precision, allocatable :: RWORK(:)
double precision, allocatable :: A(:,:),VT(:,:),WORK(:),U(:,:)


M=size(MAT,1)
N=size(MAT,2)

lda = m
ldu = m
ldvt = n
LWORK= 5*N
Allocate (a(m,n), u(ldu,m), VT(ldvt,n), WORK(LWORK), rwork(5*n))
A=MAT ! This is destoroyed


JOBU='N'
JOBVT='N' 

Call dgesvd(JOBU, JOBVT, m, n, a, lda, sval, u, ldu, vt, ldvt, WORK, lwork, &
        rwork, info)


end subroutine Matrix_Singularvalues


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gaussian random number
!!  BoxMuller(gauss)
!!
!! Here gauss has the structure gauss(:)
!! generage ensemble exp(-1/2(x^2) and exp(-1/2(y^2))
!! 
!! output N gaussian randum numbers
subroutine BoxMuller(gauss)
implicit none

double precision, intent(out) :: gauss(:)
double precision, parameter :: PI=dacos(-1d0)
double precision, allocatable :: rand(:)

integer :: N
integer i,a

N=size(gauss,1)

if( mod(N,2) == 0 ) then
  allocate( rand(1:N) )
else
  allocate( rand(1:N+1) )
endif

call random_number( rand )

do i=1,N/2
  gauss(2*i-1) = dsqrt(-2d0*dlog(rand(2*i-1)))*dsin(2d0*Pi*rand(2*i))
  gauss(2*i) = dsqrt(-2d0*dlog(rand(2*i-1)))*dcos(2d0*Pi*rand(2*i))
enddo
if( mod(N,2) /= 0 ) then
  gauss(N) = dsqrt(-2d0*dlog(rand(N)))*dcos(2d0*Pi*rand(N+1))
endif

end subroutine BoxMuller

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_random_seed
implicit none

integer seedsize
integer, allocatable :: seed(:)
integer i

call random_seed(size=seedsize) ! seedsizeを取得（コンパイラに依存）
allocate(seed(seedsize))
do i=1,seedsize
  call system_clock(count=seed(i))
enddo

call random_seed(put=seed(:)) ! seedの初期値を与える
end subroutine set_random_seed


end module subroutines
