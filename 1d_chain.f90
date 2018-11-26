!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN
program main
use global_parameters
use subroutines
implicit none


character(128) OUT_FMT,command
double precision :: Ham0, Ham1
!!!!!!!!!!!!!!!!
!! variables 
! time
double precision :: time
! X(t), P(t)
double precision, allocatable :: pos1(:) ! pos_X1(1:NumBall)
double precision, allocatable :: mom1(:) ! mom_X1(1:NumBall)
double precision, allocatable :: force1(:) ! force1(1:NumBall)
! X(t+T_diff), P(t+T_diff)
double precision, allocatable :: pos2(:) ! pos_X2(1:NumBall)
double precision, allocatable :: mom2(:) ! mom_X2(1:NumBall)
double precision, allocatable :: force2(:) ! force2(1:NumBall)
! 
double precision, allocatable :: int_pos1(:) ! \int pos1 dt
double precision, allocatable :: int_mom1(:) ! \int mom1 dt
double precision, allocatable :: int_force1(:) ! \int force1 dt
double precision, allocatable :: int_pos2(:) ! \int pos2 dt
double precision, allocatable :: int_mom2(:) ! \int mom2 dt
double precision, allocatable :: int_force2(:) ! \int force2 dt
double precision, allocatable :: int_XX(:,:) ! \int X_i(s)X_j(s+t) ds
double precision, allocatable :: int_PP(:,:) ! \int P_i(s)P_j(s+t) ds
double precision, allocatable :: int_FF(:,:) ! \int F_i(s)F_j(s+t) ds
!!!!!!!!!!!!!!!!
integer :: counter,k

! パラメータの読み込み
call set_parameters

! 動的変数の割り当て
!!
allocate( pos1(1:NumBall) )
allocate( mom1(1:NumBall) )
allocate( force1(1:NumBall) )
allocate( int_pos1(1:NumBall) )
allocate( int_mom1(1:NumBall) )
allocate( int_force1(1:NumBall) )
!!
allocate( pos2(1:NumBall) )
allocate( mom2(1:NumBall) )
allocate( force2(1:NumBall) )
allocate( int_pos2(1:NumBall) )
allocate( int_mom2(1:NumBall) )
allocate( int_force2(1:NumBall) )
!!
allocate( int_XX(1:NumBall,1:NumBall) )
allocate( int_PP(1:NumBall,1:NumBall) )
allocate( int_FF(1:NumBall,1:NumBall) )
!!
allocate( k_spring2(1:NumBall) )

! 初期条件を設定
call set_initial_values(pos1,mom1,force1,pos2,mom2,force2,time)
!write(*,*) NumBall
! 

OUT_FMT=trim('(a,I3,2X,a,f6.2,2X,a,E12.2,2X,a,E12.2,2X,a,E12.2,2X,a,E12.2)')
if( write_output == 0 ) then
  open(unit=MXX_FILE, file=MXX_FILE_NAME, status='replace', action='write')
  open(unit=MPP_FILE, file=MPP_FILE_NAME, status='replace', action='write')
  open(unit=MFF_FILE, file=MFF_FILE_NAME, status='replace', action='write')
  write(MXX_FILE,OUT_FMT) "# NBall=",NumBall,"M=",Mass,"r_NL=",r_NL,"totalT=",T_tot, "T_av=",T_av, "T_diff=",T_diff
  write(MPP_FILE,OUT_FMT) "# NBall=",NumBall,"M=",Mass,"r_NL=",r_NL,"totalT=",T_tot, "T_av=",T_av, "T_diff=",T_diff
  write(MFF_FILE,OUT_FMT) "# NBall=",NumBall,"M=",Mass,"r_NL=",r_NL,"totalT=",T_tot, "T_av=",T_av, "T_diff=",T_diff
endif


int_pos1=0d0
int_pos2=0d0
int_mom1=0d0
int_mom2=0d0
int_force1=0d0
int_force2=0d0
int_XX=0d0
int_PP=0d0
int_FF=0d0
counter=0

call calc_hamiltonian2(Ham0,pos1,mom1)
do k=0,N_tot-1
  !!! integration 
  if( write_output == 0 ) then 
    call integration(int_pos1, int_pos2, int_XX, pos1, pos2, counter)
    call integration(int_mom1, int_mom2, int_PP, mom1, mom2, counter)
    call integration(int_force1, int_force2, int_FF, force1, force2, counter)
  endif
  !!!! time evolution 
  time=time+deltaT 
  call LeapFrog(pos1,mom1,force1)
  call LeapFrog(pos2,mom2,force2)

  counter = counter + 1
  if( counter == N_av ) then
    if( write_output == 0 ) then 
      call write_SV(int_pos1, int_pos2, int_XX, time, MXX_FILE)
      call write_SV(int_mom1, int_mom2, int_PP, time, MPP_FILE)
      call write_SV(int_force1, int_force2, int_FF, time, MFF_FILE)
    endif
    !! reset integration data
    int_pos1=0d0
    int_pos2=0d0
    int_mom1=0d0
    int_mom2=0d0
    int_force1=0d0
    int_force2=0d0
    int_XX=0d0
    int_PP=0d0
    int_FF=0d0
    counter = 0
  endif
  !!!!!!!!!!!!!!!!!!!!
enddo
call calc_hamiltonian2(Ham1,pos1,mom1)

write(*,*) "# H(final)-H(initial)=",Ham1-Ham0

if( write_output == 0 ) then
  close(MXX_FILE)
  close(MPP_FILE)
  close(MFF_FILE)
endif

open(unit=Outconf_FILE, file=trim("CONFIG/" // Outconf_FILE_NAME), status='replace', form='unformatted')
write(Outconf_FILE) job_number
write(Outconf_FILE) time 
write(Outconf_FILE) pos1
write(Outconf_FILE) mom1
close(Outconf_FILE)


command="cd CONFIG; ln -sf " // trim(Outconf_FILE_NAME) //  " " // trim(LASTCONF_FILE_NAME)
call system(command)

end program main



