!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! global parameters
module global_parameters
implicit none

integer :: NumBall
double precision :: MASS
double precision :: k_spring
double precision :: k_nonlinear
double precision, allocatable :: k_spring2(:)
double precision :: r_NL


!! parameters to control the simulation
integer :: job_number 
integer :: new_config
!double precision :: time ! 現在時刻
double precision :: deltaT ! 刻み幅
double precision :: T_tot ! 全シミュレーション時間
double precision :: T_av ! この時間だけ積分してoutputする。
double precision :: T_diff ! ノートで言うところの"t"
integer :: write_output

integer :: N_tot  ! nint(totalT / deltaT)
integer :: N_av ! nint(T_av / deltaT)
integer :: matrix_size


!! variables 
! X(t), P(t)
!double precision, allocatable :: pos1(:) ! pos_X1(1:NumBall)
!double precision, allocatable :: mom1(:) ! mom_X1(1:NumBall)
!double precision, allocatable :: force1(:) ! force1(1:NumBall)
!! X(t+T_diff), P(t+T_diff)
!double precision, allocatable :: pos2(:) ! pos_X2(1:NumBall)
!double precision, allocatable :: mom2(:) ! mom_X2(1:NumBall)
!double precision, allocatable :: force2(:) ! force2(1:NumBall)
!! 
!double precision, allocatable :: int_pos1(:) ! \int pos1 dt
!double precision, allocatable :: int_mom1(:) ! \int mom1 dt
!double precision, allocatable :: int_force1(:) ! \int force1 dt
!double precision, allocatable :: int_pos2(:) ! \int pos2 dt
!double precision, allocatable :: int_mom2(:) ! \int mom2 dt
!double precision, allocatable :: int_force2(:) ! \int force2 dt
!double precision, allocatable :: int_XX(:,:) ! \int X_i(s)X_j(s+t) ds
!double precision, allocatable :: int_PP(:,:) ! \int P_i(s)P_j(s+t) ds
!double precision, allocatable :: int_FF(:,:) ! \int F_i(s)F_j(s+t) ds

character(10), parameter :: FMT_time="E15.8"
character(10), parameter :: FMT_vals="E23.15"


character(128) :: LASTCONF_FILE_NAME
character(128) :: INPUT_FILE_NAME="input.dat"
character(128) :: Inconf_FILE_NAME
character(128) :: Outconf_FILE_NAME
!!!
character(128) :: MXX_FILE_NAME
character(128) :: MPP_FILE_NAME
character(128) :: MFF_FILE_NAME

!integer, parameter :: INPUT_FILE=10
integer, parameter :: INPUT_FILE=10
integer, parameter :: Inconf_FILE=11
integer, parameter :: Outconf_FILE=12
!!!
integer, parameter :: MXX_FILE=20
integer, parameter :: MPP_FILE=21
integer, parameter :: MFF_FILE=22



end module global_parameters


