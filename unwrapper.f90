!//################# TrajXCUDA #########################
!//######### Suresh Kondati Natarajan ##################
!//##### Lehrstuhl fuer theoretische chemie ############
!//######## Ruhr Universitaet Bochum ###################



program unfold
implicit none


double precision,dimension(:,:,:),allocatable :: rf,rl
double precision, dimension(3)::r,rdiff
character(len=2),dimension(:,:),allocatable ::ELT
integer :: nstruct,natoms,i,j,k
double precision,dimension(4) :: l
double precision :: ichk
CHARACTER(len=32) :: arg

  !DO i = 1, iargc()
              CALL getarg(1, arg)
!		CALL getarg(2,nstruct)           
 !END DO
nstruct=10000


!call system('grep -c Atoms input.xyz > nstruct')
!open(1,file="nstruct",status="old")
!read(1,*)nstruct
!close(1)

open(10,file="input.xyz",status="old")
read(10,*)natoms
rewind(10)
open(20,file="input.xyz_unwrapped",status="replace")
allocate(rf(nstruct,natoms,3),rl(nstruct,natoms,3),ELT(nstruct,natoms))

l(1)=0
l(2)=0
l(3)=0
l(4)=0
do i=1,nstruct
read(10,*)
read(10,*)
do j=1,natoms
read(10,*)ELT(i,j),rf(i,j,:)
enddo
enddo
if(arg .eq. 'ortho') then
read(10,*)l(1),l(2),l(3)
elseif(arg .eq. 'mono') then
read(10,*)l(1),l(2),l(3),l(4)
endif
close(10)

if(l(4) == 0) then
do i=1,nstruct
write(20,*) natoms
write(20,*) "Atoms. ",l(:)

 do j=1,natoms
  if ( i == 1 ) then
   rl(i,j,:) = rf(i,j,:)
   write(20,*)ELT(i,j),rl(i,j,:)
  else
   do k=1,3
    r(k) = rf(i,j,k) - rf(i-1,j,k)
    if ( abs(r(k)) .gt. l(k)/2 ) then
     if (r(k) .gt. 0) then
      r(k) = abs(r(k)) - l(k)
      rl(i,j,k) = rl(i-1,j,k) - abs(r(k))
     else
      r(k) = abs(r(k)) - l(k)
      rl(i,j,k) = rl(i-1,j,k) + abs(r(k))
     endif
    else
     rl(i,j,k) = rl(i-1,j,k) + r(k)
    endif
   enddo
   write(20,*)ELT(i,j),rl(i,j,:)
  endif
 enddo
enddo
write(20,*)l(1),l(2),l(3)
close(20)
!l4 !=0
else
write(*,*) "bibli"
do i=1,nstruct
write(20,*) natoms
write(20,*) "Atoms. ",l(:)

 do j=1,natoms
  if ( i == 1 ) then
   rl(i,j,:) = rf(i,j,:)
   write(20,*)ELT(i,j),rl(i,j,:)
  else
    r(1) = rf(i,j,1) - rf(i-1,j,1)
    r(2) = rf(i,j,2) - rf(i-1,j,2)
    r(3) = rf(i,j,3) - rf(i-1,j,3)
    !r(k) = rf(i,j,k) - rf(i-1,j,k)
    if ( abs(r(2)) .gt. l(2)/2 ) then
     if (r(2) .gt. 0) then
      r(2) = abs(r(2)) - l(2)
      rl(i,j,2) = rl(i-1,j,2) - abs(r(2))
      r(1) = abs(r(1)) - l(4)
      rl(i,j,1) = rl(i-1,j,1) - abs(r(1))
     else
      r(2) = abs(r(2)) - l(2)
      rl(i,j,2) = rl(i-1,j,2) + abs(r(2))
      r(1) = abs(r(1)) - l(4)
      rl(i,j,1) = rl(i-1,j,1) + abs(r(1))
     endif
    else
     rl(i,j,2) = rl(i-1,j,2) + r(2)
    endif
    if ( abs(r(1)) .gt. l(1)/2 ) then
     if (r(1) .gt. 0) then
      r(1) = abs(r(1)) - l(1)
      rl(i,j,1) = rl(i-1,j,1) - abs(r(1))
     else
      r(1) = abs(r(1)) - l(1)
      rl(i,j,1) = rl(i-1,j,1) + abs(r(1))
     endif
    else
     rl(i,j,1) = rl(i-1,j,1) + r(1)
    endif
    if ( abs(r(3)) .gt. l(3)/2 ) then
     if (r(3) .gt. 0) then
      r(3) = abs(r(3)) - l(3)
      rl(i,j,3) = rl(i-1,j,3) - abs(r(3))
     else
      r(3) = abs(r(3)) - l(3)
      rl(i,j,3) = rl(i-1,j,3) + abs(r(3))
     endif
    else
     rl(i,j,3) = rl(i-1,j,3) + r(3)
    endif    
   write(20,*)ELT(i,j),rl(i,j,:)
  endif
 enddo
enddo
write(20,*)l(:)
close(20)  
endif

end program unfold
