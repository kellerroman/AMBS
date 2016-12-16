program calc_error
implicit none
integer, parameter :: n_sol = 19
integer :: fu,i,ios,n_sol_sim
integer :: pos_in_ref
real(kind=8) :: sol(n_sol),pos(n_sol)
real(kind=8) :: xpos, xval
real(kind=8) :: error
real(kind=8), allocatable :: sol_sim(:),pos_sim(:)
write(*,*) "calculating residual error for sod shock tube solution"
open(newunit=fu,file="str.dx")
do i = 1, n_sol
   read(fu,*) pos(i), sol(i)
end do 
close(fu)
open(newunit=fu,file="data_1d.csv")
read(fu,*)
n_sol_sim = 0
do
   read(fu,*,iostat=ios)
   if (ios /= 0) then
      exit
   end if
   n_sol_sim = n_sol_sim + 1

end do
n_sol_sim = n_sol_sim - 1
!write(*,*) n_sol_sim
allocate(sol_sim(n_sol_sim))
allocate(pos_sim(n_sol_sim))
rewind(fu)
read(fu,*)
do i = 1, n_sol_sim
   read(fu,*) pos_sim(i), sol_sim(i)
end do
close(fu)
pos_in_ref = 1
error = 0.0D0
open(newunit=fu,file="error.csv")
do i = 1, n_sol_sim
   if (pos_sim(i) >= pos(pos_in_ref+1)) then
      pos_in_ref = pos_in_ref + 1
      if (pos(pos_in_ref+1) == pos(pos_in_ref)) then
         pos_in_ref = pos_in_ref + 1
      end if
      !write(*,*) pos_sim(i), pos(pos_in_ref)
   end if
   xpos = (pos(pos_in_ref+1) - pos_sim(i)) / (pos(pos_in_ref+1) - pos(pos_in_ref))
   xval = sol(pos_in_ref) * xpos + sol(pos_in_ref+1) * (1.0d0-xpos)
   write(fu,*) pos_sim(i), sol_sim(i), xval, sol_sim(i) - xval
   error = error + (sol_sim(i) - xval) * (sol_sim(i) - xval)
end do
error = sqrt(error / n_sol_sim)
write(*,*) error
close(fu)
end program
