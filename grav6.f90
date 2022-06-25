program gravimetry
implicit none

integer :: h,i,j,k
integer, parameter :: ndim=3, n=10, nprobe=1, nsteps=n**ndim
real, parameter :: dx=0.1, L=n*dx, dt=0.01, pi=4.D0*DATAN(1.D0)
real*4, dimension(n,n,n) :: rho_in !match with "ndim"
real*4 :: X(nprobe,nsteps,ndim), V(nprobe,nsteps,ndim), a(nprobe,nsteps,ndim), inst_acc(ndim), v_tang(ndim)
real :: ratio, rad1,rad2, radius, theta, phi, raggio,psi1,psi2

call gen_rho(ndim,n,dx,L,rho_in)
print *, 'the generated rho, sliced at n/2, is'
call print_matrix(n,n,rho_in(:,n/2,:)) !CLARIFY WORKING OF SLICING!! VERY COMFY!

!launch probes

            do h=1,nprobe
call random_number(ratio)
call random_number(theta)
theta=theta*2*pi
call random_number(phi)
phi=phi*pi
print *, 'ratio,theta,phi=', ratio, theta, phi

    rad1= L*(3**0.5)
    rad2= L*(3**0.5)*3
    radius= rad1+(rad2-rad1)*ratio
    X(h,1,1)= radius*sin(phi)*cos(theta)
    X(h,1,2)= radius*sin(phi)*sin(theta)
    X(h,1,3)= radius*cos(phi)

    call a_tang_vect(X(h,1,:),ndim,v_tang)
do k=1,ndim
V(h,1,k)= v_tang(k)
end do  !how to introduce adimensionalization to get nice trajs?

        do i=2,nsteps
            do k=1,ndim
        inst_acc(k)=0
            end do

call inst_acc_comp(X(h,i-1,:),rho_in,ndim,n,dx,inst_acc, i)
    do k=1,ndim !update pos and speed
    X(h,i,k)=X(h,i-1,k)+V(h,i-1,k)*dt
    V(h,i,k)=V(h,i-1,k)+inst_acc(k)*dt !+ or -? 
    a(h,i-1,k)=inst_acc(k)
    end do

call compute_polars(X(h,i-1,:),ndim,psi1,psi2,raggio)

 !   if (mod(i,10)==0) then  !PRINTIIING
        print *, 'cycle_____________n', i
	        80 FORMAT ('',10F6.1)
            print *, 'X='
            write(*,80) (X(1,i-1,j),j=1,ndim )
            print *, 'raggio,psi1,psi2='
            write(*,80) raggio,psi1,psi2           
            ! print *, '----'
            ! print *, 'V='
            ! write(*,80) (V(1,i-1,j),j=1,ndim )
            ! print *, '----'
            ! print *, 'inst_acc='
            ! write(*,80) (inst_acc(j),j=1,ndim )
            ! print *, '--------------------------'
 !       end if
        end do
            end do 

!when hits the domain quit? STABILITY ISSUES!
!launch another body UNTIL: 
!number of steps done along all bodies= n^ndim

end program

!________________________________________________
SUBROUTINE gen_rho (ndim,n,dx, L,rho_in)
implicit none
!intent IN
integer,intent(in) :: ndim,n
real, intent(in) :: dx , L
!intent OUT
real*4, dimension(n,n,n) :: rho_in !match with "ndim"
!declare local
integer :: h,i,j
real*4 :: radius
!computations

do h=1,n
do i=1,n
do j=1,n
    radius=((h*dx-L/2)**2+(i*dx-L/2)**2+(j*dx-L/2)**2)**0.5
    if (radius<L/3) then
    rho_in(h,i,j)= 2- radius/L
    else
    rho_in(h,i,j)=0
    end if
end do
end do
end do

END SUBROUTINE gen_rho
!_______________________________________________

SUBROUTINE inst_acc_comp(Xp,rho_in,ndim,n,dx,acc, count)
implicit none
!intent IN
integer,intent(in) :: ndim,n, count
real,intent(in) :: dx
real*4,intent(in) :: Xp(ndim),rho_in(n,n,n)
!intent OUT
    !real*4, dimension(ndim),intent(out) :: acc
real*4,intent(out) :: acc(ndim)
!declare local
integer :: h,i,j,k
real*4 :: r2, radius(ndim)
real :: Xe(ndim)


    do k=1,3
        acc(k)=0
    end do

!computations
do h=1,n !along x
do i=1,n !along y
do j=1,n !along z
r2=0
Xe(1)=h*dx  !optimize syntax?
Xe(2)=i*dx
Xe(3)=j*dx

    do k=1,ndim
r2=r2+(Xp(k)-Xe(k))**2
radius(k)=Xe(k)-Xp(k)
    end do

    do k=1,ndim 
acc(k)=acc(k)+ 0.3*rho_in(h,i,j)*(radius(k))/(r2**(1.5)) !check vectorial gravity
    end do  

end do
end do
end do

END SUBROUTINE inst_acc_comp
!________________________________________________
SUBROUTINE a_tang_vect(X,ndim,vn_tang)
implicit none
!intent IN
integer,intent(in) :: ndim
real*4, intent(in) :: X(ndim)
!intent OUT
real, intent(out) :: vn_tang(ndim)
!locals
integer :: i,j,k
real*4 :: v_tang(ndim)
    call random_number(v_tang(1))
    call random_number(v_tang(2))
v_tang(3)= -( X(1)*v_tang(1) + X(2)*v_tang(2) )/X(3)
call vect_normalizer(v_tang,ndim,vn_tang)
END SUBROUTINE a_tang_vect
!_______________________________________________
SUBROUTINE vect_normalizer(vector,ndim,normalized)
implicit none
!intent IN
integer,intent(in) :: ndim
real*4, intent(in) :: vector(ndim)
!intent OUT
real, intent(out) :: normalized(ndim)
!locals
integer :: i,j,k
real :: norm2

norm2=0
do i=1,ndim
norm2=norm2+vector(i)**2
end do

do i=1,ndim
normalized(i)=vector(i)/norm2**0.5
end do

END SUBROUTINE vect_normalizer
!_______________________________________________
SUBROUTINE compute_polars(X,ndim,psi1,psi2,r)
implicit none
!intent IN
integer,intent(in) :: ndim
real*4, intent(in) :: X(ndim)
!intent OUT
real, intent(out) :: psi1,psi2,r
!locals
integer :: i,j,k
real*4 :: pi=4.D0*DATAN(1.D0)

r=(X(1)**2+X(2)**2+X(3)**2)**0.5
psi2=180/pi*acos(X(3)/r)
psi1= 180/pi*sign(X(2),real(1))* acos(X(1)/(r*sin(psi2)))

END SUBROUTINE compute_polars
!________________________________________________
!________________________________________________
SUBROUTINE print_matrix (nr,nc,mat)
implicit none
integer,intent(in) :: nr,nc
real, intent(in) :: mat(nr,nc)
integer :: i,j
do i=1,nr
	!print *, (mat(i,j), j=1,nc )
	write(*,80) (mat(i,j),j=1,nc )
	80 FORMAT ('',10F6.1)
end do
print *, '-------------------------------------------'
END SUBROUTINE print_matrix
!_______________________________________________