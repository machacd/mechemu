!this file contains functions specific to the emulator
MODULE lin_mod

USE service_functions
IMPLICIT NONE

type :: linear_model_data

!input and design data of the emulator
real,allocatable :: input(:), output(:,:), parameters(:,:), parameters_physical(:)
real,allocatable :: F(:,:), g(:), H(:),I(:,:),cor_factor_multi(:),real_world(:),indices(:)
!hyperparameters of the emulator
real,allocatable :: k_lam(:),k_loss(:),k_delay(:)
real ::  gamma, cor_factor,delta_t,k_level,V_ini_inp,E_ini_inp
INTEGER :: t_max, n_used, n_test_sets, m, dim_obs, no_of_pars,mode,&
  lambda_dim,input_dim

!c_emu
real :: T_indices_store(3),det
real,allocatable :: T_store(:,:),hT_store(:,:),sigma_tilde(:,:,:,:,&
                :,:),sigma(:,:),z_prime(:),z(:),y(:),&
                g_store(:,:,:,:),g_store_var(:,:,:,:),estimation_vector(:)
real :: comp_time_1, comp_time_2, comp_time_total

!estimation storage variables
real,allocatable :: estimation_variance(:,:,:,:), estimation_mean(:,:), results(:)

!kalman filter storage variables 
real,allocatable :: states_means(:,:,:), states_variances(:,:,:,:),&
states_means_obs(:,:,:),states_variances_obs(:,:,:,:),e_obs_total(:,:,:,:),&
var_e_obs_total(:,:,:,:)

!initial matrices
real,allocatable :: sigma_ini(:,:),v_ini(:,:),e_ini(:)

end type linear_model_data
!linear model

contains

subroutine read_data(this)
type(linear_model_data):: this
  INTEGER :: i,j
  integer :: tm,dim,nu,nt,m
  tm=this%t_max
  dim=this%dim_obs
  nu=this%n_used
  nt=this%n_test_sets
  m=this%m

 this%n_test_sets=1
 this%delta_t=1
 this%gamma=2.0
  !the folowing is used, if we look at correlation lengths separately for each parameter
  !in the kalman filter emulator library
  ! IF (ALLOCATED(this%cor_factor_multi)) THEN
  !   DEALLOCATE(this%cor_factor_multi)
  ! ENDIF
  ! ALLOCATE(this%cor_factor_multi(this%no_of_pars))
  ! this%cor_factor_multi=1
! DO i=1,this%no_of_pars
!     READ(2068,*) this%cor_factor_multi(i)   
! END DO
! REWIND(2068)
  ! this%cor_factor_multi=this%cor_factor

  IF (ALLOCATED(this%input)) THEN
      DEALLOCATE(this%input)
  ENDIF
  allocate(this%input(tm))
  this%input=0

!also quite provblem specific
  IF (ALLOCATED(this%parameters_physical)) THEN
      DEALLOCATE(this%parameters_physical)
  ENDIF
  allocate(this%parameters_physical(5)) !just a simplification, I don't care about all the pars
  this%parameters_physical=0

!read auxiliary parameters like lambdas and losses, from the file lambda.cfg

IF (ALLOCATED(this%k_lam)) THEN
  DEALLOCATE(this%k_lam)
ENDIF
ALLOCATE(this%k_lam(this%lambda_dim))
this%k_lam=0

IF (ALLOCATED(this%k_loss)) THEN
  DEALLOCATE(this%k_loss)
ENDIF
ALLOCATE(this%k_loss(this%input_dim))
this%k_loss=0

IF (ALLOCATED(this%k_delay)) THEN
  DEALLOCATE(this%k_delay)
ENDIF
ALLOCATE(this%k_delay(this%input_dim))
this%k_delay=0


!next loop is also PROBLEM SPECIFIC
 ! this%input=this%input/(60.0/this%delta_t)*0.001 ! convert to m/(min*delta_t)
  this%input=this%input*this%parameters_physical(5)*10000 ! multiply by the area
  ! (in ha)

! create unitary matrix to save few lines...

  IF (ALLOCATED(this%I)) THEN
    DEALLOCATE(this%I)
  ENDIF
  ALLOCATE(this%I(m,m))
  this%I=0
  do i=1,m 
    this%I(i,i)=1
  end do

  this%comp_time_1= 0
  this%comp_time_2= 0

IF (ALLOCATED(this%z_prime)) THEN
    DEALLOCATE(this%z_prime)
ENDIF
allocate(this%z_prime(this%n_used*this%t_max*this%m))
this%z_prime=0

end subroutine

function getF(this, alpha) result(out)
type(linear_model_data) :: this
INTEGER :: alpha
integer :: i
real :: out(this%m,this%m),lambda_store(this%lambda_dim)
out=0
lambda_store=lambda(this,alpha)
DO i=1,this%m   
    out(i,i)=-lambda_store(i)
    !this condition is only for the purpose of writing a paper...
    IF (i>this%input_dim) THEN 
         out(i,i-1)=lambda_store(i-1)
        IF (this%input_dim==2) THEN 
         out(3,1)=lambda_store(1)
         out(3,2)=lambda_store(2)
        endif
    ENDIF
END DO
End function

subroutine allocate_initial(this)
IMPLICIT none
type(linear_model_data) :: this
integer :: i

IF (ALLOCATED(this%E_ini)) THEN
  DEALLOCATE(this%E_ini)
ENDIF
ALLOCATE(this%E_ini(this%m))
this%E_ini=this%E_ini_inp

IF (ALLOCATED(this%V_ini)) THEN
  DEALLOCATE(this%V_ini)
ENDIF
ALLOCATE(this%V_ini(this%m,this%m))
this%V_ini=0
do i=1,this%m
  this%V_ini(i,i)=this%V_ini_inp
end do

IF (ALLOCATED(this%sigma_ini)) THEN
  DEALLOCATE(this%sigma_ini)
ENDIF
ALLOCATE(this%sigma_ini(this%m,this%m))
this%sigma_ini=0.0
do i=1,this%m
  this%sigma_ini(i,i)=1
end do

end subroutine

function getFexp(this, alpha) result(out)
type(linear_model_data) :: this
INTEGER :: alpha
real :: out(this%m,this%m)
out=0
out=expo_mat(getF(this,alpha),this%delta_t,this%m)
End function

function getk(this, alpha,time,Fexp) result(out)
type(linear_model_data) :: this
INTEGER :: alpha, time,i,a
real :: out(this%m)
real :: Fexp(this%m,this%m)
out=0
! input is only on the first of the reservoirs.
do i =1,this%input_dim
  if (time-this%k_delay(i) .gt. 1) then
    out(i)=this%input(time-int(this%k_delay(i)))*this%k_loss(i)*this%parameters(alpha,1)
  else
    out(i)=0
  end if
end do
out=-matmul(matmul(inv_mat(getF(this,alpha))&
    ,this%I-Fexp),out)
End function

function getH(this, alpha) result(out)
type(linear_model_data) :: this
INTEGER :: alpha
real :: out(this%dim_obs,this%m),lambda_store(this%lambda_dim)
out=0
lambda_store=lambda(this,alpha)
! currently max 2 observations points are implemented, feel free to add more
 if (this%dim_obs>1) then
   out(this%dim_obs-1,1)=this%k_level ! [m3s-1]
 end if
out(this%dim_obs,this%m)=lambda_store(this%lambda_dim) ! [m3s-1]
End function

function getdelta(this)
type(linear_model_data) :: this
real :: getdelta(this%no_of_pars)
integer :: i
do i=1,this%no_of_pars      
  getdelta(i) = maxval(this%parameters(1:this%n_used,i))-&
    minval(this%parameters(1:this%n_used,i))
end do

end function

!this subroutine makes sure that we condition with design data randomly drawn
!from our generated sets
subroutine resample_pars(this)
  type(linear_model_data):: this
  integer :: nu,tm,dim
  tm=this%t_max
  nu=this%n_used
  dim=this%dim_obs

    IF (ALLOCATED(this%output)) THEN
      DEALLOCATE(this%output)
    ENDIF
    ALLOCATE(this%output(nu,tm*dim))
    this%output=0

    IF (ALLOCATED(this%parameters)) THEN
      DEALLOCATE(this%parameters)
    ENDIF
    ALLOCATE(this%parameters(nu+1,this%no_of_pars))
    this%parameters=0


end subroutine


FUNCTION lambda(this, alpha) result(out)
type(linear_model_data) :: this
INTEGER :: alpha
real :: width,slope,n_cat,n_pipe
real :: A0(this%lambda_dim),out(this%lambda_dim)
A0 = this%k_lam
width = this%parameters_physical(1)*this%parameters(alpha,2)
slope = this%parameters_physical(2)*this%parameters(alpha,3)
n_cat = this%parameters_physical(3)*this%parameters(alpha,4)
n_pipe = this%parameters_physical(4)*this%parameters(alpha,7)
!the following governs the linear resrvoirs representing the surface inputs
out=A0/n_pipe
out(1:this%input_dim)=A0(1:this%input_dim)*width/n_cat*sqrt(slope)
END FUNCTION lambda


function rho(Sigma, beta, gamma, par1, par2)
real :: par1(:), par2(:)
real :: beta(:),gamma
real :: Sigma(:,:),rho(size(Sigma,1),size(Sigma,1))
real :: r
integer :: D,q,j

rho=0
! D=1
! q=1

r=sum((beta*(par1-par2))**gamma)
! j=floor(D/2.0)+1+q
! rho = Sigma*max((1-r),0.0)**(j+1)*((j+1)*r+1)
! rho = Sigma*max((1-sum((beta*(par1-par2))**gamma)),0.0)**2
rho = Sigma*exp(-r)
! rho = Sigma*(1-sum((beta*(par1-par2))**gamma))
! write (*,*) (1-sum((beta*(par1-par2))**gamma))
end function rho




!functions for the c_emu implementation only, not restricted to linear model,
!no reason why it is in this file, just laziness

FUNCTION g(this,a,b) result(out)
IMPLICIT NONE
type(linear_model_data):: this
real :: out(this%m,this%m),integral(this%m,this%m)
INTEGER :: a,b
real :: beta(this%no_of_pars)

integral = matmul(-inv_mat(getF(this,a)+transpose(getF(this,b))),&
  (-expo_mat(getF(this,a)+transpose(getF(this,b)),this%delta_t,this%m)+this%I))

beta=1.0/(getdelta(this)*this%cor_factor_multi)
out=rho(this%sigma_ini,beta,this%gamma,this%parameters(a,:),this%parameters(b,:))*integral
END FUNCTION g


SUBROUTINE write_g_to_file_var(this)
IMPLICIT NONE
type(linear_model_data):: this
INTEGER :: a,b

IF (ALLOCATED(this%g_store_var)) THEN
    DEALLOCATE(this%g_store_var)
ENDIF
ALLOCATE(this%g_store_var(this%n_used+1,this%n_used+1,&
    this%m,this%m))
this%g_store_var=0

do a=1,this%n_used+1
  do b=1,this%n_used+1
      this%g_store_var(a,b,:,:)=g(this,a,b)
  end do
end do
END SUBROUTINE write_g_to_file_var

SUBROUTINE write_g_to_file(this)
IMPLICIT NONE
type(linear_model_data):: this
INTEGER :: a,b

IF (ALLOCATED(this%g_store)) THEN
    DEALLOCATE(this%g_store)
ENDIF
ALLOCATE(this%g_store(this%n_used,this%n_used,&
    this%m,this%m))
this%g_store=0

do a=1,this%n_used
  do b=1,this%n_used
      this%g_store(a,b,:,:)=g(this,a,b)
  end do
end do
END SUBROUTINE write_g_to_file

FUNCTION T(this,a,i,j)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,i,j
real :: T(this%m,this%m)

T=this%I
if (j .ge. i+2) then
  if (a .ne. this%T_indices_store(1)) then      
        IF (ALLOCATED(this%hT_store)) THEN
          DEALLOCATE(this%hT_store)
           ENDIF
          allocate(this%hT_store(this%m,this%m))
          this%hT_store=getFexp(this,a)
  end if
          this%T_store=matmul(this%T_store,transpose(this%hT_store))
          T=this%T_store
else if (j .lt. i+1) then
  T=0
    IF (ALLOCATED(this%T_store)) THEN
        DEALLOCATE(this%T_store)
    ENDIF
    allocate(this%T_store(this%m,this%m))
    this%T_store=this%I
end if

end function t
  

SUBROUTINE set_z(this)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,i,amount
real,allocatable :: z_til(:,:)
real :: k_store(this%m,this%t_max)
real :: h_store(this%m,this%m)
integer :: tm,dim,nu,m,nt

tm=this%t_max
dim=this%dim_obs
nu=this%n_used
nt=this%n_test_sets
m=this%m

  IF (ALLOCATED(z_til)) THEN
    DEALLOCATE(z_til)
  ENDIF
  ALLOCATE(z_til((nu+nt)*m,tm))
  z_til=0

  IF (ALLOCATED(this%z)) THEN
      DEALLOCATE(this%z)
  ENDIF
  allocate(this%z((nu)*tm*dim))
  this%z=0

  DO a=1,(nu)
    h_store=getFexp(this,a)
    do i=1,tm
      k_store(:,i)=getk(this,a,i,h_store)
    end do
    z_til(((a-1)*m+1):(a*m),1)=this%e_ini
    this%z((a-1)*dim*tm+1:(a-1)*dim*tm+dim)=&
      matmul(getH(this,a),z_til((a-1)*m+1:a*m,1))
    do i=2,tm
      z_til(((a-1)*m+1):(a*m),i)=matmul(h_store,&
        z_til((a-1)*m+1:a*m,i-1))+k_store(:,i)
      this%z((a-1)*tm*dim+(i-1)*dim+1:(a-1)*tm*dim+i*dim)=matmul(getH(this,a),&
        z_til((a-1)*m+1:a*m,i))
    end do
  END DO
END SUBROUTINE set_z

SUBROUTINE set_y(this)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,a_test,b,i
real :: y_til(this%n_test_sets*this%m,this%t_max)
real :: temp(this%m)
real :: k_store(this%m,this%t_max)
real :: h_store(this%m,this%m)
real :: g_store_local(this%n_used,this%m,this%m)
integer :: tm,dim,nu,m

tm=this%t_max
dim=this%dim_obs
nu=this%n_used
m=this%m

IF (ALLOCATED(this%y)) THEN
    DEALLOCATE(this%y)
ENDIF
allocate(this%y(this%n_test_sets*tm*dim))
this%y=0

do a=nu+1,nu+this%n_test_sets
    do b=1,this%n_used
        g_store_local(b,:,:)=g(this,a,b)
    end do
    a_test=a-nu
    h_store=getFexp(this,a)
    do i=1,tm
      k_store(:,i)=getk(this,a,i,h_store)
    end do
    temp=0
    
    do b=1,nu
        temp=temp+matmul(g_store_local(b,:,:),&
                this%z_prime((b-1)*tm*m+1:(b-1)*tm*m+m))
    end do
       ! y_til(((a_test-1)*m+1):(a_test*m),1)=temp+k_store(:,1)
        y_til(((a_test-1)*m+1):(a_test*m),1)=this%e_ini 
    this%y((a_test-1)*tm*dim+1:(a_test-1)*tm*dim+dim)=matmul(&
        (getH(this,a)),&
        y_til(((a_test-1)*m+1):(a_test*m),1))
    do i=2,tm
        temp=0
        do b=1,nu
          temp=temp+matmul(g_store_local(b,:,:),&
            this%z_prime((b-1)*tm*m+(i-1)*m+1:(b-1)*tm*m+i*m))
        end do
        y_til((a_test-1)*m+1:a_test*m,i)=&
            matmul(h_store,y_til((a_test-1)*m+1:&
                (a_test*m),i-1))&
            +k_store(:,i)+temp
        this%y((a_test-1)*tm*dim+(i-1)*dim+1:(a_test-1)*tm*dim+i*dim)=&
        matmul((geth(this,a)),y_til(((a_test-1)*m+1):(a_test*m),i))
    end do
end do
END SUBROUTINE set_y

SUBROUTINE set_z_prime(this)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,i,j
real :: z_back(this%n_used*this%t_max*this%dim_obs)
real :: output_long(this%t_max*this%n_used*this%dim_obs)
integer :: tm,dim,nu,m

tm=this%t_max
dim=this%dim_obs
nu=this%n_used
m=this%m
this%T_indices_store(3)=0

do a=1,this%n_used
  output_long((a-1)*tm*dim+1:a*tm*dim)=this%output(a,:)
end do
! z_back=matmul((this%sigma),&
!     output_long(1:nu*tm*dim)-&
!     this%z(1:tm*nu*dim))

! z_back=solve(this,this%sigma,output_long-this%z(1:tm*nu*dim),.false.)
z_back=solve_sparse(this,this%sigma,output_long-this%z(1:tm*nu*dim))
do a=1,this%n_used
    do i=1,this%t_max
        do j=1,this%t_max
            this%z_prime(m*tm*(a-1)+m*(i-1)+1:&
                         m*tm*(a-1)+m*i)=&
            this%z_prime(m*tm*(a-1)+m*(i-1)+1:&
                         m*tm*(a-1)+m*i)+&
                matmul(matmul(T(this,a,i,j),transpose(getH(this,a))),&
                z_back((a-1)*tm*dim+(j-1)*dim+1:(a-1)*tm*dim+j*dim))
       end do
    end do
end do
END SUBROUTINE set_z_prime



SUBROUTINE set_sigma_tilde(this)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,b,i,j
real :: h_store_a(this%m,this%m)
real :: h_store_b(this%m,this%m)
real :: h_store_b_transp(this%m,this%m)
real :: g_store(this%m,this%m)


IF (ALLOCATED(this%sigma_tilde)) THEN
DEALLOCATE(this%sigma_tilde)
ENDIF
ALLOCATE(this%sigma_tilde(this%n_used,&
    this%n_used,this%t_max,this%t_max,this%m,this%m))
this%sigma_tilde=0

DO a=1,this%n_used
  h_store_a=getFexp(this,a)
  DO b=1,this%n_used
  h_store_b=getFexp(this,b)
  h_store_b_transp=transpose(h_store_b)
    g_store=this%g_store(a,b,:,:)
    this%sigma_tilde(a,b,1,1,:,:)=g_store
    DO i=2,this%t_max
      DO j=1,i
        if (this%m .gt. 1) then
          IF(i .gt. j) THEN
            this%sigma_tilde(a,b,i,j,:,:)=&
              MATMUL(h_store_a,&
              this%sigma_tilde(a,b,i-1,j,:,:))
          ELSE
            this%sigma_tilde(a,b,i,j,:,:)=&
              (MATMUL(this%sigma_tilde(a,b,i,j-1,:,:),&
              h_store_b_transp)+&
              g_store)
         END IF  
        else !this takes care of the case, where I have m=1 (much faster)
          IF(i .gt. j) THEN
            this%sigma_tilde(a,b,i,j,1,1)=&
              h_store_a(1,1)*this%sigma_tilde(a,b,i-1,j,1,1)
          ELSE
            this%sigma_tilde(a,b,i,j,1,1)=&
              this%sigma_tilde(a,b,i,j-1,1,1)*h_store_b(1,1)+g_store(1,1)
         END IF 
        end if
     END DO
    END DO
  END DO
END DO

DO a=1,(this%n_used)
  DO b=1,(this%n_used)
    DO i=1,this%t_max
      DO j=(i+1),this%t_max
              this%sigma_tilde(a,b,i,j,:,:)=&
                  TRANSPOSE(this%sigma_tilde(b,a,j,i,:,:))
      END DO
    END DO
  ENDDO
ENDDO
END SUBROUTINE set_sigma_tilde

SUBROUTINE set_sigma(this)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,b,i,j
integer :: tm,dim,nu

tm=this%t_max
dim=this%dim_obs
nu=this%n_used

IF (ALLOCATED(this%sigma)) THEN
DEALLOCATE(this%sigma)
ENDIF
allocate(this%sigma(nu*tm*dim,nu*tm*dim))
this%sigma=0

DO a=1,nu
  DO b=1,a
    DO i=1,tm
      DO j=1,tm
        this%sigma((a-1)*tm*dim+(i-1)*dim+1:(a-1)*tm*dim+i*dim,&
          (b-1)*tm*dim+(j-1)*dim+1:(b-1)*tm*dim+j*dim)=&
          MATMUL(MATMUL(&
          getH(this,a),&
          this%sigma_tilde(a,b,i,j,:,:)),&
          transpose(getH(this,b)))
      END DO
    END DO
  END DO
END DO

DO a=1,nu
  DO b=a+1,nu
    DO i=1,tm
      DO j=1,tm
        this%sigma((a-1)*tm*dim+(i-1)*dim+1:(a-1)*tm*dim+i*dim,&
          (b-1)*tm*dim+(j-1)*dim+1:(b-1)*tm*dim+j*dim)=&
        (this%sigma((b-1)*tm*dim+(j-1)*dim+1:(b-1)*tm*dim+j*dim,&
          (a-1)*tm*dim+(i-1)*dim+1:(a-1)*tm*dim+i*dim))
      END DO
    END DO
  ENDDO
ENDDO

! this%sigma=inv_mat(this%sigma)

END SUBROUTINE set_sigma

FUNCTION solve(this,matrix,vector,calc_determinant) result(out)
type(linear_model_data):: this
INTEGER :: i,m,j
real :: matrix(:,:),vector(:),out(size(matrix,1))
INTEGER,ALLOCATABLE :: ipiv(:)
real,ALLOCATABLE :: work(:)
INTEGER :: info=0
logical :: calc_determinant

m=size(matrix,1)
out = vector
i =INT(m)
IF (ALLOCATED(work)) THEN
 DEALLOCATE(work)
ENDIF
ALLOCATE(work(i))
work=0
IF (ALLOCATED(ipiv)) THEN
 DEALLOCATE(ipiv)
ENDIF
ALLOCATE(ipiv(i))
ipiv=0

CALL ssytrf('U',i,matrix,i,ipiv,work,i,info)
IF ( info .NE. 0 ) THEN
 write(*,*) 'not OK'
END IF

if (calc_determinant .eqv. .true.) then
        this%det=0
        do j=1,m
                this%det=this%det+log(abs(matrix(j,j)))
        end do
end if
CALL ssytrs('U',i,1,matrix,i,ipiv,out,i,info)
IF ( info .NE. 0 ) THEN
 write(*,*) 'not OK'
END IF


end function

FUNCTION solve_sparse(this,matrix,vector) result(out)
implicit           real (a-h,o-z)
type(linear_model_data):: this
real :: matrix(:,:),vector(:),out(size(matrix,1))
logical ::            goodb, normal, precon

!     ------------------------------------------------------------------
!     test   solves sets up and solves a system (A - shift * I)x = b,
!     using Aprod to define A and Msolve to define a preconditioner.
!     ------------------------------------------------------------------

common    /mshift/ shiftm, pertm
                                                        
intrinsic          abs
external           aprod
external           msolve
logical ::            checkA
real ::   b(size(matrix,1)), r1(size(matrix,1)), r2(size(matrix,1)), v(size(matrix,1)), w(size(matrix,1))
real ::   x(size(matrix,1)), y(size(matrix,1)),  xtrue(size(matrix,1))
integer ::            n, nout, itnlim, istop, itn
parameter        ( one = 1.0,  two = 2.0 )

write(*,*) "start of CG"
precon = .FALSE. 
goodb = .FALSE.
shift  = 0.0
pertbn = 0.0
n=size(matrix,1)


shiftm = shift
pertm  = abs( pertbn )
nout   = 6
b=vector

checkA = .FALSE. 
itnlim = n * 2
! itnlim = 10
rtol   = 1.0E-6

call symmlq( n, b, r1, r2, v, w, x, y, &
aprod, msolve, checkA, goodb, precon, shift, &
nout , itnlim, rtol, &
istop, itn, anorm, acond, rnorm, ynorm,matrix )
out=x
write(*,*) "end of CG"

end function


FUNCTION solve_multidim(this,matrix,vector) result(out)
type(linear_model_data):: this
INTEGER :: i,m,j
real :: matrix(:,:),vector(:,:),out(size(vector,1),size(vector,2))
INTEGER,ALLOCATABLE :: ipiv(:)
real,ALLOCATABLE :: work(:)
INTEGER :: info=0

m=size(matrix,1)
out = vector
i =INT(m)
IF (ALLOCATED(work)) THEN
 DEALLOCATE(work)
ENDIF
ALLOCATE(work(i))
work=0
IF (ALLOCATED(ipiv)) THEN
 DEALLOCATE(ipiv)
ENDIF
ALLOCATE(ipiv(i))
ipiv=0

CALL ssytrf('U',i,matrix,i,ipiv,work,i,info)
IF ( info .NE. 0 ) THEN
 write(*,*) 'not OK'
END IF

CALL ssytrs('U',i,this%dim_obs*this%t_max,matrix,i,ipiv,out,i,info)
IF ( info .NE. 0 ) THEN
 write(*,*) 'not OK'
END IF


end function

! very inefficient calculation of variance, for plotting purposes
! a recursive formula should be created for a serious usage
function calc_variance(this) result(out)
IMPLICIT NONE
type(linear_model_data) :: this
INTEGER :: a,b,i,j
real :: h_store_a(this%m,this%m)
real :: h_store_b(this%m,this%m)
real :: h_store_b_transp(this%m,this%m)
real :: g_store(this%m,this%m)
real :: sigma_tilde_column_left(this%n_used,&
    this%t_max,this%t_max,this%m,this%m)
real :: sigma_tilde_column_right(this%n_used,&
    this%t_max,this%t_max,this%m,this%m)
real :: sigma_tilde_n1(this%t_max,this%t_max,this%m,this%m)
real :: sigma_left(this%t_max*this%dim_obs,this%t_max*this%dim_obs*this%n_used)
real :: sigma_right(this%t_max*this%dim_obs*this%n_used,this%t_max*this%dim_obs)
real :: sigma_n1(this%t_max*this%dim_obs,this%t_max*this%dim_obs)
real :: out(this%t_max*this%dim_obs,this%t_max*this%dim_obs)
integer :: tm,dim,nu


DO a=1,this%n_used
  h_store_a=getFexp(this,a)
  b=this%n_used+1
  h_store_b=getFexp(this,b)
  h_store_b_transp=transpose(h_store_b)
    g_store=this%g_store_var(a,b,:,:)
    sigma_tilde_column_right(a,1,1,:,:)=g_store
    DO i=2,this%t_max
      DO j=1,i
        if (this%m .gt. 1) then
          IF(i .gt. j) THEN
            sigma_tilde_column_right(a,i,j,:,:)=&
              MATMUL(h_store_a,&
              sigma_tilde_column_right(a,i-1,j,:,:))
          ELSE
            sigma_tilde_column_right(a,i,j,:,:)=&
              (MATMUL(sigma_tilde_column_right(a,i,j-1,:,:),&
              h_store_b_transp)+&
              g_store)
         END IF  
        else !this takes care of the case, where I have m=1 (much faster)
          IF(i .gt. j) THEN
            sigma_tilde_column_right(a,i,j,1,1)=&
              h_store_a(1,1)*sigma_tilde_column_right(a,i-1,j,1,1)
          ELSE
            sigma_tilde_column_right(a,i,j,1,1)=&
              sigma_tilde_column_right(a,i,j-1,1,1)*h_store_b(1,1)+g_store(1,1)
         END IF 
        end if
     END DO
    END DO
END DO

DO b=1,this%n_used
  h_store_a=getFexp(this,a)
  a=this%n_used+1
  h_store_b=getFexp(this,b)
  h_store_b_transp=transpose(h_store_b)
    g_store=this%g_store_var(a,b,:,:)
    sigma_tilde_column_left(b,1,1,:,:)=g_store
    DO i=2,this%t_max
      DO j=1,i
        if (this%m .gt. 1) then
          IF(i .gt. j) THEN
            sigma_tilde_column_left(b,i,j,:,:)=&
              MATMUL(h_store_a,&
              sigma_tilde_column_left(b,i-1,j,:,:))
          ELSE
            sigma_tilde_column_left(b,i,j,:,:)=&
              (MATMUL(sigma_tilde_column_left(b,i,j-1,:,:),&
              h_store_b_transp)+&
              g_store)
         END IF  
        else !this takes care of the case, where I have m=1 (much faster)
          IF(i .gt. j) THEN
            sigma_tilde_column_left(b,i,j,1,1)=&
              h_store_a(1,1)*sigma_tilde_column_left(b,i-1,j,1,1)
          ELSE
            sigma_tilde_column_left(b,i,j,1,1)=&
              sigma_tilde_column_left(b,i,j-1,1,1)*h_store_b(1,1)+g_store(1,1)
         END IF 
        end if
     END DO
    END DO
END DO


DO a=1,(this%n_used)
    DO i=1,this%t_max
      DO j=(i+1),this%t_max
              sigma_tilde_column_right(a,i,j,:,:)=&
                  sigma_tilde_column_left(a,j,i,:,:)
              sigma_tilde_column_left(a,i,j,:,:)=&
                  sigma_tilde_column_right(a,j,i,:,:)
      END DO
    END DO
ENDDO
  a=this%n_used+1
  b=this%n_used+1
  h_store_a=getFexp(this,a)
  h_store_b=getFexp(this,b)
  h_store_b_transp=transpose(h_store_b)
    g_store=this%g_store_var(a,b,:,:)
    sigma_tilde_n1(1,1,:,:)=g_store
    DO i=2,this%t_max
      DO j=1,i
        if (this%m .gt. 1) then
          IF(i .gt. j) THEN
            sigma_tilde_n1(i,j,:,:)=&
              MATMUL(h_store_a,&
              sigma_tilde_n1(i-1,j,:,:))
          ELSE
            sigma_tilde_n1(i,j,:,:)=&
              (MATMUL(sigma_tilde_n1(i,j-1,:,:),&
              h_store_b_transp)+&
              g_store)
         END IF  
        else !this takes care of the case, where I have m=1 (much faster)
          IF(i .gt. j) THEN
            sigma_tilde_n1(i,j,1,1)=&
              h_store_a(1,1)*sigma_tilde_n1(i-1,j,1,1)
          ELSE
            sigma_tilde_n1(i,j,1,1)=&
              sigma_tilde_n1(i,j-1,1,1)*h_store_b(1,1)+g_store(1,1)
         END IF 
        end if
     END DO
    END DO
    DO i=1,this%t_max
      DO j=(i+1),this%t_max
              sigma_tilde_n1(i,j,:,:)=&
                  sigma_tilde_n1(j,i,:,:)
      END DO
    END DO
tm=this%t_max
dim=this%dim_obs
nu=this%n_used
DO a=1,nu
    b=nu+1
    DO i=1,tm
      DO j=1,tm
        sigma_left((i-1)*dim+1:i*dim,(a-1)*tm*dim+(j-1)*dim+1:(a-1)*tm*dim+j*dim)=&
          MATMUL(MATMUL(&
          getH(this,b),&
          sigma_tilde_column_left(a,i,j,:,:)),&
          transpose(getH(this,a)))
      END DO
    END DO
END DO

DO a=1,nu
    b=nu+1
    DO i=1,tm
      DO j=1,tm
        sigma_right((a-1)*tm*dim+(i-1)*dim+1:(a-1)*tm*dim+i*dim,(j-1)*dim+1:j*dim)=&
          MATMUL(MATMUL(&
          getH(this,a),&
          sigma_tilde_column_right(a,i,j,:,:)),&
          transpose(getH(this,b)))
      END DO
    END DO
END DO
    b=nu+1
    a=nu+1
    DO i=1,tm
      DO j=1,tm
        sigma_n1((i-1)*dim+1:i*dim,(j-1)*dim+1:j*dim)=&
          MATMUL(MATMUL(&
          getH(this,a),&
          sigma_tilde_n1(i,j,:,:)),&
          transpose(getH(this,b)))
      END DO
    END DO
   out=sigma_n1

END function calc_variance



END MODULE
