PROGRAM CST_MULTI

IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Declear this part first !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER*4, PARAMETER :: n_layer = 1+0
INTEGER*4, PARAMETER :: GQ_Node = 20
REAL*8, PARAMETER :: inf_range = 5000.0d0
INTEGER*4, PARAMETER :: gauss_p = 200000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER*8 :: i, j, k, a, l, i_xi, j_xi, i_x, c_at, r, &
             INFO, point, IOSTAT, row, col, Max_x, Max_y, a_i, run_time, Mode, M_i
INTEGER*8, DIMENSION(1:3*n_layer)  :: IPIV
INTEGER*8, DIMENSION(1:6) :: IPIV02

REAL*8 :: Fact01, Fact02, DUM_X, h1, h2, step_xi, x, step_x, step_y, &
          Surf_F_range, y, DUM_y, ref_l, abs_x, m_y, x0, y0,&
          ref_lambda, s_eq, fr, intv, evalu_line
REAL*8, ALLOCATABLE :: con_x(:,:), con_y(:,:), con_z(:,:), a_run(:)
REAL*8, DIMENSION(1:n_layer,(Gauss_p-1)*GQ_Node,1) :: xi, zeta
REAL*8, DIMENSION(1,1:n_layer) :: Dept, y1, y2, omega_ini, omega,&
                                  lambda_ini, lambda, rho, h, y3, &
                                  length_s, mat_l
REAL*8, DIMENSION(1:25,1:25) :: C_iterN, X_iterN
REAL*16, PARAMETER :: PI=3.1415926535897932384626433832795028841971
COMPLEX*16, DIMENSION(1,1:n_layer) :: mu

COMPLEX*16, DIMENSION(1,1:gauss_p) :: F_p
COMPLEX*16, DIMENSION(n_layer,1:3,1) :: Yn
COMPLEX*16, DIMENSION(1:3,1:3) :: DUM_M
COMPLEX*16, DIMENSION(1:n_layer,1:3,1:3) :: M, N, M_inv, k_layer
COMPLEX*16, DIMENSION(3*n_layer,3*n_layer) :: K_assem
COMPLEX*16, DIMENSION(1:3*n_layer,1) :: Surf_F, DUM_Surf_F, Surf_F_set
COMPLEX*16, DIMENSION(n_layer,1:3,1) :: Disp_xi, DUM_Disp_xi
COMPLEX*16, DIMENSION((Gauss_p-1)*GQ_Node,n_layer,1:3,1) :: Cons_k
COMPLEX*16, DIMENSION((Gauss_p-1)*GQ_Node,1) :: u_x, u_y, w_z, s_xx, s_xy, s_yx, s_yy, c_xz, c_yz, s_zz
COMPLEX*16 :: GQ_u_x, GQ_u_y, GQ_w_z, GQ_s_xx, GQ_s_xy, GQ_s_yx, GQ_s_yy, GQ_c_xz, GQ_c_yz, GQ_s_zz
COMPLEX*16 :: SUM_u_x, SUM_u_y, SUM_w_z, SUM_s_xx, SUM_s_xy, SUM_s_yx, SUM_s_yy, SUM_c_xz, SUM_c_yz, SUM_s_zz
COMPLEX*16 :: u_x_phys, u_y_phys, w_z_phys, s_xx_phys, s_xy_phys, s_yx_phys, s_yy_phys, c_xz_phys, c_yz_phys, s_zz_phys, s_yy_physs

COMPLEX*16 ::zx1, zx2, zx3


OPEN(3, file = '5000x200000w1.0dlambtry.txt', status='UNKNOWN')

DO i=1,3*n_layer
    Surf_F(i,1) = 0.0
    Surf_F_set(i,1) = 0.0
END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!IF (y3(1,n_layer) < m_y) WRITE(*,*)
!IF (y3(1,n_layer) < m_y) WRITE(*,*) "WARNING! CHECK YOUR DEPT CAREFULLY."
!IF (y3(1,n_layer) < m_y) WRITE(*,*)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Change according to variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UNIFORM DISTRIBUTED LOAD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!run_time = 1
!ALLOCATE(a_run(run_time*2))

point = 201
fr = -20.0d0

evalu_line = 0.0d0

!intv = abs(fr)/(point-1)
    IF (point==1) THEN
        intv = 0.1d0
    ELSE
        intv = abs(fr)*2/(point-1)
    END IF
    

DO M_i=1,1
    IF (M_i == 1) Then
        Surf_F(2,1) = -1.0     ! Normal load   , p(xi)      !! Normalized by 2*mu
        Surf_F(1,1) =  0.0     ! Shear load    , q(xi)
        Surf_F(3,1) =  0.0     ! Moment        , m(xi)
        Surf_F_range = 1.0     ! Force acting on [-a,a] range
    ELSE IF (M_i == 2) Then
        Surf_F(2,1) =  0.0     ! Normal load   , p(xi)      !! Normalized by 2*mu
        Surf_F(1,1) = -1.0     ! Shear load    , q(xi)
        Surf_F(3,1) =  0.0     ! Moment        , m(xi)
        Surf_F_range = 1.0     ! Force acting on [-a,a] range
    ELSE IF (M_i == 3) Then
        Surf_F(2,1) =  0.0     ! Normal load   , p(xi)      !! Normalized by 2*mu
        Surf_F(1,1) =  0.0     ! Shear load    , q(xi)
        Surf_F(3,1) = -1.0     ! Moment        , m(xi)
        Surf_F_range = 1.0      ! Force acting on [-a,a] range
    ELSE
        STOP
    END IF
    

        !! ~~ Bi-mat
        !! mu_ini(1,2*i-1)     = a_run(a_i*2-1)
        !! mu_ini(1,2*i)       = a_run(a_i*2)



    DO i=1,n_layer
    Dept(1,i)           = 1.0d0/n_layer
    END DO
    
    !DO i=1,n_layer
    !Dept(1,1)           = 0.1d0
    !Dept(1,2)           = 0.1d0
    !Dept(1,3)           = 0.5d0
    !Dept(1,4)           = 0.5d0
    !Dept(1,n_layer)     = 1000.0d0
    !END DO
    
    DO i=1,n_layer
    lambda_ini(1,i) = 0.0d0
    END DO

    DO i=1,n_layer
    omega_ini(1,i)   = 1.0d0
    END DO

    DO i=1,n_layer
    length_s(1,i)   = 0.0001d0
    END DO
    
    DO i=1,n_layer
    rho(1,i)   =  1.0d0
    END DO
    
    DO i=1,n_layer
    mat_l(1,i)   =  0.0d0
    END DO
    
    DO i=1,n_layer
    mu(1,i)   = (1.0, 0.1d0)
    END DO
    
    !DO i=1,n_layer
    !mu(1,1)           = 2.0d0
    !mu(1,2)           = 2.0d0
    !mu(1,n_layer)     = 1.0d0
    !END DO
    
    ref_l = 1.00
    ref_lambda= 1.00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -------------------------------------------------------

    IF (gauss_p  <= 1) CALL ABORT

    y3(1,1) = 0.0
    DO i=1,n_layer
        y1(1,i) = 0.0
        y2(1,i) = Dept(1,i)
        IF (i/=1) THEN
            y3(1,i) = y3(1,i-1) + Dept(1,i)
        ELSE
            y3(1,i) = y3(1,i) + Dept(1,i)
        END IF
    END DO
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ALWAYS EVALUATE AT THE INTERFACE !!!!
    !evalu_line =y3(1,n_layer-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*)
x0 = 1.0
y0 = 1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO i=1,n_layer
        h(1,i) = Dept(1,i) / y0
    END DO

    step_xi   = (inf_range*2)/(gauss_p-1)
    
       
    CALL Gauss_Quad(GQ_Node,C_iterN,X_iterN)


Do i=1,n_layer
    omega(1,i)  = omega_ini(1,i)
    !rho(1,i) = (length_s(1,i)/ref_l)
    lambda(1,i)= (lambda_ini(1,i)/ref_lambda) 
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


WRITE(*,*)  '        rho                        Omega_bar                &
                                                 Total_dept                &
                    lambda_bar'
DO i=1,n_layer
    WRITE(*,*) rho(1,i), omega(1,i) , y3(1,i), lambda(1,i)
END DO
WRITE(*,*)
WRITE(*,*)

! ---------------------------------------------------------------------
!       start here !!
! ---------------------------------------------------------------------

    !#--------------------------------------#!
    !#------------- xi session -------------#!
    !#--------------------------------------#!


    DO i_xi = 1, gauss_p-1


            DO j_xi = 1, GQ_Node

                DO r = 1,3*n_layer
                    DO j = 1,3*n_layer
                        K_assem(r,j) = 0.0
                    END DO
                END DO

            DO i=1,n_layer !! Layer's here
            h1 = (-inf_range+(i_xi-1)*step_xi) * x0
            h2 = (-inf_range+(i_xi)*step_xi) * x0
            Fact01 = (h2-h1)/2
            Fact02 = (h2+h1)/2
            
            !xi(i,(i_xi-1)*GQ_Node+j_xi,1)  = h1

                xi(i,(i_xi-1)*GQ_Node+j_xi,1)  = (Fact01 * X_iterN(j_xi,GQ_Node) + Fact02)
                !zeta(i,(i_xi-1)*GQ_Node+j_xi,1)= sqrt(((xi(i,(i_xi-1)*GQ_Node+j_xi,1)**2)*(rho(1,i)**2))+1)

                IF (i==1) THEN
                    ! Uniformly distributed load
                    !Surf_F_set = (Surf_F/(CMPLX(0,1)*xi(i,(i_xi-1)*GQ_Node+j_xi,1)))* &
!(exp(CMPLX(0,1)*xi(i,(i_xi-1)*GQ_Node+j_xi,1)*(Surf_F_range)) - &
!exp(CMPLX(0,1)*xi(i,(i_xi-1)*GQ_Node+j_xi,1)*(-Surf_F_range)))
                    !point load
                    Surf_F_set = Surf_F 
                END IF
                    DUM_Surf_F = Surf_F_set
                    
                    !write(3,*) DUM_Surf_F


                    ! ---------------------------------------------------------------------
                    ! Generate M(xi) and N(xi) trhough the subroutine
                    ! ---------------------------------------------------------------------
                    ! BEWARE OF TRUCATIONAL ERROR OF CODING
          			
                        CALL M_N(n_layer, xi(i,(i_xi-1)*GQ_Node+j_xi,1), &
                                 y , h(1,i), rho(1,i),&
                                 omega(1,i), lambda(1,i), length_s(1,i), &
                                 mat_l(1,i), mu(1,i), Yn(i,1:3,1), &
                                 M(i,1:3,1:3), N(i,1:3,1:3))
                                 
                        !WRITE(3,*) xi(i,(i_xi-1)*GQ_Node+j_xi,1), Yn(i,1:3,1), M(i,1:3,1:3), N(i,1:3,1:3)
                        !WRITE(3,*) xi(i,(i_xi-1)*GQ_Node+j_xi,1), M(i,1:3,1:3)
                        
                    ! ---------------------------------------------------------------------
                    ! Matrix inversion of M(xi)
                    ! ---------------------------------------------------------------------
                             	                               
                                DUM_M(1:3,1:3) = M(i,1:3,1:3)
                                
                                
                                
                                M_inv(i,1:3,1:3) = inv(DUM_M(1:3,1:3))

				!write(3,*)  xi(i,(i_xi-1)*GQ_Node+j_xi,1), M_inv(i,1:3,1:3)
                                
                    ! ---------------------------------------------------------------------
                    ! Mutiplication of N(xi)*M_inv(xi)
                    ! Note that, K(xi)=N(xi)*(inv(M(xi)))
                    ! ---------------------------------------------------------------------

                            DO j = 1, 3
                                DO k = 1, 3
                                    k_layer(i,j,k) = 0.0           ! set zero in each element fisrt
                                    DO l =1, 3                     ! (row i of inv_M)*(col j of N)
                                        k_layer(i,j,k) = k_layer(i,j,k) + N(i,j,l)*M_inv(i,l,k)
                                        !write(3,*)  i, j, k,xi(i,(i_xi-1)*GQ_Node+j_xi,1), k_layer(i,j,k)
                                    END DO
                                END DO
                            END DO
                            
                            
                            

                    ! ---------------------------------------------------------------------
                    ! ASSEMBLE EACH k_layer
                    ! ---------------------------------------------------------------------
                    ! set zero to all element

                        ! Assembly !
                            IF (i == 1) THEN
                                IF (n_layer == 1) THEN
                                    DO j = 1, 3
                                        DO k = 1, 3
                                        K_assem(3*(i-1)+j,3*(i-1)+k) = K_assem(3*(i-1)+j,3*(i-1)+k) + k_layer(i,j,k)
                                        END DO
                                    END DO
                                ELSE
                                    DO j = 1, 6
                                        DO k = 1, 6
                                        K_assem(3*(i-1)+j,3*(i-1)+k) = K_assem(3*(i-1)+j,3*(i-1)+k) + k_layer(i,j,k)
                                        END DO
                                    END DO
                                END IF

                            ELSE IF (i < n_layer) THEN
                                DO j = 1, 6
                                    DO k = 1, 6
                                        IF (j <= 3) THEN
                                            K_assem(3*(i-1)+j,3*(i-1)+k) = K_assem(3*(i-1)+j,3*(i-1)+k) - k_layer(i,j,k)
                                        ELSE
                                            K_assem(3*(i-1)+j,3*(i-1)+k) = K_assem(3*(i-1)+j,3*(i-1)+k) + k_layer(i,j,k)
                                        END IF
                                    END DO
                                END DO

                            ELSE
                                DO j = 1, 3
                                    DO k = 1, 3
                                    K_assem(3*(i-1)+j,3*(i-1)+k) = K_assem(3*(i-1)+j,3*(i-1)+k) - k_layer(i,j,k)
                                    END DO
                                END DO
                            END IF

                    END DO ! n_layer
                    
                    
                                
                      
                    ! in some case of expression use:
                    ! REALPART(expr) -- IMAGPART(expr)
                    ! ---------------------------------------------------------------------
                    ! Solve the system of equation to get an unknow coefficients
                    ! SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
                    ! ---------------------------------------------------------------------
                    ! *** Need to set a A / B as well
                    ! A = K_assem
                    ! B = Surf_F
                    ! Ax=B

                    ! Get displacement at each interface
                    CALL ZGESV(SIZE(K_assem,1), SIZE(DUM_Surf_F,2),&
                     K_assem, SIZE(K_assem,1), &
                     IPIV, DUM_Surf_F, SIZE(DUM_Surf_F,1), INFO)
                     
                     !write(3,*) K_assem


                    DO i=1,n_layer

                        DO j=1,3
                            disp_xi(i,j,1) = 0.0
                        END DO
                       

                        IF (i==n_layer) THEN
                            DO j=1,3
                                disp_xi(i,j,1) = DUM_Surf_F(j+(i-1)*3,1)
                            END DO
                        ELSE
                            DO j=1,6
                                disp_xi(i,j,1) = DUM_Surf_F(j+(i-1)*3,1)
                            END DO
                            
                          
                        END IF
                        DUM_Disp_xi = disp_xi
                        
                        !write(3,*) DUM_Disp_xi
                        			

                        ! Get unknown coefficient for each layer
                        CALL ZGESV(SIZE(M(i,1:3,1:3),1),SIZE(DUM_Disp_xi,3),M(i,1:3,1:3),SIZE(M(i,1:3,1:3),1),&
                                    IPIV02,DUM_Disp_xi(i,:,:),SIZE(DUM_Disp_xi(i,:,:),1),INFO)

                        DO j=1,3
                            Cons_k((i_xi-1)*GQ_Node+j_xi,i,j,1) = DUM_Disp_xi(i,j,1)
                         
                        END DO

                    END DO

            END DO

    END DO


!#---------------------------------------------#!
!#------------- END of xi session -------------#!
!#---------------------------------------------#!
! Note that it is okay to bring this part outside the loop
! and using the same loop for the next part as well.
! Since it will be thae same, I will just continue using it.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X AND Y HERE
! ---------------------------------------------------------------------
! Get all field quatities at each (x,y) point.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(UNIT=10, FILE='output_data.csv', STATUS='REPLACE')
! เขียน Header ให้ Python รู้จัก
WRITE(10, *) "x,u_y_real,u_y_imag,u_x_real,u_x_imag"

DO i_x = 1, point

u_x_phys = 0.0
u_y_phys = 0.0
w_z_phys = 0.0
s_xx_phys = 0.0
s_xy_phys = 0.0
s_yx_phys = 0.0
s_yy_phys = 0.0
c_xz_phys = 0.0
c_yz_phys = 0.0
s_eq = 0.0


    !x =  1
    x = (fr + (i_x-1)*intv) /x0
    DUM_y = evalu_line

    DO i=1,n_layer
        IF  (DUM_y <= y3(1,i)) THEN
            C_at=i
            EXIT
        ELSE
            C_at=i
        END IF
    END DO

    IF (C_at==1) THEN
        y = DUM_y / y0
    ELSE
        y = (DUM_y - y3(1,C_at-1)) / y0
    END IF

    DO i_xi = 1, gauss_p-1
        h1 = (-inf_range+(i_xi-1)*step_xi) * x0
        h2 = (-inf_range+(i_xi)*step_xi) * x0
        Fact01 = (h2-h1)/2
        Fact02 = (h2+h1)/2

            SUM_u_x = 0.0
            SUM_u_y = 0.0
            SUM_w_z = 0.0
            SUM_s_xx = 0.0
            SUM_s_xy = 0.0
            SUM_s_yx = 0.0
            SUM_s_yy = 0.0
            SUM_c_xz = 0.0
            SUM_c_yz = 0.0
      

            DO j_xi = 1, GQ_Node
            
            

            CALL Field_q(xi(C_at,(i_xi-1)*GQ_Node+j_xi,1), y, &
                 omega(1,C_at), lambda(1,C_at),length_s(1,C_at), &
                 mat_l(1,C_at), mu(1,C_at),&
                 Cons_k((i_xi-1)*GQ_Node+j_xi,C_at,:,:), &
                 h(1,C_at), rho(1,C_at), &
                 Yn(i,1:3,1), &
                 u_x((i_xi-1)*GQ_Node+j_xi,1), &
                 u_y((i_xi-1)*GQ_Node+j_xi,1), &
                 w_z((i_xi-1)*GQ_Node+j_xi,1), &
                 s_xx((i_xi-1)*GQ_Node+j_xi,1), &
                 s_xy((i_xi-1)*GQ_Node+j_xi,1), &
                 s_yx((i_xi-1)*GQ_Node+j_xi,1), &
                 s_yy((i_xi-1)*GQ_Node+j_xi,1), &
                 c_xz((i_xi-1)*GQ_Node+j_xi,1), &
                 c_yz((i_xi-1)*GQ_Node+j_xi,1) )
                 !WRITE(3,*) Cons_k((i_xi-1)*GQ_Node+j_xi,C_at,:,:)
                 !WRITE(3,*) xi(C_at,(i_xi-1)*GQ_Node+j_xi,1), u_y((i_xi-1)*GQ_Node+j_xi,1), u_x((i_xi-1)*GQ_Node+j_xi,1)
		
            ! ---------------------------------------------------------------------
            ! Gauss-Quadrature.
            ! ---------------------------------------------------------------------
            GQ_u_x = (1.0/(2.0*pi)) * Fact01 * u_x((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_u_x = SUM_u_x + GQ_u_x

            GQ_u_y  = (1.0/(2.0*pi)) * Fact01 * u_y((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_u_y  = SUM_u_y + GQ_u_y

            GQ_w_z  = (1.0/(2.0*pi)) * Fact01 * w_z((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_w_z  = SUM_w_z + GQ_w_z

            GQ_s_xx = (1.0/(2.0*pi)) * Fact01 * s_xx((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_s_xx = SUM_s_xx + GQ_s_xx

            GQ_s_xy = (1.0/(2.0*pi)) * Fact01 * s_xy((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_s_xy = SUM_s_xy + GQ_s_xy

            GQ_s_yx = (1.0/(2.0*pi)) * Fact01 * s_yx((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_s_yx = SUM_s_yx + GQ_s_yx

            GQ_s_yy = (1.0/(2.0*pi)) * Fact01 * s_yy((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_s_yy = SUM_s_yy + GQ_s_yy

            GQ_c_xz = (1.0/(2.0*pi)) * Fact01 * c_xz((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_c_xz = SUM_c_xz + GQ_c_xz

            GQ_c_yz = (1.0/(2.0*pi)) * Fact01 * c_yz((i_xi-1)*GQ_Node+j_xi,1) * C_iterN(j_xi,GQ_Node) &
                    * EXP(-xi(C_at,(i_xi-1)*GQ_Node+j_xi,1) * x * CMPLX(0,1))
            SUM_c_yz = SUM_c_yz + GQ_c_yz

           

            END DO ! END of j_xi

        u_x_phys = u_x_phys + SUM_u_x
        u_y_phys = u_y_phys + SUM_u_y
        w_z_phys = w_z_phys + SUM_w_z
        s_xx_phys = s_xx_phys + SUM_s_xx
        s_xy_phys = s_xy_phys + SUM_s_xy
        s_yx_phys = s_yx_phys + SUM_s_yx
        s_yy_phys = s_yy_phys + SUM_s_yy
        c_xz_phys = c_xz_phys + SUM_c_xz
        c_yz_phys = c_yz_phys + SUM_c_yz
         
     

    END DO

s_eq = sqrt( 3.0 * ( 0.5 * (  &
       ( (REALPART(s_xx_phys) - (1.0/3.0)* (REALPART(s_xx_phys)+REALPART(s_yy_phys)) )**2.0 )+   &
       ( (REALPART(s_yy_phys) - (1.0/3.0)* (REALPART(s_xx_phys)+REALPART(s_yy_phys)) )**2.0 )+   &
       (1.0/4.0)*(REALPART(s_xy_phys)**2) + &
       (1.0/2.0)*REALPART(s_xy_phys)*REALPART(s_yx_phys) +    &
       (1.0/4.0)*(REALPART(s_yx_phys)**2) + &
       (1.0/2.0)*(rho(1,C_at)**2.0)*((REALPART(c_xz_phys)**2.0)+(REALPART(c_yz_phys)**2.0)))))

WRITE(*,*) int(C_at), int(i_x), real(x) , real(DUM_y), REALPART(s_yy_phys)


! เขียนข้อมูลในลูปของคุณ
! สมมติว่า x, u_y_phys, u_x_phys อยู่ในลูปนี้
WRITE(10, '(F10.4, ",", 4(E15.7, ","))') x, REAL(u_y_phys), AIMAG(u_y_phys), REAL(u_x_phys), AIMAG(u_x_phys)

END DO  ! จบลูป i_x


!WRITE(3,*) x, REALPART(u_y_phys), IMAGPART(u_y_phys), REALPART(u_x_phys), IMAGPART(u_x_phys)
!WRITE(3,*)  x, REALPART(s_yy_phys)
!REALPART(s_yx_phys), REALPART(c_xz_phys), REALPART(c_yz_phys)    
!WRITE(3,*)  x, , DUM_y, REALPART(u_x_phys), &
           !REALPART(u_y_phys), REALPART(w_z_phys), &
           !REALPART(s_xx_phys), REALPART(s_xy_phys), &
           !REALPART(s_yx_phys), REALPART(s_yy_phys), &
           !REALPART(c_xz_phys), REALPART(c_yz_phys), &
           ! s_eq

END DO
END DO


CLOSE(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                     !!
!!            Functions and Subroutines :)             !!
!!                                                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CONTAINS

! ---------------------------------------------------------------------
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.

    FUNCTION inv(A) result(Ainv)

    COMPLEX*16, DIMENSION(:,:), INTENT(IN) :: A
    COMPLEX*16, DIMENSION(size(A,1),size(A,2)) :: Ainv

    COMPLEX*16, DIMENSION(size(A,1)) :: work  ! work array for LAPACK
    INTEGER, DIMENSION(size(A,1)) :: ipiv   ! pivot indices
    INTEGER :: n, info

    ! External procedures defined in LAPACK
    EXTERNAL ZGETRF
    EXTERNAL ZGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    CALL ZGETRF(n, n, Ainv, n, ipiv, info)

    IF (info /= 0) THEN
    STOP 'Matrix is numerically singular!'
    END IF

    ! ZGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    CALL ZGETRI(n, Ainv, n, ipiv, work, n, info)

    IF (info /= 0) THEN
    STOP 'Matrix inversion failed!'
    END IF

    END FUNCTION inv

! ---------------------------------------------------------------------
    
! ---------------------------------------------------------------------
    SUBROUTINE M_N(n_layer, xi, y, h, rho, omega, lambda , length_s, &
    mat_l, mu, Yn, M, N)

    INTEGER*4, INTENT(IN) :: n_layer
    INTEGER*8 :: i, j
    REAL*8, INTENT(IN) :: xi, y, h, rho, omega, lambda ,length_s , &
    mat_l 
    COMPLEX*16, INTENT(IN) :: mu
    COMPLEX*16 ::zx1, zx2, zx3
    COMPLEX*16, DIMENSION(1:3,1), INTENT(OUT) :: Yn
    COMPLEX*16, DIMENSION(1:3,1:3), INTENT(OUT) :: M, N
    

    ! ------------------------------------------------
    !       Clear !
    ! ------------------------------------------------

    DO i=1,3
        DO j=1,3
            M(i,j) = 0.0
            N(i,j) = 0.0
        END DO
    END DO
    
    DO i=1,3
        Yn(i,1) = 0.0
    END DO
   
    zx1 = (xi**2-((omega**2)/(lambda+(mu*2))))
    Yn(1,1) = -sqrt(zx1)
    
    zx2 = (xi**2+(((1-(((mat_l*omega)**2)/(12*mu)))+(((1-&
     ((mat_l*omega)**2)/(12*mu))**2)+4*((omega*length_s)**2/mu))**(0.5))&
     /(2*(length_s**2))))
    Yn(2,1) = -sqrt(zx2)
    
    zx3 = (xi**2+(((1-(((mat_l*omega)**2)/(12*mu)))-(((1-&
     ((mat_l*omega)**2)/(12*mu))**2)+4*((omega*length_s)**2/mu))**(0.5))&
     /(2*(length_s**2))))
    Yn(3,1) = -sqrt(zx3)
    
    !WRITE(3,*) xi, omega, lambda, mu, zx1 , zx2, zx3, zx4, zx5, zx6
    !WRITE(3,*) xi, omega, lambda, mu, zx1, Yn(1,1) , Yn(2,1), Yn(3,1), Yn(4,1), Yn(5,1), Yn(6,1)
    
    ! ------------------------------------------------
    !       M(Xi) session
    ! ------------------------------------------------

    
    
    M(1,1) = CMPLX(0,1) * xi 
    
    M(1,2) = CMPLX(0,1) * Yn(2,1) 
    
    M(1,3) = CMPLX(0,1) * Yn(3,1) 
    
    
    M(2,1) = - Yn(1,1) 
    
    M(2,2) = - xi 
 
    M(2,3) = - xi 

    
    M(3,1) =  0.0
    
    M(3,2) = ( CMPLX(0,1) / 2 ) * (xi**2 - Yn(2,1)**2) 
    
    M(3,3) = ( CMPLX(0,1) / 2 ) * (xi**2 - Yn(3,1)**2) 
   
	
    ! ------------------------------------------------
    !       N(Xi) session
    ! ------------------------------------------------

    
    N(1,1) = CMPLX(0,1) * mu * xi * Yn(1,1) 
    
    N(1,2) = (CMPLX(0,1) / 2 ) * ( mu * ( Yn(2,1)**2 + xi**2 )&
- mu * (length_s**2) *( xi**2 - Yn(2,1)**2 )**2 &
+ ((((mat_l*omega)**2) / 12) * (xi**2 - Yn(2,1)**2) ))
    
    N(1,3) = (CMPLX(0,1) / 2 ) * ( mu * ( Yn(3,1)**2 + xi**2 ) - mu &
    * (length_s**2) *( xi**2 - Yn(3,1)**2 )**2 + &
    ((((mat_l*omega)**2) / 12) * (xi**2 - Yn(3,1)**2) ))

    
    N(2,1) =  ((lambda/2) * (xi**2) - (lambda/2 + mu) * (Yn(1,1)**2)) 
    
    N(2,2) = -mu * xi * Yn(2,1) 
    
    N(2,3) = -mu * xi * Yn(3,1) 

    
    N(3,1) = 0.0
    
    N(3,2) = CMPLX(0,1) * mu * (length_s**2) * Yn(2,1) * &
    ( xi**2 - Yn(2,1)**2 )
    
    N(3,3) = CMPLX(0,1) * mu * (length_s**2) * Yn(3,1) * &
    ( xi**2 - Yn(3,1)**2 )
    
    
    
    !write(3,*) xi, omega,lambda,mu, zx1
    !WRITE(3,*) xi, omega, lambda, mu, M(2,1) , M(2,2), M(2,3), M(2,4), M(2,5), M(2,6)
    !WRITE(3,*) xi, N(1,1), N(1,2), N(1,3)
    END SUBROUTINE M_N

! ---------------------------------------------------------------------

    SUBROUTINE Field_q(xi, y, omega, lambda, length_s, mat_l, mu,&
    Cons_k, h, rho, Yn, u_x, u_y, w_z, s_xx, s_xy, s_yx, s_yy, c_xz, c_yz)
    INTEGER*8 :: i                   
    REAL*8, INTENT(IN) :: xi, y, h, rho, omega, lambda, length_s, mat_l
    COMPLEX*16, INTENT(IN) :: mu
    COMPLEX*16, DIMENSION(1:3,1), INTENT(IN) :: Cons_k
    COMPLEX*16 ::zx1, zx2, zx3
    COMPLEX*16, DIMENSION(1:3,1), INTENT(OUT) :: Yn
    COMPLEX*16, INTENT(OUT) :: u_x, u_y, w_z, s_xx, s_xy, s_yx, s_yy, c_xz, c_yz
    
    DO i=1,3
        Yn(i,1) = 0.0
    END DO
    
    zx1 = (xi**2-((omega**2)/(lambda+(mu*2))))
    Yn(1,1) = -sqrt(zx1)
    
    zx2 = (xi**2+(((1-(((mat_l*omega)**2)/(12*mu)))+(((1-&
     ((mat_l*omega)**2)/(12*mu))**2)+4*((omega*length_s)**2/mu))**(0.5))&
     /(2*(length_s**2))))
    Yn(2,1) = -sqrt(zx2)
    
    zx3 = (xi**2+(((1-(((mat_l*omega)**2)/(12*mu)))-(((1-&
     ((mat_l*omega)**2)/(12*mu))**2)+4*((omega*length_s)**2/mu))**(0.5))&
     /(2*(length_s**2))))
    Yn(3,1) = -sqrt(zx3)
    
    !u_x =  CMPLX(0,1)*xi*EXP(Yn(1,1)*y)*Cons_k(1,1) &
    !+CMPLX(0,1)*Yn(2,1)*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    !+CMPLX(0,1)*Yn(3,1)*EXP(Yn(3,1)*y)*Cons_k(3,1)
    
    u_x = -(CMPLX(0,1)/mu)*xi*(2*xi**2-(omega**2*rho/mu)-2*&
    ((xi**2-(omega**2)*rho/(lambda+2*mu))**(0.5))*&
    ((xi**2-((omega**2)*rho/(mu)))**(0.5)))/((2*(xi**2)-((omega**2)*rho/(mu)))**2&
     - 4 * (xi**2)*((xi**2-(omega**2)*rho/(lambda+2*mu))**(0.5))&
    *((xi**2-((omega**2)*rho/(mu)))**(0.5)))
    
    !u_y = (-((omega**2)*rho/(mu))*(xi**2-((omega**2)*rho)/(lambda+2*mu))**(0.5)&
    !/((2*xi**2-((omega**2)*rho/(mu)))**2 - 4*(xi**2)*(((xi**2)-((omega**2)*rho)/&
    !(lambda+2*mu))*(xi**2-(omega**2)*rho/(mu)))**(0.5)))*(1/mu)
    
    
    
    u_y = -((((omega**2)/(mu))*((xi**2-(omega**2)&
    / (lambda+2*mu))**(0.5)))/(((2*(xi**2)-((omega**2)/(mu)))**2&
     - 4 * (xi**2)*((xi**2-(omega**2)/(lambda+2*mu))**(0.5))&
    *((xi**2-((omega**2)/(mu)))**(0.5)))*mu))
    
    !u_y = -(((omega**2)/(mu))*((xi**2-(omega**2)&
    !/ (lambda+2*mu))**(1/2)))/((2*(xi**2)-((omega**2)/(mu)))**2&
     !- 4 * (xi**2)*((xi**2-(omega**2)/(lambda+2*mu))**(1/2))&
    !*((xi**2-((omega**2)/(mu)))**(1/2)))
    
    !u_y = -Yn(1,1)*EXP(Yn(1,1)*y)*Cons_k(1,1) &
    !-xi*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    !-xi*EXP(Yn(3,1)*y)*Cons_k(3,1) 
    
    w_z = (CMPLX(0,1)/2)*(xi**2-Yn(2,1)**2)*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    +(CMPLX(0,1)/2)*(xi**2-Yn(3,1)**2)*EXP(Yn(3,1)*y)*Cons_k(3,1)
    
    s_xx = ((((lambda/2)+mu)*(xi**2))-((lambda/2)*(Yn(1,1)**2)))*&
    EXP(Yn(1,1)*y)*Cons_k(1,1) &
    +mu*xi*Yn(2,1)*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    +mu*xi*Yn(3,1)*EXP(Yn(3,1)*y)*Cons_k(3,1) 
    
    s_xy = CMPLX(0,1)*mu*xi*Yn(1,1)*EXP(Yn(1,1)*y)*Cons_k(1,1) &
    +(CMPLX(0,1)/2)*(mu*(xi**2+Yn(2,1)**2)+mu*(length_s**2)&
    *((xi**2-Yn(2,1)**2)**2)-((mat_l*omega)**2/12)*&
    (xi**2-Yn(2,1)**2))*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    +(CMPLX(0,1)/2)*(mu*(xi**2+Yn(3,1)**2)+mu*(length_s**2)&
    *((xi**2-Yn(3,1)**2)**2)-((mat_l*omega)**2/12)*&
    (xi**2-Yn(3,1)**2))*EXP(Yn(3,1)*y)*Cons_k(3,1) 

     s_yx = CMPLX(0,1)*mu*xi*Yn(1,1)*EXP(Yn(1,1)*y)*Cons_k(1,1)&
    +(CMPLX(0,1)/2)*(mu*(Yn(2,1)**2+xi**2)-mu*(length_s**2)&
    *((xi**2-Yn(2,1)**2)**2)+((mat_l*omega)**2/12)*&
    (xi**2-Yn(2,1)**2))*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    +(CMPLX(0,1)/2)*(mu*(Yn(3,1)**2+xi**2)-mu*(length_s**2)&
    *((xi**2-Yn(3,1)**2)**2)+((mat_l*omega)**2/12)*&
    (xi**2-Yn(3,1)**2))*EXP(Yn(3,1)*y)*Cons_k(3,1) 

    s_yy = ((lambda/2)*(xi**2)-((lambda/2)+mu)*(Yn(1,1)**2))&
    *EXP(Yn(1,1)*y)*Cons_k(1,1) &
    -mu*xi*Yn(2,1)*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    -mu*xi*Yn(3,1)*EXP(Yn(3,1)*y)*Cons_k(3,1) 
    
    c_xz = mu*(length_s**2)*xi*(xi**2-Yn(2,1)**2)*EXP(Yn(2,1)*y)*Cons_k(2,1) &
    +mu*(length_s**2)*xi*(xi**2-Yn(3,1)**2)*EXP(Yn(3,1)*y)*Cons_k(3,1) 

    c_yz = CMPLX(0,1) * mu * (length_s**2) * Yn(2,1) * (xi**2 - Yn(2,1)**2) * EXP(Yn(2,1)*y) * Cons_k(2,1) &
    +CMPLX(0,1)* mu* (length_s**2) * Yn(3,1) * (xi**2 - Yn(3,1)**2)&
    * EXP(Yn(3,1)*y) * Cons_k(3,1) 
    
    
    END SUBROUTINE Field_q

! ---------------------------------------------------------------------

    SUBROUTINE Gauss_Quad(GQ_Node,C_iterN,X_iterN)
    ! NOTE: C_iterN = Weights
    !       X_iterN = Abscissae
    ! Checked with https://keisan.casio.com/exec/system/1329114617

    INTEGER*4 :: i, j
    INTEGER*4, PARAMETER :: Max_node=25
    INTEGER*4, INTENT(IN) :: GQ_Node
    REAL*8, DIMENSION(Max_node,Max_node), INTENT(OUT) :: C_iterN, X_iterN

    IF (GQ_Node > Max_node) CALL ABORT

    DO i=1,25
    DO j=1,25
    C_iterN(i,j)    = 0.0
    X_iterN(i,j)    = 0.0
    END DO
    END DO

    !! NODE = 2
    C_iterN(1,2) =  1.000000000000000000
    C_iterN(2,2) =  1.000000000000000000
    X_iterN(1,2) = -0.577350269189625765
    X_iterN(2,2) =  0.577350269189625765

    !! NODE = 3
    C_iterN(1,3) =  0.555555555555555556
    C_iterN(2,3) =  0.888888888888888889
    C_iterN(3,3) =  0.555555555555555556
    X_iterN(1,3) = -0.774596669241483377
    X_iterN(2,3) =  0.000000000000000000
    X_iterN(3,3) =  0.774596669241483377


    !! NODE = 4
    C_iterN(1,4) =  0.34785484513745386
    C_iterN(2,4) =  0.652145154862546143
    C_iterN(3,4) =  0.652145154862546143
    C_iterN(4,4) =  0.347854845137453857
    X_iterN(1,4) = -0.861136311594052575
    X_iterN(2,4) = -0.339981043584856265
    X_iterN(3,4) =  0.339981043584856265
    X_iterN(4,4) =  0.861136311594052575

    !! NODE = 5
    C_iterN(1,5) =  0.236926885056189088
    C_iterN(2,5) =  0.478628670499366470
    C_iterN(3,5) =  0.568888888888888889
    C_iterN(4,5) =  0.478628670499366468
    C_iterN(5,5) =  0.236926885056189088
    X_iterN(1,5) = -0.906179845938663993
    X_iterN(2,5) = -0.538469310105683091
    X_iterN(3,5) =  0.000000000000000000
    X_iterN(4,5) =  0.538469310105683091
    X_iterN(5,5) =  0.906179845938663993

    !! NODE = 6
    C_iterN(1,6) =  0.171324492379170345
    C_iterN(2,6) =  0.360761573048138608
    C_iterN(3,6) =  0.46791393457269105
    C_iterN(4,6) =  0.46791393457269105
    C_iterN(5,6) =  0.360761573048138608
    C_iterN(6,6) =  0.171324492379170345
    X_iterN(1,6) = -0.932469514203152028
    X_iterN(2,6) = -0.661209386466264514
    X_iterN(3,6) = -0.238619186083196909
    X_iterN(4,6) =  0.238619186083196909
    X_iterN(5,6) =  0.66120938646626451
    X_iterN(6,6) =  0.932469514203152028

    !! NODE = 7
    C_iterN(1,7) =  0.129484966168869693
    C_iterN(2,7) =  0.279705391489276668
    C_iterN(3,7) =  0.381830050505118945
    C_iterN(4,7) =  0.417959183673469388
    C_iterN(5,7) =  0.38183005050511895
    C_iterN(6,7) =  0.279705391489276668
    C_iterN(7,7) =  0.129484966168869693
    X_iterN(1,7) = -0.949107912342758525
    X_iterN(2,7) = -0.74153118559939444
    X_iterN(3,7) = -0.405845151377397167
    X_iterN(4,7) =  0.000000000000000000
    X_iterN(5,7) =  0.405845151377397167
    X_iterN(6,7) =  0.74153118559939444
    X_iterN(7,7) =  0.949107912342758525



    !! NODE = 8
    C_iterN(1,8) =  0.101228536290376259
    C_iterN(2,8) =  0.22238103445337447
    C_iterN(3,8) =  0.313706645877887287
    C_iterN(4,8) =  0.362683783378361983
    C_iterN(5,8) =  0.36268378337836198
    C_iterN(6,8) =  0.313706645877887287
    C_iterN(7,8) =  0.222381034453374471
    C_iterN(8,8) =  0.10122853629037626
    X_iterN(1,8) = -0.960289856497536232
    X_iterN(2,8) = -0.79666647741362674
    X_iterN(3,8) = -0.525532409916328986
    X_iterN(4,8) = -0.183434642495649805
    X_iterN(5,8) =  0.183434642495649805
    X_iterN(6,8) =  0.525532409916328986
    X_iterN(7,8) =  0.79666647741362674
    X_iterN(8,8) =  0.960289856497536232

    !! NODE = 9
    C_iterN(1,9) =  0.081274388361574412
    C_iterN(2,9) =  0.1806481606948574
    C_iterN(3,9) =  0.260610696402935462
    C_iterN(4,9) =  0.31234707704000284
    C_iterN(5,9) =  0.33023935500125976
    C_iterN(6,9) =  0.31234707704000284
    C_iterN(7,9) =  0.26061069640293546
    C_iterN(8,9) =  0.180648160694857404
    C_iterN(9,9) =  0.081274388361574412
    X_iterN(1,9) = -0.96816023950762609
    X_iterN(2,9) = -0.836031107326635794
    X_iterN(3,9) = -0.613371432700590397
    X_iterN(4,9) = -0.32425342340380893
    X_iterN(5,9) =  0.000000000000000000
    X_iterN(6,9) =  0.32425342340380893
    X_iterN(7,9) =  0.613371432700590397
    X_iterN(8,9) =  0.836031107326635794
    X_iterN(9,9) =  0.96816023950762609

    !! NODE = 10
    C_iterN(1,10) =  0.066671344308688138
    C_iterN(2,10) =  0.14945134915058059
    C_iterN(3,10) =  0.219086362515982044
    C_iterN(4,10) =  0.26926671930999636
    C_iterN(5,10) =  0.29552422471475287
    C_iterN(6,10) =  0.29552422471475287
    C_iterN(7,10) =  0.269266719309996355
    C_iterN(8,10) =  0.21908636251598204
    C_iterN(9,10) =  0.14945134915058059
    C_iterN(10,10) =  0.066671344308688138
    X_iterN(1,10) = -0.97390652851717172
    X_iterN(2,10) = -0.86506336668898451
    X_iterN(3,10) = -0.679409568299024406
    X_iterN(4,10) = -0.433395394129247191
    X_iterN(5,10) = -0.14887433898163121
    X_iterN(6,10) =  0.148874338981631211
    X_iterN(7,10) =  0.433395394129247191
    X_iterN(8,10) =  0.679409568299024406
    X_iterN(9,10) =  0.865063366688984511
    X_iterN(10,10) =  0.97390652851717172


    !! NODE = 11
    C_iterN(1,11) =  0.055668567116173666
    C_iterN(2,11) =  0.125580369464904625
    C_iterN(3,11) =  0.186290210927734251
    C_iterN(4,11) =  0.23319376459199048
    C_iterN(5,11) =  0.26280454451024666
    C_iterN(6,11) =  0.27292508677790063
    C_iterN(7,11) =  0.262804544510246662
    C_iterN(8,11) =  0.23319376459199048
    C_iterN(9,11) =  0.18629021092773425
    C_iterN(10,11) =  0.12558036946490462
    C_iterN(11,11) =  0.055668567116173666
    X_iterN(1,11) = -0.978228658146056993
    X_iterN(2,11) = -0.887062599768095299
    X_iterN(3,11) = -0.730152005574049324
    X_iterN(4,11) = -0.519096129206811816
    X_iterN(5,11) = -0.26954315595234497
    X_iterN(6,11) =  0
    X_iterN(7,11) =  0.269543155952344972
    X_iterN(8,11) =  0.519096129206811816
    X_iterN(9,11) =  0.730152005574049324
    X_iterN(10,11) =  0.887062599768095299
    X_iterN(11,11) =  0.978228658146056993

    !! NODE = 12
    C_iterN(1,12) =  0.047175336386511827
    C_iterN(2,12) =  0.10693932599531843
    C_iterN(3,12) =  0.160078328543346226
    C_iterN(4,12) =  0.203167426723065922
    C_iterN(5,12) =  0.23349253653835481
    C_iterN(6,12) =  0.249147045813402785
    C_iterN(7,12) =  0.24914704581340279
    C_iterN(8,12) =  0.233492536538354809
    C_iterN(9,12) =  0.203167426723065922
    C_iterN(10,12) =  0.16007832854334623
    C_iterN(11,12) =  0.106939325995318431
    C_iterN(12,12) =  0.047175336386511827
    X_iterN(1,12) = -0.981560634246719251
    X_iterN(2,12) = -0.904117256370474857
    X_iterN(3,12) = -0.769902674194304687
    X_iterN(4,12) = -0.587317954286617447
    X_iterN(5,12) = -0.367831498998180194
    X_iterN(6,12) = -0.125233408511468916
    X_iterN(7,12) =  0.125233408511468916
    X_iterN(8,12) =  0.367831498998180194
    X_iterN(9,12) =  0.587317954286617447
    X_iterN(10,12) =  0.769902674194304687
    X_iterN(11,12) =  0.904117256370474857
    X_iterN(12,12) =  0.981560634246719251

    !! NODE = 13
    C_iterN(1,13) =  0.04048400476531588
    C_iterN(2,13) =  0.092121499837728448
    C_iterN(3,13) =  0.138873510219787239
    C_iterN(4,13) =  0.178145980761945738
    C_iterN(5,13) =  0.207816047536888502
    C_iterN(6,13) =  0.22628318026289724
    C_iterN(7,13) =  0.23255155323087391
    C_iterN(8,13) =  0.22628318026289724
    C_iterN(9,13) =  0.207816047536888502
    C_iterN(10,13) =  0.17814598076194574
    C_iterN(11,13) =  0.138873510219787239
    C_iterN(12,13) =  0.09212149983772845
    C_iterN(13,13) =  0.04048400476531588
    X_iterN(1,13) = -0.98418305471858815
    X_iterN(2,13) = -0.917598399222977965
    X_iterN(3,13) = -0.801578090733309913
    X_iterN(4,13) = -0.642349339440340221
    X_iterN(5,13) = -0.448492751036446853
    X_iterN(6,13) = -0.230458315955134794
    X_iterN(7,13) =  0
    X_iterN(8,13) =  0.23045831595513479
    X_iterN(9,13) =  0.44849275103644685
    X_iterN(10,13) =  0.642349339440340221
    X_iterN(11,13) =  0.801578090733309913
    X_iterN(12,13) =  0.917598399222977965
    X_iterN(13,13) =  0.98418305471858815


    !! NODE = 14
    C_iterN(1,14) =  0.035119460331751863
    C_iterN(2,14) =  0.08015808715976021
    C_iterN(3,14) =  0.121518570687903185
    C_iterN(4,14) =  0.157203167158193535
    C_iterN(5,14) =  0.185538397477937814
    C_iterN(6,14) =  0.205198463721295604
    C_iterN(7,14) =  0.21526385346315779
    C_iterN(8,14) =  0.21526385346315779
    C_iterN(9,14) =  0.205198463721295604
    C_iterN(10,14) =  0.18553839747793781
    C_iterN(11,14) =  0.157203167158193535
    C_iterN(12,14) =  0.121518570687903185
    C_iterN(13,14) =  0.08015808715976021
    C_iterN(14,14) =  0.035119460331751863
    X_iterN(1,14) = -0.986283808696812339
    X_iterN(2,14) = -0.928434883663573517
    X_iterN(3,14) = -0.827201315069764993
    X_iterN(4,14) = -0.68729290481168547
    X_iterN(5,14) = -0.515248636358154092
    X_iterN(6,14) = -0.31911236892788976
    X_iterN(7,14) = -0.108054948707343662
    X_iterN(8,14) =  0.108054948707343662
    X_iterN(9,14) =  0.31911236892788976
    X_iterN(10,14) =  0.515248636358154092
    X_iterN(11,14) =  0.68729290481168547
    X_iterN(12,14) =  0.827201315069764993
    X_iterN(13,14) =  0.928434883663573517
    X_iterN(14,14) =  0.986283808696812339


    !! NODE = 15
    C_iterN(1,15) =  0.0307532419961172684
    C_iterN(2,15) =  0.0703660474881081247
    C_iterN(3,15) =  0.107159220467171935
    C_iterN(4,15) =  0.139570677926154314
    C_iterN(5,15) =  0.166269205816993934
    C_iterN(6,15) =  0.18616100001556221
    C_iterN(7,15) =  0.198431485327111577
    C_iterN(8,15) =  0.202578241925561273
    C_iterN(9,15) =  0.19843148532711158
    C_iterN(10,15) =  0.18616100001556221
    C_iterN(11,15) =  0.166269205816993934
    C_iterN(12,15) =  0.13957067792615431
    C_iterN(13,15) =  0.10715922046717194
    C_iterN(14,15) =  0.07036604748810812
    C_iterN(15,15) =  0.030753241996117268
    X_iterN(1,15) = -0.987992518020485429
    X_iterN(2,15) = -0.9372733924007059
    X_iterN(3,15) = -0.848206583410427216
    X_iterN(4,15) = -0.724417731360170047
    X_iterN(5,15) = -0.570972172608538848
    X_iterN(6,15) = -0.39415134707756337
    X_iterN(7,15) = -0.201194093997434522
    X_iterN(8,15) =  0
    X_iterN(9,15) =  0.201194093997434522
    X_iterN(10,15) =  0.39415134707756337
    X_iterN(11,15) =  0.570972172608538848
    X_iterN(12,15) =  0.724417731360170047
    X_iterN(13,15) =  0.848206583410427216
    X_iterN(14,15) =  0.937273392400705904
    X_iterN(15,15) =  0.987992518020485429


    !! NODE = 16
    C_iterN(1,16) =  0.027152459411754095
    C_iterN(2,16) =  0.062253523938647893
    C_iterN(3,16) =  0.095158511682492785
    C_iterN(4,16) =  0.12462897125553387
    C_iterN(5,16) =  0.149595988816576732
    C_iterN(6,16) =  0.169156519395002538
    C_iterN(7,16) =  0.18260341504492359
    C_iterN(8,16) =  0.1894506104550685
    C_iterN(9,16) =  0.189450610455068496
    C_iterN(10,16) =  0.182603415044923589
    C_iterN(11,16) =  0.169156519395002538
    C_iterN(12,16) =  0.149595988816576732
    C_iterN(13,16) =  0.12462897125553387
    C_iterN(14,16) =  0.09515851168249278
    C_iterN(15,16) =  0.062253523938647893
    C_iterN(16,16) =  0.02715245941175409
    X_iterN(1,16) = -0.989400934991649933
    X_iterN(2,16) = -0.944575023073232576
    X_iterN(3,16) = -0.865631202387831744
    X_iterN(4,16) = -0.75540440835500303
    X_iterN(5,16) = -0.61787624440264375
    X_iterN(6,16) = -0.458016777657227386
    X_iterN(7,16) = -0.281603550779258913
    X_iterN(8,16) = -0.09501250983763744
    X_iterN(9,16) =  0.09501250983763744
    X_iterN(10,16) =  0.281603550779258913
    X_iterN(11,16) =  0.458016777657227386
    X_iterN(12,16) =  0.617876244402643748
    X_iterN(13,16) =  0.755404408355003034
    X_iterN(14,16) =  0.86563120238783174
    X_iterN(15,16) =  0.944575023073232576
    X_iterN(16,16) =  0.989400934991649933

    !! NODE = 17
    C_iterN(1,17) =  0.024148302868547932
    C_iterN(2,17) =  0.055459529373987201
    C_iterN(3,17) =  0.085036148317179181
    C_iterN(4,17) =  0.11188384719340397
    C_iterN(5,17) =  0.135136368468525473
    C_iterN(6,17) =  0.15404576107681029
    C_iterN(7,17) =  0.168004102156450045
    C_iterN(8,17) =  0.176562705366992646
    C_iterN(9,17) =  0.179446470356206526
    C_iterN(10,17) =  0.17656270536699265
    C_iterN(11,17) =  0.168004102156450045
    C_iterN(12,17) =  0.154045761076810288
    C_iterN(13,17) =  0.13513636846852547
    C_iterN(14,17) =  0.111883847193403971
    C_iterN(15,17) =  0.085036148317179181
    C_iterN(16,17) =  0.055459529373987201
    C_iterN(17,17) =  0.024148302868547932
    X_iterN(1,17) = -0.990575475314417336
    X_iterN(2,17) = -0.950675521768767761
    X_iterN(3,17) = -0.880239153726985902
    X_iterN(4,17) = -0.78151400389680141
    X_iterN(5,17) = -0.657671159216690766
    X_iterN(6,17) = -0.512690537086476968
    X_iterN(7,17) = -0.351231763453876315
    X_iterN(8,17) = -0.178484181495847856
    X_iterN(9,17) =  0
    X_iterN(10,17) =  0.178484181495847856
    X_iterN(11,17) =  0.351231763453876315
    X_iterN(12,17) =  0.512690537086476968
    X_iterN(13,17) =  0.657671159216690766
    X_iterN(14,17) =  0.781514003896801407
    X_iterN(15,17) =  0.880239153726985902
    X_iterN(16,17) =  0.950675521768767761
    X_iterN(17,17) =  0.99057547531441734


    !! NODE = 18
    C_iterN(1,18) =  0.02161601352648331
    C_iterN(2,18) =  0.0497145488949698
    C_iterN(3,18) =  0.076425730254889057
    C_iterN(4,18) =  0.100942044106287166
    C_iterN(5,18) =  0.12255520671147846
    C_iterN(6,18) =  0.140642914670650651
    C_iterN(7,18) =  0.154684675126265245
    C_iterN(8,18) =  0.164276483745832723
    C_iterN(9,18) =  0.169142382963143592
    C_iterN(10,18) =  0.169142382963143592
    C_iterN(11,18) =  0.16427648374583272
    C_iterN(12,18) =  0.15468467512626524
    C_iterN(13,18) =  0.140642914670650651
    C_iterN(14,18) =  0.12255520671147846
    C_iterN(15,18) =  0.100942044106287166
    C_iterN(16,18) =  0.07642573025488906
    C_iterN(17,18) =  0.049714548894969796
    C_iterN(18,18) =  0.02161601352648331
    X_iterN(1,18) = -0.991565168420930947
    X_iterN(2,18) = -0.955823949571397755
    X_iterN(3,18) = -0.892602466497555739
    X_iterN(4,18) = -0.803704958972523116
    X_iterN(5,18) = -0.69168704306035321
    X_iterN(6,18) = -0.559770831073947535
    X_iterN(7,18) = -0.411751161462842646
    X_iterN(8,18) = -0.25188622569150551
    X_iterN(9,18) = -0.084775013041735301
    X_iterN(10,18) =  0.084775013041735301
    X_iterN(11,18) =  0.25188622569150551
    X_iterN(12,18) =  0.411751161462842646
    X_iterN(13,18) =  0.559770831073947535
    X_iterN(14,18) =  0.691687043060353208
    X_iterN(15,18) =  0.803704958972523116
    X_iterN(16,18) =  0.892602466497555739
    X_iterN(17,18) =  0.955823949571397755
    X_iterN(18,18) =  0.991565168420930947

    !! NODE = 19
    C_iterN(1,19) =  0.019461788229726477
    C_iterN(2,19) =  0.0448142267656996
    C_iterN(3,19) =  0.069044542737641227
    C_iterN(4,19) =  0.091490021622449999
    C_iterN(5,19) =  0.111566645547333995
    C_iterN(6,19) =  0.128753962539336228
    C_iterN(7,19) =  0.142606702173606612
    C_iterN(8,19) =  0.152766042065859667
    C_iterN(9,19) =  0.15896884339395435
    C_iterN(10,19) =  0.161054449848783696
    C_iterN(11,19) =  0.158968843393954348
    C_iterN(12,19) =  0.152766042065859667
    C_iterN(13,19) =  0.14260670217360661
    C_iterN(14,19) =  0.128753962539336228
    C_iterN(15,19) =  0.111566645547333995
    C_iterN(16,19) =  0.091490021622449999
    C_iterN(17,19) =  0.06904454273764123
    C_iterN(18,19) =  0.0448142267656996
    C_iterN(19,19) =  0.019461788229726477
    X_iterN(1,19) = -0.992406843843584403
    X_iterN(2,19) = -0.960208152134830031
    X_iterN(3,19) = -0.903155903614817902
    X_iterN(4,19) = -0.82271465653714283
    X_iterN(5,19) = -0.720966177335229379
    X_iterN(6,19) = -0.600545304661681024
    X_iterN(7,19) = -0.464570741375960946
    X_iterN(8,19) = -0.316564099963629832
    X_iterN(9,19) = -0.160358645640225376
    X_iterN(10,19) =  0
    X_iterN(11,19) =  0.160358645640225376
    X_iterN(12,19) =  0.316564099963629832
    X_iterN(13,19) =  0.464570741375960946
    X_iterN(14,19) =  0.600545304661681024
    X_iterN(15,19) =  0.720966177335229379
    X_iterN(16,19) =  0.822714656537142825
    X_iterN(17,19) =  0.903155903614817902
    X_iterN(18,19) =  0.960208152134830031
    X_iterN(19,19) =  0.992406843843584403


    !! NODE = 20
    C_iterN(1,20) =  0.017614007139152118
    C_iterN(2,20) =  0.040601429800386941
    C_iterN(3,20) =  0.062672048334109064
    C_iterN(4,20) =  0.083276741576704749
    C_iterN(5,20) =  0.10193011981724044
    C_iterN(6,20) =  0.118194531961518417
    C_iterN(7,20) =  0.131688638449176627
    C_iterN(8,20) =  0.142096109318382051
    C_iterN(9,20) =  0.149172986472603747
    C_iterN(10,20) =  0.152753387130725851
    C_iterN(11,20) =  0.15275338713072585
    C_iterN(12,20) =  0.149172986472603747
    C_iterN(13,20) =  0.142096109318382051
    C_iterN(14,20) =  0.13168863844917663
    C_iterN(15,20) =  0.118194531961518417
    C_iterN(16,20) =  0.10193011981724044
    C_iterN(17,20) =  0.083276741576704749
    C_iterN(18,20) =  0.062672048334109064
    C_iterN(19,20) =  0.040601429800386941
    C_iterN(20,20) =  0.0176140071391521183
    X_iterN(1,20) = -0.993128599185094925
    X_iterN(2,20) = -0.963971927277913791
    X_iterN(3,20) = -0.912234428251325906
    X_iterN(4,20) = -0.839116971822218823
    X_iterN(5,20) = -0.746331906460150793
    X_iterN(6,20) = -0.636053680726515026
    X_iterN(7,20) = -0.510867001950827098
    X_iterN(8,20) = -0.373706088715419561
    X_iterN(9,20) = -0.227785851141645078
    X_iterN(10,20) = -0.0765265211334973338
    X_iterN(11,20) =  0.076526521133497334
    X_iterN(12,20) =  0.227785851141645078
    X_iterN(13,20) =  0.373706088715419561
    X_iterN(14,20) =  0.510867001950827098
    X_iterN(15,20) =  0.636053680726515026
    X_iterN(16,20) =  0.746331906460150793
    X_iterN(17,20) =  0.839116971822218823
    X_iterN(18,20) =  0.912234428251325906
    X_iterN(19,20) =  0.963971927277913791
    X_iterN(20,20) =  0.993128599185094925

    !! NODE = 21
    C_iterN(1,21) =  0.0160172282577743333
    C_iterN(2,21) =  0.036953789770852494
    C_iterN(3,21) =  0.057134425426857208
    C_iterN(4,21) =  0.076100113628379302
    C_iterN(5,21) =  0.093444423456033862
    C_iterN(6,21) =  0.10879729916714838
    C_iterN(7,21) =  0.121831416053728534
    C_iterN(8,21) =  0.132268938633337462
    C_iterN(9,21) =  0.139887394791073155
    C_iterN(10,21) =  0.14452440398997006
    C_iterN(11,21) =  0.146081133649690427
    C_iterN(12,21) =  0.144524403989970059
    C_iterN(13,21) =  0.139887394791073155
    C_iterN(14,21) =  0.13226893863333746
    C_iterN(15,21) =  0.12183141605372853
    C_iterN(16,21) =  0.108797299167148378
    C_iterN(17,21) =  0.093444423456033862
    C_iterN(18,21) =  0.0761001136283793
    C_iterN(19,21) =  0.05713442542685721
    C_iterN(20,21) =  0.036953789770852494
    C_iterN(21,21) =  0.0160172282577743333
    X_iterN(1,21) = -0.9937521706203895
    X_iterN(2,21) = -0.967226838566306294
    X_iterN(3,21) = -0.920099334150400829
    X_iterN(4,21) = -0.853363364583317284
    X_iterN(5,21) = -0.768439963475677909
    X_iterN(6,21) = -0.667138804197412319
    X_iterN(7,21) = -0.551618835887219807
    X_iterN(8,21) = -0.424342120207438784
    X_iterN(9,21) = -0.288021316802401097
    X_iterN(10,21) = -0.14556185416089509
    X_iterN(11,21) =  0
    X_iterN(12,21) =  0.145561854160895091
    X_iterN(13,21) =  0.288021316802401097
    X_iterN(14,21) =  0.424342120207438784
    X_iterN(15,21) =  0.551618835887219807
    X_iterN(16,21) =  0.667138804197412319
    X_iterN(17,21) =  0.768439963475677909
    X_iterN(18,21) =  0.853363364583317284
    X_iterN(19,21) =  0.920099334150400829
    X_iterN(20,21) =  0.967226838566306294
    X_iterN(21,21) =  0.9937521706203895



    !! NODE = 22

    C_iterN(1,22) =  0.014627995298272201
    C_iterN(2,22) =  0.033774901584814155
    C_iterN(3,22) =  0.0522933351526832859
    C_iterN(4,22) =  0.069796468424520488
    C_iterN(5,22) =  0.0859416062170677274
    C_iterN(6,22) =  0.100414144442880965
    C_iterN(7,22) =  0.112932296080539218
    C_iterN(8,22) =  0.123252376810512424
    C_iterN(9,22) =  0.13117350478706237
    C_iterN(10,22) =  0.136541498346015171
    C_iterN(11,22) =  0.139251872855631993
    C_iterN(12,22) =  0.139251872855631993
    C_iterN(13,22) =  0.136541498346015171
    C_iterN(14,22) =  0.131173504787062371
    C_iterN(15,22) =  0.123252376810512424
    C_iterN(16,22) =  0.112932296080539218
    C_iterN(17,22) =  0.100414144442880965
    C_iterN(18,22) =  0.085941606217067727
    C_iterN(19,22) =  0.06979646842452049
    C_iterN(20,22) =  0.052293335152683286
    C_iterN(21,22) =  0.033774901584814155
    C_iterN(22,22) =  0.014627995298272201
    X_iterN(1,22) = -0.994294585482399292
    X_iterN(2,22) = -0.970060497835428727
    X_iterN(3,22) = -0.926956772187174001
    X_iterN(4,22) = -0.865812577720300137
    X_iterN(5,22) = -0.787816805979208162
    X_iterN(6,22) = -0.69448726318668278
    X_iterN(7,22) = -0.587640403506911593
    X_iterN(8,22) = -0.469355837986757026
    X_iterN(9,22) = -0.341935820892084225
    X_iterN(10,22) = -0.207860426688221286
    X_iterN(11,22) = -0.069739273319722221
    X_iterN(12,22) =  0.069739273319722221
    X_iterN(13,22) =  0.207860426688221286
    X_iterN(14,22) =  0.341935820892084225
    X_iterN(15,22) =  0.469355837986757026
    X_iterN(16,22) =  0.58764040350691159
    X_iterN(17,22) =  0.69448726318668278
    X_iterN(18,22) =  0.787816805979208162
    X_iterN(19,22) =  0.865812577720300137
    X_iterN(20,22) =  0.926956772187174001
    X_iterN(21,22) =  0.970060497835428727
    X_iterN(22,22) =  0.994294585482399292

    !! NODE = 23
    C_iterN(1,23) =  0.013411859487141772
    C_iterN(2,23) =  0.030988005856979444
    C_iterN(3,23) =  0.0480376717310846686
    C_iterN(4,23) =  0.064232421408525852
    C_iterN(5,23) =  0.079281411776718955
    C_iterN(6,23) =  0.0929157660600351475
    C_iterN(7,23) =  0.10489209146454141
    C_iterN(8,23) =  0.114996640222411365
    C_iterN(9,23) =  0.123049084306729531
    C_iterN(10,23) =  0.12890572218808215
    C_iterN(11,23) =  0.132462039404696617
    C_iterN(12,23) =  0.13365457218610618
    C_iterN(13,23) =  0.13246203940469662
    C_iterN(14,23) =  0.12890572218808215
    C_iterN(15,23) =  0.12304908430672953
    C_iterN(16,23) =  0.114996640222411365
    C_iterN(17,23) =  0.10489209146454141
    C_iterN(18,23) =  0.09291576606003515
    C_iterN(19,23) =  0.079281411776718955
    C_iterN(20,23) =  0.064232421408525852
    C_iterN(21,23) =  0.048037671731084669
    C_iterN(22,23) =  0.030988005856979444
    C_iterN(23,23) =  0.013411859487141772
    X_iterN(1,23) = -0.994769334997552124
    X_iterN(2,23) = -0.972542471218115232
    X_iterN(3,23) = -0.932971086826016102
    X_iterN(4,23) = -0.876752358270441667
    X_iterN(5,23) = -0.804888401618839892
    X_iterN(6,23) = -0.718661363131950195
    X_iterN(7,23) = -0.619609875763646156
    X_iterN(8,23) = -0.50950147784600755
    X_iterN(9,23) = -0.390301038030290831
    X_iterN(10,23) = -0.264135680970344931
    X_iterN(11,23) = -0.133256824298466111
    X_iterN(12,23) = 0
    X_iterN(13,23) = 0.133256824298466111
    X_iterN(14,23) = 0.26413568097034493
    X_iterN(15,23) = 0.390301038030290831
    X_iterN(16,23) = 0.50950147784600755
    X_iterN(17,23) = 0.619609875763646156
    X_iterN(18,23) = 0.718661363131950195
    X_iterN(19,23) = 0.804888401618839892
    X_iterN(20,23) = 0.876752358270441667
    X_iterN(21,23) = 0.932971086826016102
    X_iterN(22,23) = 0.972542471218115232
    X_iterN(23,23) = 0.994769334997552124

    !! NODE = 24
    C_iterN(1,24) =  0.0123412297999872
    C_iterN(2,24) =  0.02853138862893366
    C_iterN(3,24) =  0.044277438817419806
    C_iterN(4,24) =  0.059298584915436781
    C_iterN(5,24) =  0.07334648141108031
    C_iterN(6,24) =  0.086190161531953276
    C_iterN(7,24) =  0.097618652104113888
    C_iterN(8,24) =  0.107444270115965635
    C_iterN(9,24) =  0.115505668053725601
    C_iterN(10,24) =  0.121670472927803391
    C_iterN(11,24) =  0.125837456346828296
    C_iterN(12,24) =  0.127938195346752157
    C_iterN(13,24) =  0.127938195346752157
    C_iterN(14,24) =  0.1258374563468283
    C_iterN(15,24) =  0.121670472927803391
    C_iterN(16,24) =  0.1155056680537256
    C_iterN(17,24) =  0.10744427011596563
    C_iterN(18,24) =  0.09761865210411389
    C_iterN(19,24) =  0.08619016153195328
    C_iterN(20,24) =  0.07334648141108031
    C_iterN(21,24) =  0.0592985849154367808
    C_iterN(22,24) =  0.04427743881741981
    C_iterN(23,24) =  0.028531388628933663
    C_iterN(24,24) =  0.0123412297999872
    X_iterN(1,24) = -0.99518721999702136
    X_iterN(2,24) = -0.974728555971309498
    X_iterN(3,24) = -0.938274552002732759
    X_iterN(4,24) = -0.886415527004401034
    X_iterN(5,24) = -0.820001985973902922
    X_iterN(6,24) = -0.740124191578554364
    X_iterN(7,24) = -0.648093651936975569
    X_iterN(8,24) = -0.545421471388839536
    X_iterN(9,24) = -0.433793507626045139
    X_iterN(10,24) = -0.31504267969616337
    X_iterN(11,24) = -0.19111886747361631
    X_iterN(12,24) = -0.064056892862605626
    X_iterN(13,24) =  0.064056892862605626
    X_iterN(14,24) =  0.191118867473616309
    X_iterN(15,24) =  0.315042679696163374
    X_iterN(16,24) =  0.433793507626045139
    X_iterN(17,24) =  0.545421471388839536
    X_iterN(18,24) =  0.648093651936975569
    X_iterN(19,24) =  0.740124191578554364
    X_iterN(20,24) =  0.820001985973902922
    X_iterN(21,24) =  0.886415527004401034
    X_iterN(22,24) =  0.938274552002732759
    X_iterN(23,24) =  0.974728555971309498
    X_iterN(24,24) =  0.99518721999702136

    !! NODE = 25
    C_iterN(1,25) =  0.011393798501026288
    C_iterN(2,25) =  0.026354986615032137
    C_iterN(3,25) =  0.04093915670130631
    C_iterN(4,25) =  0.054904695975835192
    C_iterN(5,25) =  0.068038333812356917
    C_iterN(6,25) =  0.080140700335001018
    C_iterN(7,25) =  0.09102826198296365
    C_iterN(8,25) =  0.100535949067050644
    C_iterN(9,25) =  0.108519624474263653
    C_iterN(10,25) =  0.114858259145711648
    C_iterN(11,25) =  0.11945576353578477
    C_iterN(12,25) =  0.12224244299031004
    C_iterN(13,25) =  0.123176053726715451
    C_iterN(14,25) =  0.122242442990310042
    C_iterN(15,25) =  0.11945576353578477
    C_iterN(16,25) =  0.114858259145711648
    C_iterN(17,25) =  0.10851962447426365
    C_iterN(18,25) =  0.10053594906705064
    C_iterN(19,25) =  0.09102826198296365
    C_iterN(20,25) =  0.08014070033500102
    C_iterN(21,25) =  0.068038333812356917
    C_iterN(22,25) =  0.054904695975835192
    C_iterN(23,25) =  0.040939156701306313
    C_iterN(24,25) =  0.026354986615032137
    C_iterN(25,25) =  0.011393798501026288
    X_iterN(1,25) = -0.995556969790498098
    X_iterN(2,25) = -0.976663921459517512
    X_iterN(3,25) = -0.942974571228974339
    X_iterN(4,25) = -0.894991997878275369
    X_iterN(5,25) = -0.833442628760834001
    X_iterN(6,25) = -0.759259263037357631
    X_iterN(7,25) = -0.673566368473468365
    X_iterN(8,25) = -0.57766293024122297
    X_iterN(9,25) = -0.473002731445714961
    X_iterN(10,25) = -0.361172305809387838
    X_iterN(11,25) = -0.243866883720988432
    X_iterN(12,25) = -0.122864692610710396
    X_iterN(13,25) =  0
    X_iterN(14,25) =  0.122864692610710396
    X_iterN(15,25) =  0.243866883720988432
    X_iterN(16,25) =  0.361172305809387838
    X_iterN(17,25) =  0.47300273144571496
    X_iterN(18,25) =  0.577662930241222968
    X_iterN(19,25) =  0.673566368473468365
    X_iterN(20,25) =  0.759259263037357631
    X_iterN(21,25) =  0.833442628760834001
    X_iterN(22,25) =  0.894991997878275369
    X_iterN(23,25) =  0.942974571228974339
    X_iterN(24,25) =  0.976663921459517512
    X_iterN(25,25) =  0.995556969790498098

    END SUBROUTINE Gauss_Quad


! ---------------------------------------------------------------------
CLOSE(10)

CLOSE(3)
END PROGRAM CST_MULTI
