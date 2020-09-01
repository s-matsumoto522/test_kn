module prediction_x
    implicit none
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = -80, NXmax = 80     !x方向の計算領域の形状
    integer, parameter :: NYmin = -80, NYmax = 80     !y方向の計算領域の形状
    integer, parameter :: NZmin = -80, NZmax = 80     !z方向の計算領域の形状
    double precision, parameter :: Xmax = 0.5d0*pi, Ymax = 0.5d0*pi, Zmax = 0.5d0*pi  !各方向の計算領域の最大値
    double precision, parameter :: dt = 1.0d-2      !時間刻み幅
    double precision, parameter :: Re = 10d0        !レイノルズ数
    integer, save :: Ng                             !格子点数
    double precision, save :: dX, dY, dZ            !各方向の刻み幅
    double precision, save :: ddX, ddY, ddZ    !各方向の刻み幅の逆数
    double precision, save :: ddX2, ddY2, ddZ2 !各方向の刻み幅の逆数の二乗
    double precision, save :: Xmin, Ymin, Zmin      !各方向の計算領域の最小値
contains
!********************************
!   格子を設定するサブルーチン  *
!********************************
subroutine set_grid(X, Y, Z)
    double precision, intent(out) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    integer iX, iY, iZ
    !---各方向の刻み幅を計算---
    dX = Xmax / dble(NXmax)
    dY = Ymax / dble(NYmax)
    dZ = Zmax / dble(NZmax)
    !---計算用数値の設定---
    ddX = 1.0d0 / dX
    ddY = 1.0d0 / dY
    ddZ = 1.0d0 / dZ
    ddX2 = 1.0d0 / (dX**2)
    ddY2 = 1.0d0 / (dY**2)
    ddZ2 = 1.0d0 / (dZ**2)
    Ng = (NXmax-NXmin+1)*(NYmax-NYmin+1)*(NZmax-NZmin+1)
    !---各方向の計算領域の最小値を計算---
    Xmin = dX*dble(NXmin)
    Ymin = dY*dble(NYmin)
    Zmin = dZ*dble(NZmin)
    !---格子点の設定---
    do iZ = NZmin-1, NZmax+1
        do iY = NYmin-1, NYmax+1
            do iX = NXmin-1, NXmax+1
                X(iX, iY, iZ) = dX*dble(iX)
                Y(iX, iY, iZ) = dY*dble(iY)
                Z(iX, iY, iZ) = dZ*dble(iZ)
            enddo
        enddo
    enddo
end subroutine set_grid
!*******************************************************
!   予測速度のx成分の境界条件を設定するサブルーチン    *
!*******************************************************
    subroutine set_xpvel_bc(Vx_p, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision C
        C = 0.5d0*dt/Re
    end subroutine set_xpvel_bc
!********************************************************
!   KN法において連立方程式の右辺を計算するサブルーチン  *
!********************************************************
    subroutine cal_RHSx(Vx, Vy, Vz, p, Ax2, istep, RHSx)
        integer, intent(in) :: istep
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(inout) :: Ax2(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: RHSx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision Ax1(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision Bx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision dpx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision C
        integer iX, iY, iZ
        !---計算用数値の計算---
        C = 0.5d0*dt/Re
        if(istep == 1) then
            call cal_x_advection(Vx, Vy, Vz, Ax2)
            call cal_x_viscosity(Vx, Bx)
            call dif_pressure_x(p, dpx)
            do iZ = NZmin, NZmax-1
                do iY = NYmin, NYmax-1
                    do iX = NXmin, NXmax-1
                        RHS(iX,iY,iZ) = sin((X(iX,iY,iZ)+0.5d0*dX)+Y(iX,iY,iZ)+Z(iX,iY,iZ))
                    enddo
                enddo
            enddo
        else
            Ax1(:, :, :) = Ax2(:, :, :)
            call cal_x_advection(Vx, Vy, Vz, Ax2)
            call cal_x_viscosity(Vx, Bx)
            call dif_pressure_x(p, dpx)
            do iZ = NZmin, NZmax-1
                do iY = NYmin, NYmax-1
                    do iX = NXmin, NXmax-1
                        RHS(iX,iY,iZ) = Vx(iX,iY,iZ)+C*Bx(iX,iY,iZ)-dt*dpx(iX,iY,iZ) &
                                        +0.5d0*dt*(3.0d0*Ax2(iX,iY,iZ)-Ax1(iX,iY,iZ))
                    enddo
                enddo
            enddo
        endif
    end subroutine cal_RHSx
!*******************************
!   反復計算を行うサブルーチン *
!*******************************
    subroutine cal_itr_x(RHSx, Vx_p, X, Y, Z)
        double precision, intent(in) :: RHSx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision Er(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision C, C0, Cx, Cy, Cz, dC0     !計算用数値
        integer, parameter :: itrmax = 10000        !最大反復回数
        double precision, parameter :: eqs = 1.0d-4 !誤差のしきい値
        double precision, parameter :: beta = 0.8d0 !過緩和係数
        integer itr, iX, iY, iZ, start
        double precision error, norm_RHSx, norm_er
        !---計算用数値の計算---
        C = 0.5d0*dt/Re
        C0 = (1.0d0 - 2.0d0*C*(ddX2 + ddY2 + ddZ2))
        dC0 = 1.0d0/C0
        Cx = C*ddX2
        Cy = C*ddY2
        Cz = C*ddZ2
        do itr = 1, itrmax
            !誤差の初期値の設定
            error = 0.0d0
            norm_RHSx = 0.0d0
            norm_er = 0.0d0
            !境界を除く領域で計算
            do iZ = NZmin, NZmax-1
                do iY = NYmin, NYmax-1
                    do start = NXmin, NXmin+1
                        do iX = start, NXmax-1, 2
                            Er(iX,iY,iZ) = RHSx(iX,iY,iZ)-Cx*(Vx_p(iX-1,iY,iZ)+Vx(iX+1,iY,iZ)) &
                                        -Cy*(Vx_p(iX,iY-1,iZ)+Vx_p(iX,iY+1,iZ)) &
                                        -Cz*(Vx_p(iX,iY,iZ-1)+Vx_p(iX,iY,iZ+1))-C0*Vx_p(iX,iY,iZ)
                            !---解の更新---
                            Vx_p(iX,iY,iZ) = Vx_p(iX,iY,iZ)+beta*dC0*Er(iX,iY,iZ)
                            norm_er = norm_er + Er(iX,iY,iZ)**2
                            norm_RHSx = norm_RHSx + RHSx(iX,iY,iZ)**2
                        enddo
                    enddo
                enddo
            enddo
            norm_er = sqrt(norm_er / dble(Ng))
            norm_RHSx = aqrt(norm_RHSx / dble(Ng))
            !誤差の計算
            error = norm_er / norm_RHSx
            !境界条件の設定
            call set_xpvel_bc(Vx_p, X, Y, Z)
            !収束判定
            if(error < eps) then
                write(*, *) 'converged.'
                exit
            endif
        enddo
    end subroutine cal_itr_x
!*********************************************
!   x方向の予測速度を求めるサブルーチン(KN)  *
!*********************************************
    subroutine prediction_x_kn(Vx, Vy, Vz, p, Ax2, Vx_p, istep, X, Y, Z)
        integer, intent(in) :: istep
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Ax2(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Vx_p(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision RHSx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        !---SOR法で連立方程式を解く---
        Vx_p(:, :, :) = 0.0d0       !初期値の設定
        !---連立方程式の右辺の計算---
        call cal_RHSx(Vx, Vy, Vz, p, Ax2, 1, RHSx)
        !---反復計算---
        call cal_itr_x(RHSx, Vx_p, X, Y, Z)
    end subroutine prediction_x_kn
!************************************
!   解析解を計算するサブルーチン    *
!************************************
    subroutine cal_error(RHS, RHS_th, error)
        double precision, intent(in) :: RHS(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: RHS_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: error
    end subroutine cal_error
!********************************
!   誤差を計算するサブルーチン  *
!********************************