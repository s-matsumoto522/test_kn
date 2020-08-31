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
!********************************
!   各格子点における速度の計算  *
!********************************
    subroutine set_velocity(Vx, Vy, Vz, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    Vx(iX,iY,iZ) = sin((X(iX,iY,iZ) + 0.5d0*dX) + (Y(iX,iY,iZ) + 0.5d0*dY) + (Z(iX,iY,iZ) + 0.5d0*dZ))
                    Vy(iX,iY,iZ) = cos((X(iX,iY,iZ) + 0.5d0*dX) + (Y(iX,iY,iZ) + 0.5d0*dY) + (Z(iX,iY,iZ) + 0.5d0*dZ))
                    Vz(iX,iY,iZ) = sin(2.0d0*((X(iX,iY,iZ) + 0.5d0*dX) + (Y(iX,iY,iZ) + 0.5d0*dY) + (Z(iX,iY,iZ) + 0.5d0*dZ)))
                enddo
            enddo
        enddo
    end subroutine set_velocity
!********************************************
!   移流項のx成分を計算するサブルーチン     *
!********************************************
    subroutine cal_x_advection(Vx, Vy, Vz, Ax)
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Ax(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        double precision Cx, Cy, Cz  !計算用数値
        !---計算用数値の設定---
        Cx = ddX / 4.0d0
        Cy = ddY / 4.0d0
        Cz = ddZ / 4.0d0
        !---移流項の計算---
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Ax(iX, iY, iZ) = Cx*((Vx(iX-1,iY,iZ) + Vx(iX,iY,iZ))**2 - (Vx(iX,iY,iZ) + Vx(iX+1,iY,iZ))**2) &
                                    + Cy*((Vy(iX,iY-1,iZ) + Vy(iX+1,iY-1,iZ))*(Vx(iX,iY-1,iZ) + Vx(iX,iY,iZ)) &
                                    - (Vy(iX,iY,iZ) + Vy(iX+1,iY,iZ))*(Vx(iX,iY,iZ) + Vx(iX,iY+1,iZ))) &
                                    + Cz*((Vz(iX,iY,iZ-1) + Vz(iX+1,iY,iZ-1))*(Vx(iX,iY,iZ-1) + Vx(iX,iY,iZ)) &
                                    - (Vz(iX,iY,iZ) + Vz(iX+1,iY,iZ))*(Vx(iX,iY,iZ) + Vx(iX,iY,iZ+1)))
                enddo
            enddo
        enddo
    end subroutine cal_x_advection
!********************************************
!   移流項のy成分を計算するサブルーチン     *
!********************************************
    subroutine cal_y_advection(Vx, Vy, Vz, Ay)
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Ay(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        double precision Cx, Cy, Cz  !計算用数値
        !---計算用数値の設定---
        Cx = ddX / 4.0d0
        Cy = ddY / 4.0d0
        Cz = ddZ / 4.0d0
        !---移流項の計算---
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Ay(iX, iY, iZ) = Cx*((Vx(iX-1,iY,iZ) + Vx(iX-1,iY+1,iZ))*(Vy(iX-1,iY,iZ) + Vy(iX,iY,iZ)) &
                                    - (Vx(iX,iY,iZ) + Vx(iX,iY+1,iZ))*(Vy(iX,iY,iZ) + Vy(iX+1,iY,iZ))) &
                                    + Cy*((Vy(iX,iY-1,iZ) + Vy(iX,iY,iZ))**2 - (Vy(iX,iY,iZ) + Vy(iX,iY+1,iZ))**2) &
                                    + Cz*((Vz(iX,iY,iZ-1) + Vz(iX,iY+1,iZ-1))*(Vy(iX,iY,iZ-1) + Vy(iX,iY,iZ)) &
                                    - (Vz(iX,iY,iZ) + Vz(iX,iY+1,iZ))*(Vy(iX,iY,iZ) + Vy(iX,iY,iZ+1)))
                enddo
            enddo
        enddo
    end subroutine cal_y_advection
!********************************************
!   移流項のz成分を計算するサブルーチン     *
!********************************************
    subroutine cal_z_advection(Vx, Vy, Vz, Az)
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Az(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        double precision Cx, Cy, Cz  !計算用数値
        !---計算用数値の設定---
        Cx = ddX / 4.0d0
        Cy = ddY / 4.0d0
        Cz = ddZ / 4.0d0
        !---移流項の計算---
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Az(iX, iY, iZ) = Cx*((Vx(iX-1,iY,iZ) + Vx(iX-1,iY,iZ+1))*(Vz(iX-1,iY,iZ) + Vz(iX,iY,iZ)) &
                                    - (Vx(iX,iY,iZ) + Vx(iX,iY,iZ+1))*(Vz(iX,iY,iZ) + Vz(iX+1,iY,iZ))) &
                                    + Cy*((Vy(iX,iY-1,iZ) + Vy(iX,iY-1,iZ+1))*(Vz(iX,iY-1,iZ) + Vz(iX,iY,iZ)) &
                                    - (Vy(iX,iY,iZ) + Vy(iX,iY,iZ+1))*(Vz(iX,iY,iZ) + Vz(iX,iY+1,iZ))) &
                                    + Cz*((Vz(iX,iY,iZ-1) + Vz(iX,iY,iZ))**2 - (Vz(iX,iY,iZ) + Vz(iX,iY,iZ+1))**2)
                enddo
            enddo
        enddo
    end subroutine cal_z_advection
!******************************************
!   粘性項のx成分を計算するサブルーチン   *
!******************************************
    subroutine cal_x_viscosity(Vx, Bx)
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Bx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        !---粘性項の計算---
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Bx(iX, iY, iZ) = ddX2*(Vx(iX-1,iY,iZ) - 2.0d0*Vx(iX,iY,iZ) + Vx(iX+1,iY,iZ)) &
                                    + ddY2*(Vx(iX,iY-1,iZ) - 2.0d0*Vx(iX,iY,iZ) + Vx(iX,iY+1,iZ)) &
                                    + ddZ2*(Vx(iX,iY,iZ-1) - 2.0d0*Vx(iX,iY,iZ) + Vx(iX,iY,iZ+1))
                enddo
            enddo
        enddo
    end subroutine cal_x_viscosity
!******************************************
!   粘性項のy成分を計算するサブルーチン   *
!******************************************
    subroutine cal_y_viscosity(Vy, By)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: By(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        !---粘性項の計算---
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    By(iX, iY, iZ) = ddX2*(Vy(iX-1,iY,iZ) - 2.0d0*Vy(iX,iY,iZ) + Vy(iX+1,iY,iZ)) &
                                    + ddY2*(Vy(iX,iY-1,iZ) - 2.0d0*Vy(iX,iY,iZ) + Vy(iX,iY+1,iZ)) &
                                    + ddZ2*(Vy(iX,iY,iZ-1) - 2.0d0*Vy(iX,iY,iZ) + Vy(iX,iY,iZ+1))
                enddo
            enddo
        enddo
    end subroutine cal_y_viscosity
!******************************************
!   粘性項のz成分を計算するサブルーチン   *
!******************************************
    subroutine cal_z_viscosity(Vz, Bz)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: Bz(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        !---粘性項の計算---
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Bz(iX, iY, iZ) = ddX2*(Vz(iX-1,iY,iZ) - 2.0d0*Vz(iX,iY,iZ) + Vz(iX+1,iY,iZ)) &
                                    + ddY2*(Vz(iX,iY-1,iZ) - 2.0d0*Vz(iX,iY,iZ) + Vz(iX,iY+1,iZ)) &
                                    + ddZ2*(Vz(iX,iY,iZ-1) - 2.0d0*Vz(iX,iY,iZ) + Vz(iX,iY,iZ+1))
                enddo
            enddo
        enddo
    end subroutine cal_z_viscosity
!****************************************
!   圧力項をx方向に微分するサブルーチン *
!****************************************
    subroutine dif_pressure_x(p, dpx)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: dpx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    dpx(iX, iY, iZ) = ddX*(-p(iX, iY, iZ)+p(iX+1, iY, iZ))
                enddo
            enddo
        enddo
    end subroutine dif_pressure_x
!****************************************
!   圧力項をy方向に微分するサブルーチン *
!****************************************
    subroutine dif_pressure_y(p, dpy)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: dpy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    dpy(iX, iY, iZ) = ddY*(-p(iX, iY, iZ)+p(iX, iY+1, iZ))
                enddo
            enddo
        enddo
    end subroutine dif_pressure_y
!****************************************
!   圧力項をz方向に微分するサブルーチン *
!****************************************
    subroutine dif_pressure_z(p, dpz)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: dpz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    dpz(iX, iY, iZ) = ddZ*(-p(iX, iY, iZ)+p(iX, iY, iZ+1))
                enddo
            enddo
        enddo
    end subroutine dif_pressure_z
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
    subroutine cal_RHS(Vx, Vy, Vz, p, Ax2, istep, RHS)
        integer, intent(in) :: istep
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(inout) :: Ax2(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: RHS(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
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
                        RHS(iX,iY,iZ) = Vx(iX,iY,iZ)+C*Bx(iX,iY,iZ)-dt*dpx(iX,iY,iZ)+dt*Ax2(iX,iY,iZ)
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
    end subroutine cal_RHS
!*******************************
!   反復計算を行うサブルーチン *
!*******************************
    subroutine cal_itr(RHS, Vx_p, X, Y, Z)
        double precision, intent(in) :: RHS(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
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
        double precision error, norm_rhs, norm_er
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
            norm_rhs = 0.0d0
            norm_er = 0.0d0
            !境界を除く領域で計算
            do iZ = NZmin, NZmax-1
                do iY = NYmin, NYmax-1
                    do start = NXmin, NXmin+1
                        do iX = start, NXmax-1, 2
                            Er(iX,iY,iZ) = RHS(iX,iY,iZ)-Cx*(Vx_p(iX-1,iY,iZ)+Vx(iX+1,iY,iZ)) &
                                        -Cy*(Vx_p(iX,iY-1,iZ)+Vx_p(iX,iY+1,iZ)) &
                                        -Cz*(Vx_p(iX,iY,iZ-1)+Vx_p(iX,iY,iZ+1))-C0*Vx_p(iX,iY,iZ)
                            !---解の更新---
                            Vx_p(iX,iY,iZ) = Vx_p(iX,iY,iZ)+beta*dC0*Er(iX,iY,iZ)
                            norm_er = norm_er + Er(iX,iY,iZ)**2
                            norm_rhs = norm_rhs + RHS(iX,iY,iZ)**2
                        enddo
                    enddo
                enddo
            enddo
            norm_er = sqrt(norm_er / dble(Ng))
            norm_rhs = aqrt(norm_rhs / dble(Ng))
            !誤差の計算
            error = norm_er / norm_rhs
            !境界条件の設定
            call set_xpvel_bc(Vx_p, X, Y, Z)
            !収束判定
            if(error < eps) then
                write(*, *) 'converged.'
                exit
            endif
        enddo
    end subroutine cal_itr
!**************************************
!   予測速度を求めるサブルーチン(KN)  *
!**************************************
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
        double precision RHS(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        !---SOR法で連立方程式を解く---
        Vx_p(:, :, :) = 0.0d0       !初期値の設定
        !---連立方程式の右辺の計算---
        call cal_RHS(Vx, Vy, Vz, p, Ax2, 1, RHS)
        !---反復計算---
        call cal_itr(RHS, Vx_p, X, Y, Z)
    end subroutine prediction_x_kn