module x_RHSy
    implicit none
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = 1, NXmax = 160     !x方向の計算領域の形状
    integer, parameter :: NYmin = 1, NYmax = 160     !y方向の計算領域の形状
    integer, parameter :: NZmin = 1, NZmax = 160     !z方向の計算領域の形状
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax = 2.0d0*pi, Zmax = 2.0d0*pi  !各方向の計算領域の最大値
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
                Vx(iX, iY, iZ) = sin((X(iX,iY,iZ) + 0.5d0*dX) + Y(iX,iY,iZ)  + Z(iX,iY,iZ))
                Vy(iX, iY, iZ) = cos(X(iX,iY,iZ) + (Y(iX,iY,iZ) + 0.5d0*dY) + Z(iX,iY,iZ))
                Vz(iX, iY, iZ) = sin(2.0d0*(X(iX,iY,iZ) + Y(iX,iY,iZ) + (Z(iX,iY,iZ) + 0.5d0*dZ)))
            enddo
        enddo
    enddo
end subroutine set_velocity
!********************************
!   各格子点における圧力の計算  *
!********************************
    subroutine set_pressure(p, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax+1
            do iY = NYmin-1, NZmax+1
                do iX = NXmin-1, NZmax+1
                    p(iX,iY,iZ) = cos(2.0d0*(X(iX,iY,iZ) + Y(iX,iY,iZ) + Z(iX,iY,iZ)))
                enddo
            enddo
        enddo
    end subroutine set_pressure
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
                    dpx(iX, iY, iZ) = ddX*(-p(iX, iY, iZ) + p(iX+1, iY, iZ))
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
                    dpy(iX, iY, iZ) = ddY*(-p(iX, iY, iZ) + p(iX, iY+1, iZ))
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
                    dpz(iX, iY, iZ) = ddZ*(-p(iX, iY, iZ) + p(iX, iY, iZ+1))
                enddo
            enddo
        enddo
    end subroutine dif_pressure_z
!***************************************************************
!   KN法において連立方程式の右辺(y成分)を計算するサブルーチン  *
!***************************************************************
    subroutine cal_RHSy(Vx, Vy, Vz, p, Ay2, istep, RHSy)
        integer, intent(in) :: istep
        double precision, intent(in) :: Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(inout) :: Ay2(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: RHSy(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision Ay1(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision By(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision dpy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision C
        integer iX, iY, iZ
        !---計算用数値の計算---
        C = 0.5d0*dt/Re
        if(istep == 1) then
            call cal_x_advection(Vx, Vy, Vz, Ay2)
            call cal_x_viscosity(Vx, By)
            call dif_pressure_x(p, dpy)
            do iZ = NZmin, NZmax-1
                do iY = NYmin, NYmax-1
                    do iX = NXmin, NXmax-1
                        RHSy(iX,iY,iZ) = Vx(iX,iY,iZ)+C*By(iX,iY,iZ)-dt*dpy(iX,iY,iZ)+dt*Ay2(iX,iY,iZ)
                    enddo
                enddo
            enddo
        else
            Ay1(:, :, :) = Ay2(:, :, :)
            call cal_x_advection(Vx, Vy, Vz, Ay2)
            call cal_x_viscosity(Vx, By)
            call dif_pressure_x(p, dpy)
            do iZ = NZmin, NZmax-1
                do iY = NYmin, NYmax-1
                    do iX = NXmin, NXmax-1
                        RHSy(iX,iY,iZ) = Vx(iX,iY,iZ)+C*By(iX,iY,iZ)-dt*dpy(iX,iY,iZ) &
                                        +0.5d0*dt*(3.0d0*Ay2(iX,iY,iZ)-Ay1(iX,iY,iZ))
                    enddo
                enddo
            enddo
        endif
    end subroutine cal_RHSy
!****************************************
!   右辺の解析解を計算するサブルーチン  *
!****************************************
    subroutine cal_RHSy_th(X, Y, Z, RHSy_th)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: RHSy_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        double precision C
        !---計算用数値の計算---
        C = 0.5d0*dt/Re
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    RHSy_th(iX,iY,iZ) = sin(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ)) &
                            -3.0d0*C*sin(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY) &
                            +Z(iX,iY,iZ))+2.0d0*dt*sin(2.0d0*(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ))) &
                            +dt*(-2.0d0*(sin(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ))*cos(X(iX,iY,iZ) &
                            +(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ)))-(cos(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY) &
                            +Z(iX,iY,iZ)))**2+(sin(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ)))**2 &
                            -cos(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ))*sin(2.0d0*(X(iX,iY,iZ) &
                            +(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ)))-2.0d0*sin(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY) &
                            +Z(iX,iY,iZ))*cos(2.0d0*(X(iX,iY,iZ)+(Y(iX,iY,iZ)+0.5d0*dY)+Z(iX,iY,iZ))))
                enddo
            enddo
        enddo
    end subroutine cal_RHSy_th
end module x_RHSy

program main
    use x_RHSy
    implicit none
    double precision X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Ay2(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision RHSy(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision RHSy_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision error, error_norm
    integer iX, iY, iZ
    call set_grid(X, Y, Z)
    call set_velocity(Vx, Vy, Vz, X, Y, Z)
    call set_pressure(p, X, Y, Z)
    call cal_RHSy(Vx, Vy, Vz, p, Ay2, 1, RHSy)
    call cal_RHSy_th(X, Y, Z, RHSy_th)
    error_norm = 0.0d0
    do iZ = NZmin, NZmax-1
        do iY = NYmin, NYmax-1
            do iX = NXmin, NXmax-1
                error_norm = error_norm + (RHSy(iX,iY,iZ) - RHSy_th(iX,iY,iZ))**2
            enddo
        enddo
    enddo
    error = sqrt(error_norm / dble(Ng))
    open(11, file = 'chk_RHSy_precision.dat', position = 'append')
    write(11, *) dY, error
end program main