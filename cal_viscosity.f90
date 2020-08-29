module x_vis
    implicit none
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = -80, NXmax = 80     !x方向の計算領域の形状
    integer, parameter :: NYmin = -80, NYmax = 80     !y方向の計算領域の形状
    integer, parameter :: NZmin = -80, NZmax = 80     !z方向の計算領域の形状
    double precision, parameter :: Xmax = pi, Ymax = pi, Zmax = pi  !各方向の計算領域の最大値
    double precision, parameter :: NU = 1.0d0       !動粘性係数
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
                    Vx(iX, iY, iZ) = sin((X(iX,iY,iZ) + 0.5d0*dX) + (Y(iX,iY,iZ) + 0.5d0*dY) + (Z(iX,iY,iZ) + 0.5d0*dZ))
                    Vy(iX, iY, iZ) = cos((X(iX,iY,iZ) + 0.5d0*dX) + (Y(iX,iY,iZ) + 0.5d0*dY) + (Z(iX,iY,iZ) + 0.5d0*dZ))
                    Vz(iX, iY, iZ) = sin(2.0d0*((X(iX,iY,iZ) + 0.5d0*dX) + (Y(iX,iY,iZ) + 0.5d0*dY) + (Z(iX,iY,iZ) + 0.5d0*dZ)))
                enddo
            enddo
        enddo
    end subroutine set_velocity
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
                    Bx(iX, iY, iZ) = NU*(ddX2*(Vx(iX-1,iY,iZ) - 2.0d0*Vx(iX,iY,iZ) + Vx(iX+1,iY,iZ)) &
                                    + ddY2*(Vx(iX,iY-1,iZ) - 2.0d0*Vx(iX,iY,iZ) + Vx(iX,iY+1,iZ)) &
                                    + ddZ2*(Vx(iX,iY,iZ-1) - 2.0d0*Vx(iX,iY,iZ) + Vx(iX,iY,iZ+1)))
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
                    By(iX, iY, iZ) = NU*(ddX2*(Vy(iX-1,iY,iZ) - 2.0d0*Vy(iX,iY,iZ) + Vy(iX+1,iY,iZ)) &
                                    + ddY2*(Vy(iX,iY-1,iZ) - 2.0d0*Vy(iX,iY,iZ) + Vy(iX,iY+1,iZ)) &
                                    + ddZ2*(Vy(iX,iY,iZ-1) - 2.0d0*Vy(iX,iY,iZ) + Vy(iX,iY,iZ+1)))
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
                    Bz(iX, iY, iZ) = NU*(ddX2*(Vz(iX-1,iY,iZ) - 2.0d0*Vz(iX,iY,iZ) + Vz(iX+1,iY,iZ)) &
                    + ddY2*(Vz(iX,iY-1,iZ) - 2.0d0*Vz(iX,iY,iZ) + Vz(iX,iY+1,iZ)) &
                    + ddZ2*(Vz(iX,iY,iZ-1) - 2.0d0*Vz(iX,iY,iZ) + Vz(iX,iY,iZ+1)))
                enddo
            enddo
        enddo
    end subroutine cal_z_viscosity
!********************************************
!   粘性項の解析解を計算するサブルーチン    *
!********************************************
    subroutine cal_theory(Bx_th, By_th, Bz_th, X, Y, Z)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: Bx_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: By_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: Bz_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        integer iX, iY, iZ
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Bx_th(iX, iY, iZ) = -3.0d0*NU*sin(X(iX,iY,iZ) + Y(iX,iY,iZ) + Z(iX,iY,iZ))
                    By_th(iX, iY, iZ) = -3.0d0*NU*cos(X(iX,iY,iZ) + Y(iX,iY,iZ) + Z(iX,iY,iZ))
                    Bz_th(iX, iY, iZ) = -12.0d0*NU*sin(2.0d0*(X(iX,iY,iZ) + Y(iX,iY,iZ) + Z(iX,iY,iZ)))
                enddo
            enddo
        enddo
    end subroutine cal_theory
!********************************
!   誤差を計算するサブルーチン  *
!********************************
    subroutine cal_error(errorx, errory, errorz, Bx, By, Bz, Bx_th, By_th, Bz_th)
        double precision, intent(in) :: Bx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: By(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: Bz(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: Bx_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: By_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(in) :: Bz_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision, intent(out) :: errorx, errory, errorz
        integer iX, iY, iZ, NX, NY, NZ
        double precision Erx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision Ery(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        double precision Erz(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
        errorx = 0.0d0
        errory = 0.0d0
        errorz = 0.0d0
        NX = NXmax - NXmin + 1
        NY = NYmax - NYmin + 1
        NZ = NZmax - NZmin + 1
        do iZ = NZmin, NZmax-1
            do iY = NYmin, NYmax-1
                do iX = NXmin, NXmax-1
                    Erx(iX, iY, iZ) = Bx(iX, iY, iZ) - Bx_th(iX, iY, iZ)
                    Ery(iX, iY, iZ) = By(iX, iY, iZ) - By_th(iX, iY, iZ)
                    Erz(iX, iY, iZ) = Bz(iX, iY, iZ) - Bz_th(iX, iY, iZ)
                    errorx = errorx + Erx(iX, iY, iZ)**2
                    errory = errory + Ery(iX, iY, iZ)**2
                    errorz = errorz + Erz(iX, iY, iZ)**2
                enddo
            enddo
        enddo
        errorx = errorx / (NX*NY*NZ)
        errory = errory / (NX*NY*NZ)
        errorz = errorz / (NX*NY*NZ)
    end subroutine cal_error
end module x_vis

program main
    use x_vis
    implicit none
    double precision X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Vx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Vz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision Bx(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision By(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision Bz(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision Bx_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision By_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision Bz_th(NXmin:NXmax-1, NYmin:NYmax-1, NZmin:NZmax-1)
    double precision errorx, errory, errorz
    call set_grid(X, Y, Z)
    call set_velocity(Vx, Vy, Vz, X, Y, Z)
    call cal_x_viscosity(Vx,Bx)
    call cal_y_viscosity(Vy,By)
    call cal_z_viscosity(Vz,Bz)
    call cal_theory(Bx_th, By_th, Bz_th, X, Y, Z)
    call cal_error(errorx, errory, errorz, Bx, By, Bz, Bx_th, By_th, Bz_th)
    open(11, file = 'chk_viscosity_precision.dat', position = 'append')
    write(11, *) dX, dY, dZ, errorx, errory, errorz
    close(11)
end program main
