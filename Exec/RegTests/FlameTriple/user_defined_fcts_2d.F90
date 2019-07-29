#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_2d_module

  use fuego_chemistry

implicit none

  private
  
  public :: bcfunction, zero_visc, set_Y_from_Phi

contains
  


!-----------------------

  subroutine bcfunction(x,y,dir,norm,time,u,v,rho,Yl,T,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use network,   only: nspec
      use mod_Fvar_def, only : dim
      use mod_Fvar_def, only : dv_control, tbase_control, V_in, f_flag_active_control
      use probdata_module, only : bcinit, rho_bc, Y_bc, T_bc, h_bc, v_bc, midtanh, widthtanh
      
      implicit none

      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(dim)
      integer dir, norm  ! This specify the direction and orientation of the face
      logical getuv
      REAL_T :: tanhval, vbase

      integer n

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      if ((dir == 2).and.(norm == 1)) then
        tanhval = 0.5d0*(1.0d0+TANH((x-midtanh)/widthtanh)) 
        rho = rho_bc(1,2) + tanhval*(rho_bc(1,1)-rho_bc(1,2)) 
        do n = 0, Nspec-1
          Yl(n) = Y_bc(n,2) + tanhval*(Y_bc(n,1) - Y_bc(n,2))
        end do
        T = T_bc(1)
        h = h_bc(1,2) + tanhval*(h_bc(1,1) - h_bc(1,2))
         
        if (getuv .eqv. .TRUE.) then
            
          u = zero
          if (f_flag_active_control == 1) then               
            vbase =  V_in + (time-tbase_control)*dV_control
            v = vbase + tanhval * (vbase*rho_bc(1,2)/rho_bc(1,1) - vbase)
          else 
            v = V_in + tanhval * (V_in*rho_bc(1,2)/rho_bc(1,1) - V_in)
          endif
        endif
      endif  

  end subroutine bcfunction

! ::: -----------------------------------------------------------
! ::: This routine will zero out diffusivity on portions of the
! ::: boundary that are inflow, allowing that a "wall" block
! ::: the complement aperture
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: diff      <=> diffusivity on edges
! ::: DIMS(diff) => index extent of diff array
! ::: lo,hi      => region of interest, edge-based
! ::: domlo,hi   => index extent of problem domain, edge-based
! ::: dx         => cell spacing
! ::: problo     => phys loc of lower left corner of prob domain
! ::: bc         => boundary condition flag (on orient)
! :::                   in BC_TYPES::physicalBndryTypes
! ::: idir       => which face, 0=x, 1=y
! ::: isrz       => 1 if problem is r-z
! ::: id         => index of state, 0=u
! ::: ncomp      => components to modify
! ::: 
! ::: -----------------------------------------------------------

  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
                           dx,problo,bc,idir,isrz,id,ncomp) &
                           bind(C, name="zero_visc")   

      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, LastSpec
      use mod_Fvar_def, only : domnhi, domnlo, dim
      
      implicit none
      integer DIMDEC(diff)
      integer lo(dim), hi(dim)
      integer domlo(dim), domhi(dim)
      integer bc(2*dim)
      integer idir, isrz, id, ncomp
      REAL_T  diff(DIMV(diff),*)
      REAL_T  dx(dim)
      REAL_T  problo(dim)

! Routine compiled but should be set by the user
! if there is a mix of inflox/wall at a boundary

  end subroutine zero_visc

  subroutine set_Y_from_Phi(phi,Yt)bind(C, name="set_Y_from_Phi")
  
      use mod_Fvar_def, only : maxspnml,fuelID, oxidID, bathID
      use network,   only: nspec
      use PeleLM_F,  only: pphys_get_spec_name2

      implicit none

      REAL_T, INTENT(IN)  :: phi
      REAL_T, INTENT(OUT) :: Yt(nspec)

      REAL_T :: a
      REAL_T :: Xt(nspec)
      INTEGER ::  n
      CHARACTER(LEN=maxspnml) :: name

      Xt(:) = zero
      
!     Set "a" for computing X from phi
!     hc + a.O2 -> b.CO2 + c.H2O
      
      call pphys_get_spec_name2(name,fuelID)

      a = 0.d0
      if (name .eq. 'CH4') then
         a = 2.0d0
      else if (name .eq. 'H2') then
         a = .5d0
      else if (name .eq. 'C3H8') then
         a = 5.0d0
      else if (name .eq. 'CH3OCH3') then
         a = 3.0d0
      else
         call bl_abort('setupbc: Unknown fuel type')
      end if

      Xt(oxidID) = 1.d0/(1.d0 + phi/a  + 0.79d0/0.21d0)
      Xt(fuelID) = phi * Xt(oxidID) / a
      Xt(bathID) = 1.d0 - Xt(fuelID) - Xt(oxidID)
      
      CALL CKXTY (Xt, Yt)
      
  end subroutine set_Y_from_Phi

end module user_defined_fcts_2d_module

