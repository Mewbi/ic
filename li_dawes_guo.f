      subroutine inttocart(ic, ct)
! Transform internal coordinates to cartesian coordinates
!
!     Input:
!     ic, dim=6:
!     the internal coordinates, in the order:
!     R(HO), R(OH'), R(H'F), a(HOH'), a(OH'F), d(HOH'F)
!     Angles are in degree
!
!     Output:
!     ct, dim=3, 4:
!     The cartesian coordinates, in the order HHFO
!     The coordinates will be in the same unit as the input distances
!
!
! y
!  ^
!  |
!  |      H'----F
!  |     /
!       /
!  H---O   ---------> x 
! 
      implicit none
      real*8, intent(in) :: ic(6)
      real*8, intent(out) :: ct (3, 4)
      real*8, parameter :: torad=3.141592653589793238d0/180.0d0
      real*8 Fx, Fy
      real*8 cosHOH, sinHOH, cosOHF, sinOHF, cosDIH

      cosHOH = dcos(ic(4)*torad)
      sinHOH = dsin(ic(4)*torad)
      cosOHF = dcos(ic(5)*torad)
      sinOHF = dsin(ic(5)*torad)
      cosDIH = dcos(ic(6)*torad)
      
      ! H
      ct(1, 1) = 0.0d0
      ct(2, 1) = 0.0d0
      ct(3, 1) = 0.0d0
      ! O
      ct(1, 4) = ic(1)
      ct(2, 4) = 0.0d0
      ct(3, 4) = 0.0d0
      ! H'
      ct(1, 2) = ic(1) - ic(2)*cosHOH
      ct(2, 2) =         ic(2)*sinHOH
      ct(3, 2) = 0.0d0
      ! F
      Fx =  ic(3) * (cosHOH * cosOHF - sinHOH * cosDIH * sinOHF) 
      Fy = -ic(3) * (sinHOH * cosOHF + cosHOH * cosDIH * sinOHF)
      ct(1, 3) = Fx + ct(1, 2)
      ct(2, 3) = Fy + ct(2, 2)
      ct(3, 3) = dsqrt(ic(3)**2 - Fx**2 - Fy**2)
      end subroutine inttocart

      subroutine init()
      call pes_init
      end subroutine init
      
      subroutine pes(ic, V)
!     
!     Input:
!     ic, dim=6:
!     the internal coordinates, in the order:
!     R(HO), R(OH'), R(H'F), a(HOH'), a(OH'F), d(HOH'F)
!     Angles are in degree
!
!     Output:
!     The PES in kcal/mol
!
      use nnparam
      implicit none
      real*8, intent(out) :: V
      real*8, intent(in) :: ic(6)
      
      real*8, parameter :: tokcal=23.0605478306d0
      real*8 ct(3, 4)
      real*8 vpes, vpesa, vpesb, vpesc

      call inttocart(ic, ct)
      call fh2oNN(ct,vpes,vpesa,vpesb,vpesc)
      V = vpes * tokcal
      end subroutine pes
