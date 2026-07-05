!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

!> \brief Check if two rectangular segments overlap

module overlap

   implicit none

   private
   public :: is_overlap

   interface is_overlap
      module procedure is_overlap_simple, is_overlap_per
   end interface

contains

!>
!! \brief Check if two blocks in a periodic domain are overlapping
!!
!! \details For shearing boxes divided at domain edges, an additional subroutine
!! (is_overlap_per_shear) may be needed and added to the is_overlap interface.
!<

   pure logical function is_overlap_per(this, other, periods) result(share)

      use constants,  only: xdim, ydim, zdim, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this    !< this box
      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: other   !< the other box
      integer(kind=8), dimension(xdim:zdim),        intent(in) :: periods !< where >0 then the direction is periodic with the given number of cells

      integer :: i, j, k
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: oth

      share = .false.
      do i = -1, 1
         if ((dom%has_dir(xdim) .or. periods(xdim)>0) .or. i==0) then
            oth(xdim, :) = other(xdim, :) + i*periods(xdim)
            do j = -1, 1
               if ((dom%has_dir(ydim) .or. periods(ydim)>0) .or. j==0) then
                  oth(ydim, :) = other(ydim, :) + j*periods(ydim)
                  do k = -1, 1
                     if ((dom%has_dir(zdim) .or. periods(zdim)>0) .or. k==0) then
                        oth(zdim, :) = other(zdim, :) + k*periods(zdim)
                        share = is_overlap_simple(this, oth) .or. share
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   end function is_overlap_per

!>
!! \brief Check if two blocks in a non-periodic domain are overlapping
!!
!! \details This routine does not handle periodic domain boundaries.
!! Use is_overlap_per to check overlap for boxes that cross periodic boundaries.
!<

   pure logical function is_overlap_simple(this, other) result(share)

      use constants,  only: xdim, zdim, LO, HI
      use domain,     only: dom

      implicit none

      integer(kind=8), dimension(:, :), intent(in) :: this  !< this box
      integer(kind=8), dimension(:, :), intent(in) :: other !< the other box

      integer :: d

      share = .true.
      do d = xdim, zdim
         if (dom%has_dir(d)) share = share .and. (other(d, LO) <= this(d, HI)) .and. (other(d, HI) >= this(d, LO))
      enddo

   end function is_overlap_simple

end module overlap
