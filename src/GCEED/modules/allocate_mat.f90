!
!  Copyright 2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
MODULE allocate_mat_sub

use inputoutput, only: iperiodic
use scf_data
implicit none

real(8), allocatable :: vecR(:,:,:,:)
real(8), allocatable :: vecR_tmp(:,:,:,:)

real(8), allocatable :: wk_s_h(:,:,:),wk2_s_h(:,:,:),lap_wk_s_h(:,:,:)

real(8), allocatable :: matbox_m(:,:,:),matbox_m2(:,:,:)
complex(8), allocatable :: cmatbox_m(:,:,:),cmatbox_m2(:,:,:)
real(8), allocatable :: matbox_l(:,:,:),matbox_l2(:,:,:)
complex(8), allocatable :: cmatbox_l(:,:,:),cmatbox_l2(:,:,:)

complex(8), allocatable :: zalpha2(:,:,:,:),zalpha3(:,:,:,:)

real(8),allocatable :: rgrad_wk(:,:,:,:,:,:)

complex(8),allocatable :: cgrad_wk(:,:,:,:,:,:)

real(8), allocatable :: rho_tmp(:,:,:)
real(8), allocatable :: rho_s_tmp(:,:,:,:)
real(8), allocatable :: vxc_tmp(:,:,:)
real(8), allocatable :: vxc_s_tmp(:,:,:,:)
real(8), allocatable :: eexc_tmp(:,:,:)
real(8), allocatable :: exc_dummy(:,:,:)
real(8), allocatable :: exc_dummy2(:,:,:,:)
real(8), allocatable :: exc_dummy3(:,:,:,:)

complex(8),allocatable :: rhoe_G(:)
complex(8),allocatable :: rhoe_G_tmp(:)

! FFTE routine
complex(8),allocatable :: A_FFTE(:,:,:), B_FFTE(:,:,:)
real(8),allocatable :: coef_poisson(:,:,:)
real(8),allocatable :: coef_poisson_FFTE(:,:,:)
real(8),allocatable :: A_FFTE_copy(:,:,:), A_FFTE_copy2(:,:,:)

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE allocate_mat
  implicit none

allocate (vecR(3,lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )

allocate (vecR_tmp(3,lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )

allocate (matbox_m(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )

allocate (matbox_m2(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )
allocate (cmatbox_m(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )
allocate (cmatbox_m2(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )

allocate (matbox_l(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )
allocate (matbox_l2(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )
allocate (cmatbox_l(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )
allocate (cmatbox_l2(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )

allocate (wk_s_h(ng_sta(1)-Ndh:ng_end(1)+Ndh,   &
             ng_sta(2)-Ndh:ng_end(2)+Ndh,   &
             ng_sta(3)-Ndh:ng_end(3)+Ndh))
allocate (wk2_s_h(ng_sta(1):ng_end(1),   &
              ng_sta(2):ng_end(2),   &
              ng_sta(3):ng_end(3)))
allocate (lap_wk_s_h(ng_sta(1):ng_end(1),   &
                 ng_sta(2):ng_end(2),   &
                 ng_sta(3):ng_end(3)))

if(iSCFRT==1.and.icalcforce==1)then
  select case(iperiodic)
  case(0)
    allocate(rgrad_wk(mg_sta(1):mg_end(1)+1,   &
                      mg_sta(2):mg_end(2),     &
                      mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end,3))
  case(3)
    allocate(cgrad_wk(mg_sta(1):mg_end(1)+1,   &
                      mg_sta(2):mg_end(2),     &
                      mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end,3))
  end select
else if(iSCFRT==2.and.icalcforce==1)then
  allocate(cgrad_wk(mg_sta(1):mg_end(1)+1,   &
                    mg_sta(2):mg_end(2),     &
                    mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end,3))
end if

allocate (rho_tmp(ng_num(1), ng_num(2), ng_num(3)))
allocate (rho_s_tmp(ng_num(1), ng_num(2), ng_num(3), 2))
allocate (vxc_tmp(ng_num(1), ng_num(2), ng_num(3)))
allocate (vxc_s_tmp(ng_num(1), ng_num(2), ng_num(3), 2))
allocate (eexc_tmp(ng_num(1), ng_num(2), ng_num(3)))
allocate (exc_dummy(ng_num(1), ng_num(2), ng_num(3)))
allocate (exc_dummy2(ng_num(1), ng_num(2), ng_num(3),2))
allocate (exc_dummy3(ng_num(1), ng_num(2), ng_num(3),3))

select case(iperiodic)
case(3)
  if(iSCFRT==2.and.iflag_hartree==4)then
    allocate(rhoe_G(lg_num(1)*lg_num(2)/NPUY*lg_num(3)/NPUZ))
  else
    allocate(rhoe_G(lg_num(1)*lg_num(2)*lg_num(3)))
    allocate(rhoe_G_tmp(lg_num(1)*lg_num(2)*lg_num(3)))
  end if

  if(iflag_hartree==4)then
    allocate(A_FFTE(lg_num(1),lg_num(2)/NPUY,lg_num(3)/NPUZ))
    allocate(B_FFTE(lg_num(1),lg_num(2)/NPUY,lg_num(3)/NPUZ))
    allocate(coef_poisson(lg_num(1),lg_num(2)/NPUY,lg_num(3)/NPUZ))
    allocate(A_FFTE_copy(lg_num(1),lg_num(2),lg_num(3)/NPUZ))
    allocate(A_FFTE_copy2(lg_num(1),lg_num(2),lg_num(3)/NPUZ))
  end if
end select

allocate(icoo1d(3,lg_num(1)*lg_num(2)*lg_num(3)))

END SUBROUTINE allocate_mat

!======================================================================

END MODULE allocate_mat_sub
