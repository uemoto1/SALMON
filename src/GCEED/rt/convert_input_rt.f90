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
subroutine convert_input_rt(Ntime)
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root, comm_bcast
use inputoutput
use set_numcpu
implicit none
integer :: Ntime
real(8) :: dip_spacing

if(comm_is_root(nproc_id_global))then
   open(fh_namelist, file='.namelist.tmp', status='old')
end if

!===== namelist for group_fundamental =====
if(iwrite_projection==1.and.itwproj==-1)then
  write(*,*) "Please specify itwproj when iwrite_projection=1."
  stop
end if

!===== namelist for group_propagation =====
Ntime=nt
if(dt<=1.d-10)then
  write(*,*) "please set dt."
  stop
end if
if(Ntime==0)then
  write(*,*) "please set nt."
  stop
end if

!===== namelist for group_hartree =====
if(layout_multipole<=0.or.layout_multipole>=4)then
  stop "layout_multipole must be equal to 1 or 2 or 3."
else if(layout_multipole==3)then
  if(num_multipole_xyz(1)==0.and.num_multipole_xyz(2)==0.and.num_multipole_xyz(3)==0)then
    continue
  else if(num_multipole_xyz(1)<=0.or.num_multipole_xyz(2)<=0.or.num_multipole_xyz(3)<=0)then
    stop "num_multipole_xyz must be largar than 0 when they are not default values."
  end if
  if(num_multipole_xyz(1)==0.and.num_multipole_xyz(2)==0.and.num_multipole_xyz(3)==0)then
    dip_spacing = 8.d0/au_length_aa  ! approximate spacing of multipoles 
    num_multipole_xyz(:)=int((al(:)+dip_spacing)/dip_spacing-1.d-8)
  end if
end if

!===== namelist for group_others =====

if(comm_is_root(nproc_id_global))then
  if(iwdenoption/=0.and.iwdenoption/=1)then
    write(*,*)  'iwdenoption must be equal to 0 or 1.'
    stop
  end if
end if
if(iwdenoption==0)then
  iwdenstep=0
end if

if(comm_is_root(nproc_id_global))close(fh_namelist)

end subroutine convert_input_rt
