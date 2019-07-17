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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module band
    use math_constants, only : zi
    implicit none
  
  contains

    !> Calculate kgrid production
    subroutine calc_kgrid_prod(system, rgrid_lg, rgrid_mg, wf_info, wavefunction, nk1, nk2, nk3, ndk, prod_dk)
        use structures
        use salmon_parallel
        use salmon_communication
        implicit none
        type(s_system) , intent(in) :: system
        type(s_rgrid), intent(in) :: rgrid_lg, rgrid_mg
        type(s_wf_info), intent(in) :: wf_info
        type(s_wavefunction), intent(in) :: wavefunction
        integer, intent(in) :: nk1, nk2, nk3, ndk
        complex(8), intent(out) :: prod_dk(nk1, nk2, nk3, 1:3, 0:ndk, system%no, system%no)
        integer, parameter :: nrep = 2

        integer :: ik_tbl(nk1*nrep, nk2*nrep, nk3*nrep)
        integer :: im, ik1, ik2, ik3, jdk, ik, io
        complex(8), allocatable :: zwf(:, :, :, :, :)


        ! Create ik_tbl with periodic boundary condition:
        call create_ik_tbl()

        ! Retrieve entire wavefunction:
        call retrieve_entire_zwf()

        ! Calculate production <k,io|k+dk,jo> for all k,io,jo:
        do ik3 = 1, nk3
            do ik2 = 1, nk2
                do ik1 = 1, nk1
                    call calc_prod( &
                        & ik1, ik2, ik3, &
                        & ik1, ik2, ik3, &
                        & prod_dk(ik1, ik2, ik3, 1, 0, :, :))
                    prod_dk(ik1, ik2, ik3, 2, 0, :, :) = prod_dk(ik1, ik2, ik3, 1, 0, :, :)
                    prod_dk(ik1, ik2, ik3, 3, 0, :, :) = prod_dk(ik1, ik2, ik3, 1, 0, :, :)
                    do jdk = 1, ndk
                        call calc_prod( &
                            & ik1, ik2, ik3, &
                            & ik1 + jdk, ik2, ik3, &
                            & prod_dk(ik1, ik2, ik3, 1, jdk, :, :))
                        call calc_prod( &
                            & ik1, ik2, ik3, &
                            & ik1, ik2 + jdk, ik3, &
                            & prod_dk(ik1, ik2, ik3, 2, jdk, :, :))
                        call calc_prod( &
                            & ik1, ik2, ik3, &
                            & ik1, ik2, ik3 + jdk, &
                            & prod_dk(ik1, ik2, ik3, 3, jdk, :, :))
                    end do
                end do
            end do
        end do

        deallocate(zwf)

        return

    contains

    subroutine create_ik_tbl()
        implicit none
        integer :: ik_count
        integer :: ik1_o, ik2_o, ik3_o
        integer :: ik1_r, ik2_r, ik3_r
        integer :: ir1, ir2, ir3

        ik_count = 0

        do ik3_o = 1, nk3
            do ik2_o = 1, nk2
                do ik1_o = 1, nk1
                    ik_count = ik_count + 1
                    ! Replicate coordinates:
                    do ir3 = 1, nrep
                        ik3_r = ik3_o + ir3 * nk3
                        do ir2 = 1, nrep
                            ik2_r = ik2_o + ir2 * nk2
                            do ir1 = 1, nrep
                                ik1_r = ik1_o + ir1 * nk1
                                ik_tbl(ik1_r, ik2_r, ik3_r) = ik_count
                            end do
                        end do
                    end do
                end do
            end do
        end do

        return
    end subroutine create_ik_tbl


    subroutine retrieve_entire_zwf()
        use pack_unpack
        implicit none
        integer, parameter :: im = 1, ispin = 1
        complex(8), allocatable :: zwf_tmp(:, :, :, :, :)

        allocate(zwf_tmp( &
            & rgrid_lg%is(1):rgrid_lg%ie(1), &
            & rgrid_lg%is(2):rgrid_lg%ie(2), &
            & rgrid_lg%is(3):rgrid_lg%ie(3), &
            & system%no, system%nk))
        zwf_tmp = 0d0
        
        do ik = wf_info%ik_s, wf_info%ik_e
            do io = wf_info%io_s, wf_info%io_e
                call copy_data( &
                    wavefunction%zwf( &
                        & rgrid_mg%is(1):rgrid_mg%ie(1), &
                        & rgrid_mg%is(2):rgrid_mg%ie(2), &
                        & rgrid_mg%is(3):rgrid_mg%ie(3), &
                        & ispin, io, ik, im), &
                    zwf_tmp( &
                        & rgrid_mg%is(1):rgrid_mg%ie(1), &
                        & rgrid_mg%is(2):rgrid_mg%ie(2), &
                        & rgrid_mg%is(3):rgrid_mg%ie(3), &
                        & io, ik))
            end do
        end do
        
        allocate(zwf( &
            & rgrid_lg%is(1):rgrid_lg%ie(1), &
            & rgrid_lg%is(2):rgrid_lg%ie(2), &
            & rgrid_lg%is(3):rgrid_lg%ie(3), &
            & system%no, system%nk))

        call comm_summation( &
            & zwf_tmp, zwf, &
            & system%ngrid*system%no*system%nk, &
            & wf_info%icomm_rko)

        deallocate(zwf_tmp)

        return
    end subroutine retrieve_entire_zwf


    subroutine calc_prod(iik1, iik2, iik3, jjk1, jjk2, jjk3, prod)
        implicit none
        integer, intent(in) :: iik1, iik2, iik3
        integer, intent(in) :: jjk1, jjk2, jjk3
        complex(8), intent(out) :: prod(system%no, system%no)

        complex(8) zdotc ! From BLAS
        
        integer :: iik, jjk
        integer :: iio, jjo

        iik = ik_tbl(iik1, iik2, iik3)
        jjk = ik_tbl(jjk1, jjk2, jjk3)

        do iio = 1, system%no
            do jjo = 1, iio
                ! Compute inner product: <iik,iio|jjk,jjo>
                prod(iio, jjo) = system%Hvol * ZDOTC( &
                    & system%ngrid, &
                    & zwf(:, :, :, iio, iik), 1, &
                    & zwf(:, :, :, jjo, jjk), 1)
            end do
        end do

        do iio = 1, system%no
            do jjo = iio+1, system%no
                prod(iio, jjo) = conjg(prod(jjo, iio))
            end do
        end do

        return
    end subroutine calc_prod


    end subroutine calc_kgrid_prod
end module














  