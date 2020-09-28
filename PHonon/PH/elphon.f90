!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE elphon()
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry, RY_TO_THZ, RY_TO_CMM1
  USE cell_base, ONLY : celldm, omega, ibrav, at, bg
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, amass
  USE gvecs, ONLY: doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : fft_interpolate
  USE noncollin_module, ONLY : nspin_mag, noncolin, m_loc
  USE lsda_mod, ONLY : nspin
  USE uspp,  ONLY: okvan
  USE paw_variables, ONLY : okpaw
  USE el_phon,  ONLY : done_elph
  USE dynmat, ONLY : dyn, w2
  USE modes,  ONLY : npert, nirr, u
  USE uspp_param, ONLY : nhm
  USE control_ph, ONLY : trans, xmldyn
  USE output,     ONLY : fildyn,fildvscf
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf, iuint3paw, lint3paw
  USE dfile_star,    ONLY : dvscf_star
  USE mp_bands,  ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE mp,        ONLY : mp_bcast
  USE io_global, ONLY: stdout
  USE lrus,   ONLY : int3, int3_nc, int3_paw
  USE qpoint, ONLY : xq
  USE dvscf_interpolate, ONLY : ldvscf_interpolate, dvscf_r2q
  USE ahc,    ONLY : elph_ahc
  !
  IMPLICIT NONE
  !
  INTEGER :: irr, imode0, ipert, is, npe
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  INTEGER :: i,j
  COMPLEX(DP), ALLOCATABLE :: dvscfin_all(:, :, :)
  !! dvscfin for all modes. Used when doing dvscf_r2q interpolation.
  COMPLEX(DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)
  COMPLEX(DP), allocatable :: phip (:, :, :, :)
  
  INTEGER :: ntyp_, nat_, ibrav_, nspin_mag_, mu, nu, na, nb, nta, ntb, nqs_
  REAL(DP) :: celldm_(6), w1
  CHARACTER(LEN=3) :: atm(ntyp)
   
  CALL start_clock ('elphon')

  if(dvscf_star%basis.eq.'cartesian') then
     write(stdout,*) 'Setting patterns to identity'
     u=(0.d0,0.d0)
     do irr=1,3*nat
        u(irr,irr)=(1.d0,0.d0)
     enddo
  endif
  !
  ! If ldvscf_interpolate, use Fourier interpolation instead of reading dVscf
  !
  IF (ldvscf_interpolate) THEN
    !
    WRITE (6, '(5x,a)') "Fourier interpolating dVscf"
    ALLOCATE(dvscfin_all(dfftp%nnr, nspin_mag, 3 * nat))
    CALL dvscf_r2q(xq, u, dvscfin_all)
    !
  ELSE
    WRITE (6, '(5x,a)') "Reading dVscf from file "//trim(fildvscf)
  ENDIF
  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     npe=npert(irr)
     ALLOCATE (dvscfin (dfftp%nnr, nspin_mag , npe) )
     IF (okvan) THEN
        ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, npe))
        IF (okpaw) ALLOCATE (int3_paw (nhm, nhm, nat, nspin_mag, npe))
        IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, npe))
     ENDIF
     DO ipert = 1, npe
        IF (ldvscf_interpolate) THEN
          dvscfin(:, :, ipert) = dvscfin_all(:, :, imode0 + ipert)
        ELSE
          CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
                             imode0 + ipert,  -1 )
        ENDIF
        IF (okpaw .AND. me_bgrp==0) &
             CALL davcio( int3_paw(:,:,:,:,ipert), lint3paw, &
                                          iuint3paw, imode0 + ipert, - 1 )
     END DO
     IF (okpaw) CALL mp_bcast(int3_paw, root_bgrp, intra_bgrp_comm)
     IF (doublegrid) THEN
        ALLOCATE (dvscfins (dffts%nnr, nspin_mag , npert(irr)) )
        DO is = 1, nspin_mag
           DO ipert = 1, npe
              CALL fft_interpolate (dfftp, dvscfin(:,is,ipert), dffts, dvscfins(:,is,ipert))
           ENDDO
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     CALL newdq (dvscfin, npert(irr))
     CALL elphel (irr, npert (irr), imode0, dvscfins)
     !
     imode0 = imode0 + npe
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
     IF (okvan) THEN
        DEALLOCATE (int3)
        IF (okpaw) DEALLOCATE (int3_paw)
        IF (noncolin) DEALLOCATE(int3_nc)
     ENDIF
  ENDDO
  !
  IF (ldvscf_interpolate) DEALLOCATE(dvscfin_all)
  !
  ! In AHC calculation, we do not need the dynamical matrix. So return here.
  IF (elph_ahc) THEN
     CALL stop_clock('elphon')
     RETURN
  ENDIF
  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
  !
  IF (.NOT.trans) THEN
     IF (.NOT. xmldyn) THEN
        WRITE (6, '(5x,a)') "Reading dynamics matrix from file "//trim(fildyn)
        CALL readmat (iudyn, ibrav, celldm, nat, ntyp, &
                      ityp, omega, amass, tau, xq, w2, dyn)
     ELSE
        allocate( phip(3,3,nat,nat) )
        CALL read_dyn_mat_param(fildyn, ntyp_, nat_)
        IF ( ntyp_ /= ntyp .OR. nat_ /= nat ) &
           CALL errore('elphon','uncorrect nat or ntyp',1)
          
        CALL read_dyn_mat_header(ntyp, nat, ibrav_, nspin_mag_, &
                 celldm_, at, bg, omega, atm, amass, tau, ityp, &
                 m_loc, nqs_)

        IF (ibrav_.NE.ibrav .OR. ABS ( celldm_ (1) - celldm (1) ) > 1.0d-5 &
             .OR. (nspin_mag_ /= nspin_mag ) ) CALL errore ('elphon', &
             'inconsistent data', 1)

        CALL read_dyn_mat(nat,1,xq,phip)
        !
        !  Diagonalize the dynamical matrix
        !

        
        DO i=1,3
           do na=1,nat
              nta = ityp (na)
              mu=3*(na-1)+i
              do j=1,3
                 do nb=1,nat
                   nu=3*(nb-1)+j
                   ntb = ityp (nb)
                   dyn (mu, nu) = phip (i, j, na, nb) / &
                     sqrt( amass(nta)*amass(ntb))/amu_ry
                 enddo
              enddo
           enddo
        enddo

        !
        CALL cdiagh (3 * nat, dyn, 3 * nat, w2, dyn)
        !
        ! divide by sqrt(mass) to get displacements
        !
        DO nu = 1, 3 * nat
           DO mu = 1, 3 * nat
              na = (mu - 1) / 3 + 1
              dyn (mu, nu) = dyn (mu, nu) / SQRT ( amu_ry * amass (ityp (na) ) )
           ENDDO
        ENDDO

        print*, "Entered here"
        print*, dyn(1,:)
        
        CALL read_dyn_mat_tail(nat)
  
        deallocate( phip )
     ENDIF
     !
     ! Write phonon frequency to stdout
     !
     WRITE( stdout, 8000) (xq (i), i = 1, 3)
     !
     DO nu = 1, 3 * nat
        w1 = SQRT( ABS( w2(nu) ) )
        if (w2(nu) < 0.d0) w1 = - w1
        WRITE( stdout, 8010) nu, w1 * RY_TO_THZ, w1 * RY_TO_CMM1
     ENDDO
     !
     WRITE( stdout, '(1x,74("*"))')
     !
  ENDIF ! .NOT. trans
  !
  CALL stop_clock ('elphon')
  !
8000 FORMAT(/,5x,'Diagonalizing the dynamical matrix', &
       &       //,5x,'q = ( ',3f14.9,' ) ',//,1x,74('*'))
8010 FORMAT   (5x,'freq (',i5,') =',f15.6,' [THz] =',f15.6,' [cm-1]')
  !
  RETURN
END SUBROUTINE elphon
!
!-----------------------------------------------------------------------
SUBROUTINE readmat (iudyn, ibrav, celldm, nat, ntyp, ityp, omega, &
     amass, tau, q, w2, dyn)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : amu_ry
  IMPLICIT NONE
  ! Input
  INTEGER :: iudyn, ibrav, nat, ntyp, ityp (nat)
  REAL(DP) :: celldm (6), amass (ntyp), tau (3, nat), q (3), &
       omega
  ! output
  REAL(DP) :: w2 (3 * nat)
  COMPLEX(DP) :: dyn (3 * nat, 3 * nat)
  ! local (control variables)
  INTEGER :: ntyp_, nat_, ibrav_, ityp_
  REAL(DP) :: celldm_ (6), amass_, tau_ (3), q_ (3)
  ! local
  REAL(DP) :: dynr (2, 3, nat, 3, nat)
  CHARACTER(len=80) :: line
  CHARACTER(len=3)  :: atm
  INTEGER :: nt, na, nb, naa, nbb, nu, mu, i, j
  !
  !
  REWIND (iudyn)
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
  IF ( ntyp.NE.ntyp_ .OR. nat.NE.nat_ .OR.ibrav_.NE.ibrav .OR. &
       ABS ( celldm_ (1) - celldm (1) ) > 1.0d-5) &
          CALL errore ('readmat', 'inconsistent data', 1)
  IF ( ibrav_ == 0 ) THEN
     READ (iudyn, '(a)') line
     READ (iudyn, '(a)') line
     READ (iudyn, '(a)') line
     READ (iudyn, '(a)') line
  END IF
  DO nt = 1, ntyp
     READ (iudyn, * ) i, atm, amass_
     IF ( nt.NE.i .OR. ABS (amass_ - amu_ry*amass (nt) ) > 1.0d-5) &
        CALL errore ( 'readmat', 'inconsistent data', 1 + nt)
  ENDDO
  DO na = 1, nat
     READ (iudyn, * ) i, ityp_, tau_
     IF (na.NE.i.OR.ityp_.NE.ityp (na) ) CALL errore ('readmat', &
          'inconsistent data', 10 + na)
  ENDDO
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (iudyn, '(a)') line
  READ (line (11:80), * ) (q_ (i), i = 1, 3)
  READ (iudyn, '(a)') line
  DO na = 1, nat
     DO nb = 1, nat
        READ (iudyn, * ) naa, nbb
        IF (na.NE.naa.OR.nb.NE.nbb) CALL errore ('readmat', 'error reading &
             &file', nb)
        READ (iudyn, * ) ( (dynr (1, i, na, j, nb), dynr (2, i, na, j, nb) &
             , j = 1, 3), i = 1, 3)
     ENDDO
  ENDDO
  !
  ! divide the dynamical matrix by the (input) masses (in amu)
  !
  DO nb = 1, nat
     DO j = 1, 3
        DO na = 1, nat
           DO i = 1, 3
              dynr (1, i, na, j, nb) = dynr (1, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) ) / amu_ry
              dynr (2, i, na, j, nb) = dynr (2, i, na, j, nb) / SQRT (amass ( &
                   ityp (na) ) * amass (ityp (nb) ) ) / amu_ry
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  ! solve the eigenvalue problem.
  ! NOTA BENE: eigenvectors are overwritten on dyn
  !
  CALL cdiagh (3 * nat, dynr, 3 * nat, w2, dyn)
  !
  ! divide by sqrt(mass) to get displacements
  !
  DO nu = 1, 3 * nat
     DO mu = 1, 3 * nat
        na = (mu - 1) / 3 + 1
        dyn (mu, nu) = dyn (mu, nu) / SQRT ( amu_ry * amass (ityp (na) ) )
     ENDDO
  ENDDO
  !
  !
  RETURN
END SUBROUTINE readmat
!
!-----------------------------------------------------------------------
SUBROUTINE elphel (irr, npe, imode0, dvscfins)
  !-----------------------------------------------------------------------
  !
  !      Calculation of the electron-phonon matrix elements el_ph_mat
  !         <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !      Original routine written by Francesco Mauri
  !      Modified by A. Floris and I. Timrov to include Hubbard U (01.10.2018)
  !
  USE kinds,      ONLY : DP
  USE fft_base,   ONLY : dffts
  USE ions_base,  ONLY : nat, ityp
  USE control_flags,  ONLY : iverbosity
  USE wavefunctions,  ONLY : evc
  USE buffers,    ONLY : get_buffer, save_buffer
  USE klist,      ONLY : xk, ngk, igk_k
  USE lsda_mod,   ONLY : lsda, current_spin, isk, nspin
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE wvfct,      ONLY : nbnd, npwx
  USE uspp,       ONLY : vkb
  USE el_phon,    ONLY : el_ph_mat, el_ph_mat_rec, el_ph_mat_rec_col, &
                         comp_elph, done_elph, elph_nbnd_min, elph_nbnd_max
  USE modes,      ONLY : u, nmodes
  USE units_ph,   ONLY : iubar, lrbar, iundnsscf, iudvpsi, lrdvpsi
  USE units_lr,   ONLY : iuwfc, lrwfc
  USE control_ph, ONLY : trans, current_iq
  USE ph_restart, ONLY : ph_writefile
  USE spin_orb,   ONLY : domag
  USE mp_bands,   ONLY : intra_bgrp_comm, ntask_groups
  USE mp_pools,   ONLY : npool
  USE mp,         ONLY : mp_sum, mp_bcast
  USE mp_world,   ONLY : world_comm
  USE elph_tetra_mod, ONLY : elph_tetra
  USE eqv,        ONLY : dvpsi, evq
  USE qpoint,     ONLY : nksq, ikks, ikqs, nksqtot
  USE control_lr, ONLY : lgamma
  USE fft_helper_subroutines
  USE ldaU,       ONLY : lda_plus_u, Hubbard_lmax
  USE ldaU_ph,    ONLY : dnsscf_all_modes, dnsscf
  USE io_global,  ONLY : ionode, ionode_id
  USE io_files,   ONLY : seqopn
  USE lrus,       ONLY : becp1
  USE phus,       ONLY : alphap
  USE ahc,        ONLY : elph_ahc, ib_ahc_gauge_min, ib_ahc_gauge_max

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: irr, npe, imode0
  COMPLEX(DP), INTENT(IN) :: dvscfins (dffts%nnr, nspin_mag, npe)
  ! LOCAL variables
  INTEGER :: npw, npwq, nrec, ik, ikk, ikq, ipert, mode, ibnd, jbnd, ir, ig, &
             ipol, ios, ierr, nrec_ahc
  COMPLEX(DP) , ALLOCATABLE :: aux1 (:,:), elphmat (:,:,:), tg_dv(:,:), &
                               tg_psic(:,:), aux2(:,:)
  INTEGER :: v_siz, incr
  LOGICAL :: exst
  COMPLEX(DP), EXTERNAL :: zdotc
  integer :: ibnd_fst, ibnd_lst
  !
  CALL start_clock('elphel')
  !
  IF (elph_ahc) THEN
     ibnd_fst = ib_ahc_gauge_min
     ibnd_lst = ib_ahc_gauge_max
  elseif(elph_tetra == 0) then
     ibnd_fst = 1
     ibnd_lst = nbnd
  else
     ibnd_fst = elph_nbnd_min
     ibnd_lst = elph_nbnd_max
  end if
  !
  IF (.NOT. comp_elph(irr) .OR. done_elph(irr)) RETURN

  ALLOCATE (aux1    (dffts%nnr, npol))
  ALLOCATE (elphmat ( nbnd , nbnd , npe))
  ALLOCATE( el_ph_mat_rec (nbnd,nbnd,nksq,npe) )
  el_ph_mat_rec=(0.0_DP,0.0_DP)
  ALLOCATE (aux2(npwx*npol, nbnd))
  incr=1
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz =  dffts%nnr_tg
     ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = fftx_ntgrp(dffts)
     !
  ENDIF
  !
  ! DFPT+U case
  !
  IF (lda_plus_u .AND. .NOT.trans) THEN
     !
     ! Allocate and read dnsscf_all_modes from file 
     !
     ALLOCATE (dnsscf_all_modes(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nmodes))
     dnsscf_all_modes = (0.d0, 0.d0)
     !
     IF (ionode) READ(iundnsscf,*) dnsscf_all_modes
     CALL mp_bcast(dnsscf_all_modes, ionode_id, world_comm)
     REWIND(iundnsscf)
     !  
     ! Check whether the re-read is correct
     !
     IF (iverbosity==1) CALL elphel_read_dnsscf_check() 
     !
     ! Allocate dnsscf
     !
     ALLOCATE (dnsscf(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, npe))
     dnsscf = (0.d0, 0.d0)
     !
  ENDIF
  !
  !  Start the loops over the k-points
  !
  DO ik = 1, nksq
     !
     !  ik = counter of k-points with vector k
     !  ikk= index of k-point with vector k
     !  ikq= index of k-point with vector k+q
     !       k and k+q are alternated if q!=0, are the same if q=0
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     IF (lsda) current_spin = isk (ikk)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     !
     CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
     !
     ! read unperturbed wavefuctions psi(k) and psi(k+q)
     !
     IF (nksq.GT.1) THEN
        IF (lgamma) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikk)
        ELSE
           CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           CALL get_buffer (evq, lrwfc, iuwfc, ikq)
        ENDIF
     ENDIF
     !
     DO ipert = 1, npe
        nrec = (ipert - 1) * nksq + ik
        !
        !  dvbare_q*psi_kpoint is read from file (if available) or recalculated
        !
        IF (trans) THEN
           CALL get_buffer (dvpsi, lrbar, iubar, nrec)
        ELSE
           mode = imode0 + ipert
           ! FIXME: .false. or .true. ???
           CALL dvqpsi_us (ik, u (1, mode), .FALSE., becp1, alphap)
           !
           ! DFPT+U: calculate the bare derivative of the Hubbard potential in el-ph
           !
           IF (lda_plus_u) CALL dvqhub_barepsi_us (ik, u(1,mode)) 
           !
        ENDIF
        !
        ! calculate dvscf_q*psi_k
        !
        IF ( dffts%has_task_groups ) THEN
           IF (noncolin) THEN
              CALL tg_cgather( dffts, dvscfins(:,1,ipert), tg_dv(:,1))
              IF (domag) THEN
                 DO ipol=2,4
                    CALL tg_cgather( dffts,  dvscfins(:,ipol,ipert), tg_dv(:,ipol))
                 ENDDO
              ENDIF
           ELSE
              CALL tg_cgather( dffts, dvscfins(:,current_spin,ipert), tg_dv(:,1))
           ENDIF
        ENDIF
        aux2=(0.0_DP,0.0_DP)
        DO ibnd = ibnd_fst, ibnd_lst, incr
           IF ( dffts%has_task_groups ) THEN
              CALL cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, nbnd )
              CALL apply_dpot(v_siz, tg_psic, tg_dv, 1)
              CALL cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, nbnd)
           ELSE
              CALL cft_wave (ik, evc(1, ibnd), aux1, +1)
              CALL apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipert), current_spin)
              CALL cft_wave (ik, aux2(1, ibnd), aux1, -1)
           ENDIF
        ENDDO
        dvpsi=dvpsi+aux2
        !
        CALL adddvscf (ipert, ik)
        !
        ! DFPT+U: add to dvpsi the scf part of the perturbed Hubbard potential 
        !
        IF (lda_plus_u) THEN
           IF (.NOT.trans) dnsscf(:,:,:,:,ipert) = dnsscf_all_modes(:,:,:,:,mode)
           CALL adddvhubscf (ipert, ik)
        ENDIF
        !
        ! If doing Allen-Heine-Cardona (AHC) calculation, we need dvpsi
        ! later. So, write to buffer.
        !
        IF (elph_ahc) THEN
           nrec_ahc = (ik - 1) * nmodes + ipert + imode0
           CALL save_buffer(dvpsi(1, ibnd_fst), lrdvpsi, iudvpsi, nrec_ahc)
           !
           ! If elph_ahc, the matrix elements are computed in ahc.f90
           CYCLE
           !
        ENDIF
        !
        ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
        !
        DO ibnd = ibnd_fst, ibnd_lst
           DO jbnd = ibnd_fst, ibnd_lst
              elphmat (jbnd, ibnd, ipert) = zdotc (npwq, evq (1, jbnd), 1, &
                   dvpsi (1, ibnd), 1)
              IF (noncolin) &
                 elphmat (jbnd, ibnd, ipert) = elphmat (jbnd, ibnd, ipert)+ &
                   zdotc (npwq, evq(npwx+1,jbnd),1,dvpsi(npwx+1,ibnd), 1)
           ENDDO
        ENDDO
     ENDDO
     !
     ! If elph_ahc, the matrix elements are computed in ahc.f90
     IF (elph_ahc) CYCLE
     !
     CALL mp_sum (elphmat, intra_bgrp_comm)
     !
     !  save all e-ph matrix elements into el_ph_mat
     !
     DO ipert = 1, npe
        DO jbnd = ibnd_fst, ibnd_lst
           DO ibnd = ibnd_fst, ibnd_lst
              el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = elphmat (ibnd, jbnd, ipert)
              el_ph_mat_rec (ibnd, jbnd, ik, ipert ) = elphmat (ibnd, jbnd, ipert)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  done_elph(irr)=.TRUE.
  if(elph_tetra == 0 .AND. .NOT. elph_ahc) then
     IF (npool>1) THEN
        ALLOCATE(el_ph_mat_rec_col(nbnd,nbnd,nksqtot,npe))
        CALL el_ph_collect(npe,el_ph_mat_rec,el_ph_mat_rec_col,nksqtot,nksq)
     ELSE
        el_ph_mat_rec_col => el_ph_mat_rec
     ENDIF
     CALL ph_writefile('el_phon',current_iq,irr,ierr)
     IF (npool > 1) DEALLOCATE(el_ph_mat_rec_col)
  end if
  DEALLOCATE(el_ph_mat_rec)
  !
  DEALLOCATE (elphmat)
  DEALLOCATE (aux1)
  DEALLOCATE (aux2)
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE( tg_dv )
     DEALLOCATE( tg_psic )
  ENDIF
  !
  IF (lda_plus_u .AND. .NOT.trans) THEN
     DEALLOCATE (dnsscf_all_modes)
     DEALLOCATE (dnsscf)
  ENDIF
  !
  CALL stop_clock('elphel')
  !
  RETURN
  !
END SUBROUTINE elphel
!
!------------------------------------------------------------------------
SUBROUTINE elphel_read_dnsscf_check()
  !
  ! DFPT+U: This subroutine checks whether dnsscf_all_modes was 
  !         read correctly from file.
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ityp
  USE modes,      ONLY : u, nmodes
  USE lsda_mod,   ONLY : nspin
  USE ldaU,       ONLY : Hubbard_l, is_hubbard, Hubbard_lmax
  USE ldaU_ph,    ONLY : dnsscf_all_modes
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: dnsscf_all_modes_cart(:,:,:,:,:)
  INTEGER :: na_icart, nah, is, m1, m2, na, icart, nt, na_icar, imode
  !
  ALLOCATE(dnsscf_all_modes_cart (2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nmodes))
  dnsscf_all_modes_cart = (0.d0, 0.d0)
  !
  ! Transform dnsscf_all_modes from pattern to cartesian coordinates
  !
  DO na_icart = 1, 3*nat
     DO imode = 1, nmodes
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 DO m1 = 1, 2*Hubbard_l(nt) + 1
                    DO m2 = 1, 2*Hubbard_l(nt) + 1
                       !
                       dnsscf_all_modes_cart (m1, m2, is, nah, na_icart) = &
                              dnsscf_all_modes_cart (m1, m2, is, nah, na_icart) + &
                              dnsscf_all_modes (m1, m2, is, nah, imode) * &
                              CONJG(u(na_icart,imode))
                       !
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  ! Write dnsscf in cartesian coordinates
  !
  WRITE(stdout,*)
  WRITE(stdout,*) 'DNS_SCF SYMMETRIZED IN CARTESIAN COORDINATES'
  !
  DO na = 1, nat
     DO icart = 1, 3
        WRITE( stdout,'(a,1x,i2,2x,a,1x,i2)') 'displaced atom L =', na, 'ipol=', icart
        na_icart = 3*(na-1) + icart
        DO nah = 1, nat
           nt = ityp(nah)
           IF (is_hubbard(nt)) THEN
              DO is = 1, nspin
                 WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                 DO m1 = 1, 2*Hubbard_l(nt) + 1
                    WRITE(stdout,'(14(f15.10,1x))') dnsscf_all_modes_cart (m1,:,is,nah,na_icart)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  WRITE(stdout,*)
  ! 
  DEALLOCATE(dnsscf_all_modes_cart) 
  !
  RETURN
  !
END SUBROUTINE elphel_read_dnsscf_check
!------------------------------------------------------------------------

!------------------------------------------------------------------------
SUBROUTINE elphsum ( )
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri, modified by PG
  !      New version by  Malgorzata Wierzbowska
  !
  USE kinds,       ONLY : DP
  USE constants,   ONLY : pi, rytoev, ry_to_cmm1, ry_to_ghz, degspin
  USE ions_base,   ONLY : nat, ityp, tau
  USE cell_base,   ONLY : at, bg
  USE lsda_mod,    ONLY: isk, nspin
  USE klist,       ONLY: nks, nkstot, xk, wk, nelec
  USE start_k,     ONLY: nk1, nk2, nk3
  USE symm_base,   ONLY: s, irt, nsym, invs
  USE noncollin_module, ONLY: nspin_lsda, nspin_mag
  USE wvfct,       ONLY: nbnd, et
  USE parameters,  ONLY : npk
  USE el_phon,     ONLY : el_ph_mat, done_elph, el_ph_nsigma, el_ph_ngauss, &
                          el_ph_sigma
  USE modes,       ONLY : u, nirr
  USE dynmat,      ONLY : dyn, w2
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE parallel_include
  USE mp_pools,    ONLY : my_pool_id, npool, kunit
  USE mp_images,   ONLY : intra_image_comm, me_image, nproc_image
  USE mp,          ONLY : mp_bcast
  USE control_ph,  ONLY : tmp_dir_phq, xmldyn, current_iq
  USE save_ph,     ONLY : tmp_dir_save
  USE io_files,    ONLY : prefix, tmp_dir, seqopn, create_directory
  !
  USE lr_symm_base, ONLY : minus_q, nsymq, rtau
  USE qpoint,       ONLY : xq, nksq
  USE control_lr,   ONLY : lgamma
  !
  IMPLICIT NONE
  ! epsw = 20 cm^-1, in Ry
  REAL(DP), PARAMETER :: epsw = 20.d0 / ry_to_cmm1
  REAL(DP), PARAMETER :: eps = 1.0d-6
  !
  INTEGER :: iuna2Fsave  = 40
  !
  REAL(DP), allocatable :: gam(:,:), lamb(:,:)
  !
  ! Quantities ending with "fit" are relative to the "dense" grid
  !
  REAL(DP), allocatable :: xkfit(:,:)
  REAL(DP), allocatable, target :: etfit(:,:), wkfit(:)
  INTEGER :: nksfit, nk1fit, nk2fit, nk3fit, nkfit, nksfit_real
  INTEGER, allocatable :: eqkfit(:), eqqfit(:), sfit(:)
  !
  integer :: nq, isq (48), imq
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)
  real(DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  ! workspace used for symmetrisation
  !
  COMPLEX(DP), allocatable :: g1(:,:,:), g2(:,:,:), g0(:,:), gf(:,:,:)
  COMPLEX(DP), allocatable :: point(:), noint(:)
  COMPLEX(DP) :: dyn22(3*nat,3*nat), ctemp
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, iuelph, ios, i,k,j, ii, jj
  INTEGER :: nkBZ, nti, ntj, ntk, nkr, itemp1, itemp2, nn, &
       qx,qy,qz,iq,jq,kq
  INTEGER, ALLOCATABLE :: eqBZ(:), sBZ(:)
  REAL(DP) :: weight, wqa, w0g1, w0g2, degauss1, effit1, dosef, &
       ef1, lambda, gamma
  REAL(DP), ALLOCATABLE :: deg(:), effit(:), dosfit(:)
  REAL(DP) :: etk, etq
  REAL(DP), EXTERNAL :: dos_ef, efermig, w0gauss
  character(len=80) :: name
  LOGICAL  :: exst, xmldyn_save
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat)

  COMPLEX(DP), POINTER :: el_ph_mat_collect(:,:,:,:)
  REAL(DP), ALLOCATABLE :: xk_collect(:,:)
  REAL(DP), POINTER :: wkfit_dist(:), etfit_dist(:,:)
  INTEGER :: nksfit_dist, rest, kunit_save
  INTEGER :: nks_real, ispin, nksqtot, irr, ierr
  CHARACTER(LEN=256) :: elph_dir
  CHARACTER(LEN=6) :: int_to_char
  !
  !
  !
  !  If the electron phonon matrix elements have not been calculated for
  !  all representations this routine exit
  !
  DO irr=1,nirr
     IF (.NOT.done_elph(irr)) RETURN
  ENDDO

  CALL start_clock('elphsum')

  elph_dir='elph_dir/'
  IF (ionode) INQUIRE(file=TRIM(elph_dir), EXIST=exst)
  CALL mp_bcast(exst, ionode_id, intra_image_comm) 
  IF (.NOT.exst) CALL create_directory( elph_dir )
  WRITE (6, '(5x,"electron-phonon interaction  ..."/)')
  ngauss1 = 0

  ALLOCATE(xk_collect(3,nkstot))

  ALLOCATE(deg(el_ph_nsigma))
  ALLOCATE(effit(el_ph_nsigma))
  ALLOCATE(dosfit(el_ph_nsigma))

  IF (npool==1) THEN
!
!  no pool, just copy old variables on the new ones
!
     nksqtot=nksq
     xk_collect(:,1:nks) = xk(:,1:nks)
     el_ph_mat_collect => el_ph_mat
  ELSE  
!
!  pools, allocate new variables and collect the results. All the rest
!  remain unchanged.
!
     IF (lgamma) THEN
        nksqtot=nkstot
     ELSE
        nksqtot=nkstot/2
     ENDIF
     CALL poolcollect(3, nks, xk, nkstot, xk_collect)
     ALLOCATE(el_ph_mat_collect(nbnd,nbnd,nksqtot,3*nat))
     ! FIXME: this routine should be replaced by a generic routine
     CALL el_ph_collect(3*nat,el_ph_mat,el_ph_mat_collect,nksqtot,nksq)
  ENDIF
  !
  ! read eigenvalues for the dense grid
  ! FIXME: this should be done from the xml file, not from a specialized file
  ! parallel case: only first node reads
  !
  IF ( ionode ) THEN
     tmp_dir=tmp_dir_save
     CALL seqopn( iuna2Fsave, 'a2Fsave', 'FORMATTED', exst )
     tmp_dir=tmp_dir_phq
     READ(iuna2Fsave,*) ibnd, nksfit
  END IF
  !
  CALL mp_bcast (ibnd, ionode_id, intra_image_comm)
  CALL mp_bcast (nksfit, ionode_id, intra_image_comm)
  if ( ibnd /= nbnd ) call errore('elphsum','wrong file read',iuna2Fsave)
  allocate (etfit(nbnd,nksfit), xkfit(3,nksfit), wkfit(nksfit))
  !
  IF ( ionode ) THEN
     READ(iuna2Fsave,*) etfit
     READ(iuna2Fsave,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
     READ(iuna2Fsave,*) wkfit
     READ(iuna2Fsave,*) nk1fit, nk2fit, nk3fit
     CLOSE( UNIT = iuna2Fsave, STATUS = 'KEEP' )
  END IF
  !
  ! broadcast all variables read
  !
  CALL mp_bcast (etfit, ionode_id, intra_image_comm)
  CALL mp_bcast (xkfit, ionode_id, intra_image_comm)
  CALL mp_bcast (wkfit, ionode_id, intra_image_comm)
  CALL mp_bcast (nk1fit, ionode_id, intra_image_comm)
  CALL mp_bcast (nk2fit, ionode_id, intra_image_comm)
  CALL mp_bcast (nk3fit, ionode_id, intra_image_comm)
  !
  nkfit=nk1fit*nk2fit*nk3fit
  !
  ! efermig and dos_ef require scattered points and eigenvalues
  ! isk is neither read nor used. phonon with two Fermi energies is
  ! not yet implemented.
  !
  nksfit_dist  = ( nksfit / npool )
  rest = ( nksfit - nksfit_dist * npool ) 
  IF ( ( my_pool_id + 1 ) <= rest ) nksfit_dist = nksfit_dist + 1
  kunit_save=kunit
  kunit=1

#if defined(__MPI)
  ALLOCATE(etfit_dist(nbnd,nksfit_dist))
  ALLOCATE(wkfit_dist(nksfit_dist))
  CALL poolscatter( 1, nksfit, wkfit, nksfit_dist, wkfit_dist )
  CALL poolscatter( nbnd, nksfit, etfit, nksfit_dist, etfit_dist )
#else
   wkfit_dist => wkfit
   etfit_dist => etfit
#endif
  !
  do isig=1,el_ph_nsigma
     !
     ! recalculate Ef = effit and DOS at Ef N(Ef) = dosfit using dense grid
     ! for value "deg" of gaussian broadening
     !
     deg(isig) = isig * el_ph_sigma
     !
     effit(isig) = efermig &
          ( etfit_dist, nbnd, nksfit_dist, nelec, wkfit_dist, &
              deg(isig), ngauss1, 0, isk)
     dosfit(isig) = dos_ef ( ngauss1, deg(isig), effit(isig), etfit_dist, &
          wkfit_dist, nksfit_dist, nbnd) / 2.0d0
  enddo
#if defined(__MPI)
  DEALLOCATE(etfit_dist)
  DEALLOCATE(wkfit_dist)
#endif
  kunit=kunit_save
  allocate (eqkfit(nkfit), eqqfit(nkfit), sfit(nkfit))
  !
  ! map k-points in the IBZ to k-points in the complete uniform grid
  !
  nksfit_real=nksfit/nspin_lsda
  call lint ( nsym, s, .true., at, bg, npk, 0,0,0, &
       nk1fit,nk2fit,nk3fit, nksfit_real, xkfit, 1, nkfit, eqkfit, sfit)
  deallocate (sfit, xkfit, wkfit)
  !
  ! find epsilon(k+q) in the dense grid
  !
  call cryst_to_cart (1, xq, at, -1)
  qx = nint(nk1fit*xq(1))
  if (abs(qx-nk1fit*xq(1)) > eps) &
       call errore('elphsum','q is not a vector in the dense grid',1)
  if (qx < 0) qx = qx + nk1fit
  if (qx > nk1fit) qx = qx - nk1fit
  qy = nint(nk2fit*xq(2))
  if (abs(qy-nk2fit*xq(2)) > eps) &
       call errore('elphsum','q is not a vector in the dense grid',2)
  if (qy < 0) qy = qy + nk2fit
  if (qy > nk2fit) qy = qy - nk2fit
  qz = nint(nk3fit*xq(3))
  if (abs(qz-nk3fit*xq(3)) > eps) &
       call errore('elphsum','q is not a vector in the dense grid',3)
  if (qz < 0) qz = qz + nk3fit
  if (qz > nk3fit) qz = qz - nk3fit
  call cryst_to_cart (1, xq, bg, 1)
  !
  eqqfit(:) = 0
  do i=1,nk1fit
     do j=1,nk2fit
        do k=1,nk3fit
           ik = k-1 + (j-1)*nk3fit + (i-1)*nk2fit*nk3fit + 1
           iq = i+qx
           if (iq > nk1fit) iq = iq - nk1fit
           jq = j+qy
           if (jq > nk2fit) jq = jq - nk2fit
           kq = k+qz
           if (kq > nk3fit) kq = kq - nk3fit
           nn = (kq-1)+(jq-1)*nk3fit+(iq-1)*nk2fit*nk3fit + 1
           eqqfit(ik) = eqkfit(nn)
        enddo
     enddo
  enddo
  !
  ! calculate the electron-phonon coefficient using the dense grid
  !
  nti  = nk1fit/nk1
  ntj  = nk2fit/nk2
  ntk  = nk3fit/nk3
  nkBZ  = nk1*nk2*nk3
  allocate (eqBZ(nkBZ), sBZ(nkBZ))
  !
  nks_real=nkstot/nspin_lsda
  IF ( lgamma ) THEN
     call lint ( nsymq, s, minus_q, at, bg, npk, 0,0,0, &
          nk1,nk2,nk3, nks_real, xk_collect, 1, nkBZ, eqBZ, sBZ)
  ELSE
     call lint ( nsymq, s, minus_q, at, bg, npk, 0,0,0, &
          nk1,nk2,nk3, nks_real, xk_collect, 2, nkBZ, eqBZ, sBZ)
  END IF
  !
  allocate (gf(3*nat,3*nat,el_ph_nsigma))
  gf = (0.0d0,0.0d0)
  !
  wqa  = 1.0d0/nkfit
  IF (nspin==1) wqa=degspin*wqa
  !
#if defined(__MPI)
  do ibnd = me_image+1, nbnd, nproc_image
#else
  do ibnd = 1, nbnd
#endif
     do jbnd = 1, nbnd
        allocate (g2(nkBZ*nspin_lsda,3*nat,3*nat))
        allocate (g1(nksqtot,3*nat,3*nat))
        do ik = 1, nksqtot
           do ii = 1, 3*nat
              do jj = 1, 3*nat
                 g1(ik,ii,jj)=CONJG(el_ph_mat_collect(jbnd,ibnd,ik,ii))* &
                      el_ph_mat_collect(jbnd,ibnd,ik,jj)
              enddo    ! ipert
           enddo    !jpert
        enddo   ! ik
        !
        allocate (g0(3*nat,3*nat))
        do i=1,nk1
           do j=1,nk2
              do k=1,nk3
                 do ispin=1,nspin_lsda
                    nn = k-1 + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                    itemp1 = eqBZ(nn)
                    if (ispin==2) itemp1=itemp1+nksqtot/2
                    g0(:,:) = g1(itemp1,:,:)
                    itemp2 = sBZ(nn)
                    call symm ( g0, u, xq, s, itemp2, rtau, irt, &
                         at, bg, nat)
                    if (ispin==2) nn=nn+nkBZ 
                    g2(nn,:,:) = g0(:,:)
                 enddo
              enddo ! k
           enddo !j
        enddo !i
        deallocate (g0)
        deallocate (g1)
        !
        allocate ( point(nkBZ), noint(nkfit) )
        do jpert = 1, 3 * nat
           do ipert = 1, 3 * nat
              do ispin=1,nspin_lsda
              !
              point(1:nkBZ) = &
                  g2(1+nkBZ*(ispin-1):nkBZ+nkBZ*(ispin-1),ipert,jpert)
              !
              CALL clinear(nk1,nk2,nk3,nti,ntj,ntk,point,noint)
              !
              do isig = 1,el_ph_nsigma
                 degauss1 = deg(isig)
                 effit1   = effit(isig)
                 ctemp    = 0
                 do ik=1,nkfit
                    etk   = etfit(ibnd,eqkfit(ik)+nksfit*(ispin-1)/2)
                    etq   = etfit(jbnd,eqqfit(ik)+nksfit*(ispin-1)/2)
                    ctemp = ctemp &
                          + exp(-((effit1-etk)**2 + (effit1-etq)**2)/degauss1**2)*noint(ik)
                 enddo
                 gf(ipert,jpert,isig) = gf(ipert,jpert,isig) + &
                      ctemp * wqa / (degauss1**2) / pi 
              enddo ! isig
              enddo ! ispin
           enddo    ! ipert
        enddo    !jpert
        deallocate (point, noint)
        deallocate (g2)
        !
     enddo    ! ibnd
  enddo    ! jbnd

#if defined(__MPI)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,gf,3*nat*3*nat*el_ph_nsigma, &
                     MPI_DOUBLE_COMPLEX,MPI_SUM,intra_image_comm,ierr)
#endif

  deallocate (eqqfit, eqkfit)
  deallocate (etfit)
  deallocate (eqBZ, sBZ)
!
  allocate (gam(3*nat,el_ph_nsigma), lamb(3*nat,el_ph_nsigma))
  lamb(:,:) = 0.0d0
  gam (:,:) = 0.0d0
  do isig= 1, el_ph_nsigma
     do nu = 1,3*nat
        gam(nu,isig) = 0.0d0
        do mu = 1, 3 * nat
           do vu = 1, 3 * nat
              gam(nu,isig) = gam(nu,isig) + DBLE(conjg(dyn(mu,nu)) * &
                   gf(mu,vu,isig) * dyn(vu,nu))
           enddo
        enddo
        gam(nu,isig) = gam(nu,isig) *  pi/2.0d0
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        if (sqrt(abs(w2(nu))) > epsw) then
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lamb(nu,isig) = gam(nu,isig)/pi/w2(nu)/dosfit(isig)
        else
           lamb(nu,isig) = 0.0d0
        endif
        gam(nu,isig) = gam(nu,isig)*ry_to_ghz
     enddo  !nu
  enddo  ! isig
  !
  do isig= 1,el_ph_nsigma
     WRITE (6, 9000) deg(isig), ngauss1
     WRITE (6, 9005) dosfit(isig), effit(isig) * rytoev
     do nu=1,3*nat
        WRITE (6, 9010) nu, lamb(nu,isig), gam(nu,isig)
     enddo
  enddo
  ! Isaev: save files in suitable format for processing by lambda.x
   name=TRIM(elph_dir)// 'elph.inp_lambda.' //TRIM(int_to_char(current_iq))
                                             
  IF (ionode) THEN
     open(unit=12, file=TRIM(name), form='formatted', status='unknown', &
                                    iostat=ios)

     write(12, "(5x,3f14.6,2i6)") xq(1),xq(2),xq(3), el_ph_nsigma, 3*nat
     write(12, "(6e14.6)") (w2(i), i=1,3*nat)

     do isig= 1, el_ph_nsigma
        WRITE (12, 9000) deg(isig), ngauss1
        WRITE (12, 9005) dosfit(isig), effit(isig) * rytoev
        do nu=1,3*nat
           WRITE (12, 9010) nu, lamb(nu,isig), gam(nu,isig)
        enddo
     enddo
     close (unit=12,status='keep')
  ENDIF
  ! Isaev end
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) CALL errore('elphsum','problem opening file'//TRIM(name),1)
  deallocate (gam)
  deallocate (lamb)
  write(stdout,*)
  !
  !    Prepare interface to q2r and matdyn
  !
  call star_q (xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .TRUE. )
  !
  do isig=1,el_ph_nsigma
     name=TRIM(elph_dir)//'a2Fq2r.'// TRIM(int_to_char(50 + isig)) &
                                  //'.'//TRIM(int_to_char(current_iq))
     if (ionode) then
        iuelph = 4
        open(iuelph, file=name, STATUS = 'unknown', FORM = 'formatted', &
                     iostat=ios)
     else
        !
        ! this node doesn't write: unit 6 is redirected to /dev/null
        !
        iuelph =6
     end if
     call mp_bcast(ios, ionode_id, intra_image_comm)
     IF (ios /= 0) call errore('elphsum','opening output file '// TRIM(name),1)
     dyn22(:,:) = gf(:,:,isig)
     write(iuelph,*) deg(isig), effit(isig), dosfit(isig)
     IF ( imq == 0 ) THEN
        write(iuelph,*) 2*nq
     ELSE
        write(iuelph,*) nq
     ENDIF
     xmldyn_save=xmldyn
     xmldyn=.FALSE.
     call q2qstar_ph (dyn22, at, bg, nat, nsym, s, invs, &
          irt, rtau, nq, sxq, isq, imq, iuelph)
     xmldyn=xmldyn_save
     if (ionode) CLOSE( UNIT = iuelph, STATUS = 'KEEP' )
  enddo
  deallocate (gf)
  DEALLOCATE( deg )
  DEALLOCATE( effit )
  DEALLOCATE( dosfit )
  DEALLOCATE( xk_collect )
  IF (npool /= 1) DEALLOCATE(el_ph_mat_collect)

  call stop_clock('elphsum')

  !
9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i5,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  RETURN
END SUBROUTINE elphsum

!-----------------------------------------------------------------------
SUBROUTINE elphsum_simple
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri
  !      Rewritten by Matteo Calandra
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE constants, ONLY : pi, ry_to_cmm1, ry_to_ghz, rytoev
  USE ions_base, ONLY : nat
  USE cell_base, ONLY : at, bg
  USE symm_base, ONLY : s, irt, nsym, invs
  USE klist, ONLY : xk, nelec, nks, wk
  USE wvfct, ONLY : nbnd, et
  USE el_phon, ONLY : el_ph_mat, el_ph_nsigma, el_ph_ngauss, el_ph_sigma
  USE mp_pools,  ONLY : inter_pool_comm, npool
  USE mp_images, ONLY : intra_image_comm
  USE output, ONLY : fildyn
  USE dynmat, ONLY : dyn, w2
  USE modes, ONLY : u, nirr
  USE control_ph, only : current_iq, qplot
  USE lsda_mod, only : isk
  USE el_phon,   ONLY : done_elph, gamma_disp
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp,        ONLY: mp_sum, mp_bcast

  USE lr_symm_base, ONLY : rtau, nsymq, irotmq, minus_q
  USE qpoint, ONLY : xq, nksq, ikks, ikqs
  !
  IMPLICIT NONE
  REAL(DP), PARAMETER :: eps = 20_dp/ry_to_cmm1 ! eps = 20 cm^-1, in Ry
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, iuelph, ios, irr
  INTEGER :: nmodes
  REAL(DP) :: weight, w0g1, w0g2, w0gauss, wgauss, degauss1, dosef, &
       ef1, phase_space, lambda, gamma
  REAL(DP), EXTERNAL :: dos_ef, efermig
  character(len=80) :: filelph
  CHARACTER(len=256) ::  file_elphmat
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat), dyn_corr(3*nat,3*nat)

  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char


  DO irr=1,nirr
     IF (.NOT.done_elph(irr)) RETURN
  ENDDO

  nmodes=3*nat

  filelph=TRIM(fildyn)//'.elph.'//TRIM(int_to_char(current_iq))

  ! parallel case: only first node writes
  IF ( ionode ) THEN
     !
     iuelph = find_free_unit()
     OPEN (unit = iuelph, file = TRIM(filelph), status = 'unknown', err = &
          100, iostat = ios)
     REWIND (iuelph)
  ELSE
     iuelph = 0
     !
  END IF
100 CONTINUE
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  CALL errore ('elphsum_simple', 'opening file '//filelph, ABS (ios) )

  IF (ionode) THEN
     WRITE (iuelph, '(3f15.8,2i8)') xq, el_ph_nsigma, 3 * nat
     WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
  ENDIF
  

  ngauss1=0
  DO isig = 1, el_ph_nsigma
     !     degauss1 = 0.01 * isig
     degauss1 = el_ph_sigma * isig
     el_ph_sum(:,:) = (0.d0, 0.d0)
     phase_space = 0.d0
     !
     ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
     ! for this gaussian broadening
     !
     ! Note that the weights of k+q points must be set to zero for the
     ! following call to yield correct results
     !
      
     ef1 = efermig (et, nbnd, nks, nelec, wk, degauss1, el_ph_ngauss, 0, isk)
     dosef = dos_ef (el_ph_ngauss, degauss1, ef1, et, wk, nks, nbnd)
     ! N(Ef) is the DOS per spin, not summed over spin
     dosef = dosef / 2.d0
     !
     ! Sum over bands with gaussian weights
     !
     
     DO ik = 1, nksq
        
        !
        ! see subroutine elphel for the logic of indices
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        DO ibnd = 1, nbnd
           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
                / degauss1
           DO jbnd = 1, nbnd
              w0g2 = w0gauss ( (ef1 - et (jbnd, ikq) ) / degauss1, ngauss1) &
                   / degauss1
              ! note that wk(ikq)=wk(ikk)
              weight = wk (ikk) * w0g1 * w0g2
              DO jpert = 1, 3 * nat
                 DO ipert = 1, 3 * nat
                    el_ph_sum (ipert, jpert) = el_ph_sum (ipert, jpert)  +  weight * &
                         CONJG (el_ph_mat (jbnd, ibnd, ik, ipert) ) * &
                         el_ph_mat (jbnd, ibnd, ik, jpert)
                 ENDDO
              ENDDO
              phase_space = phase_space+weight
           ENDDO
        ENDDO
        
     ENDDO
     
     ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
     !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
     !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
     !
     ! collect contributions from all pools (sum over k-points)
     !

     call mp_sum ( el_ph_sum , inter_pool_comm )
     call mp_sum ( phase_space , inter_pool_comm )

     !
     ! symmetrize el_ph_sum(mu,nu) : it transforms as the dynamical matrix
     !
     
     CALL symdyn_munu_new (el_ph_sum, u, xq, s, invs, rtau, irt,  at, &
          bg, nsymq, nat, irotmq, minus_q)
     !
     WRITE (stdout, *)
     WRITE (stdout, 9000) degauss1, ngauss1
     WRITE (stdout, 9005) dosef, ef1 * rytoev
     WRITE (stdout, 9006) phase_space
     IF (ionode) THEN
        WRITE (iuelph, 9000) degauss1, ngauss1
        WRITE (iuelph, 9005) dosef, ef1 * rytoev
     ENDIF
     
     DO nu = 1, nmodes
        gamma = 0.d0
        DO mu = 1, 3 * nat
           DO vu = 1, 3 * nat
              gamma = gamma + DBLE (CONJG (dyn (mu, nu) ) * el_ph_sum (mu, vu)&
                   * dyn (vu, nu) )
           ENDDO
        ENDDO
        gamma = pi * gamma / 2.d0
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        IF (SQRT (ABS (w2 (nu) ) ) > eps) THEN
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lambda = gamma / pi / w2 (nu) / dosef
        ELSE
           lambda = 0.d0
        ENDIF
        WRITE (stdout, 9010) nu, lambda, gamma * ry_to_gHz
        IF (ionode) WRITE (iuelph, 9010) nu, lambda, gamma * ry_to_gHz
        IF (qplot) gamma_disp(nu,isig,current_iq) = gamma * ry_to_gHz
     ENDDO
  ENDDO
  

9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i5,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  !
  IF (ionode) CLOSE (unit = iuelph)
  RETURN
  

     
     !          call star_q(x_q(1,iq), at, bg, nsym , s , invs , nq, sxq, &
     !               isq, imq, .FALSE. )
     

END SUBROUTINE elphsum_simple
  !
  !------------------------------------------------------------
  !
  subroutine elphfil_phoebe(iq, xk_collect, ikks_collect, et_collect, &
       el_ph_mat_collect)

    USE constants, ONLY : ry_to_thz, ry_to_cmm1, amu_ry
    
    use constants, only: tpi
    use io_files, only: tmp_dir, prefix
    USE cell_base, ONLY : alat, at, bg
    USE disp, ONLY : nq1, nq2, nq3, x_q
    USE io_global, ONLY : ionode
    USE ions_base, ONLY : nat, ityp, tau, amass
    USE kinds, ONLY : DP
    USE klist, ONLY : nkstot
    USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
    USE lsda_mod, ONLY : nspin
    USE wvfct, ONLY : nbnd
    USE klist, ONLY : nelec
    USE qpoint, ONLY : nksqtot
    USE symm_base, ONLY : s, sr, nsym, t_rev, irt, time_reversal, ft
    ! ft: FractionalTranslations(3,48) in crystal coords
    USE modes, ONLY : nmodes, u
    USE dynmat, ONLY : dyn, w2
    implicit none
    ! et: KS energies in Ry
    ! xk: kpoints in cartesian coords
    ! ikks: index of computed k in the full list of ks
    ! iq: index of current irreducible qpoint
    ! el_ph: g coupling on the irreducible set of k/q
    !        for now, we are using the full set of k
    real(dp), intent(in) :: et_collect(nbnd,nkstot), xk_collect(3,nkstot)
    integer, intent(in) :: ikks_collect(nksqtot)
    integer, intent(in) :: iq
    complex(dp), intent(in) :: el_ph_mat_collect(nbnd, nbnd, nksqtot, nmodes)
    !
    logical, save :: first = .true.
    integer, allocatable, save :: k_equiv_symmetry(:), q_equiv_symmetry(:), &
         k_equiv(:), q_equiv(:), wk(:), wq(:)
    integer, save :: nk, nq, nq_irr
    real(dp), allocatable, save :: kgrid_full(:,:), qgrid_full(:,:)
    !
    integer :: iqq, iq_phonon, ik_phonon, iq_rotated, ik_rotated, &
         isym, iq_star, iq_full, iksq, i, j, n_star, unit_phoebe, &
         ii, jj, ios, kk, ll, kappa, k, ix, iy, iz, i_cart, j_cart, &
         i_eigenmode, i_pattern, i_mode
    integer, allocatable :: kappa_to_k(:)
    real(dp) :: vec(3), rotation(3,3), inv_rotation(3,3), translation(3), arg, &
         a1(3), a2(3), a3(3), k_rotated(3), q_rotated(3), q_phonon(3), &
         k_phonon(3), qdiff(3)
    real(dp), allocatable :: q_star(:,:), energies_unfolded(:,:), &
         q_star_cartesian(:,:), rotated_tau(:,:), ph_frequencies(:)
    complex(dp) :: phase
    complex(dp), allocatable :: gq_coupling(:,:,:,:,:), ph_eigenvector(:,:,:), &
         ph_star_eigenvector(:,:,:,:), tmp_el_ph_mat(:,:,:,:), &
         el_ph_mat_cartesian(:,:,:,:)
    character(len=4) :: iq_char
    character(len=200) :: file_phoebe, file_xml, line
    logical :: q_in_the_list, k_in_the_list, found
    integer, external :: find_free_unit

    ! Note: this subroutine is called n_q^irr times
    ! i.e. the number of irreducible q-points used during the phonon run
    
    if ( .not. ionode ) return

    if ( first ) then

      !-------------------------------------------------
      ! the first time this is subroutine called,
      ! we do some sanity checks that the calculation has been setup as we want
      ! we also load the symmetries of the k/q meshes
      
      first = .false.
      nk = nk1 * nk2 * nk3

      if ( nk == 0 ) then
        call errore("phoebe", &
             "We must run a nscf calculation before ph.x with automatic grid", 1)
      end if
      
      nq = nq1 * nq2 * nq3
      allocate(kgrid_full(3,nk))
      allocate(qgrid_full(3,nq))
      allocate(k_equiv(nk))
      allocate(q_equiv(nq))
      allocate(k_equiv_symmetry(nk))
      allocate(q_equiv_symmetry(nq))
      allocate(wk(nk))
      allocate(wq(nq))
      call find_irreducible_grid(nk1, nk2, nk3, k1, k2, k3, kgrid_full, &
           k_equiv, wk, k_equiv_symmetry)
      call find_irreducible_grid(nq1, nq2, nq3, 0, 0, 0, qgrid_full, &
           q_equiv, wq, q_equiv_symmetry)

      do i = 1,nk
        q_phonon = qgrid_full(:,i)
        CALL cryst_to_cart(1, q_phonon, at, -1)
        isym = q_equiv_symmetry(i)
      end do
      
      nq_irr = 0
      do i = 1,nq
        if ( q_equiv(i) == i ) then
          nq_irr = nq_irr + 1
        end if
      end do

      ! check commensurability of grids
      if ( nk1 >= nq1 ) then
        if ( mod(nk1,nq1) /= 0 ) then
          call errore("phoebe", "k and q grids are incommensurate", 1)
          end if
      else
        call errore("phoebe", "q grid larger than k grid", 1)
      end if
      if ( nk2 >= nq2 ) then
        if ( mod(nk2,nq2) /= 0 ) then
          call errore("phoebe", "k and q grids are incommensurate", 1)
        end if
      else
        call errore("phoebe", "q grid larger than k grid", 1)
      end if
      if ( nk3 >= nq3 ) then
        if ( mod(nk3,nq3) /= 0 ) then
          call errore("phoebe", "k and q grids are incommensurate", 1)
        end if
      else
        call errore("phoebe", "q grid larger than k grid", 1)
      end if

      if ( k1 /= 0 .or. k2 /= 0 .or. k3 /= 0 ) then
        call errore("phoebe", "k grid must not use Gamma offset", 1)
      end if

      if ( nk1 == 0 .or. nk2 == 0 .or. nk3 == 0 .or. &
           nq1 == 0 .or. nq2 == 0 .or. nq3 == 0 ) then
        call errore("phoebe", "Must use automatic grids in quantum espresso", 1)
      end if

      ! Here I check that we restarted from a nscf run

      file_xml = trim(tmp_dir) // "../" // trim(prefix) // ".phoebe.scf.0000.dat"
      ios = 0
      open(unit=unit_phoebe, file = TRIM(file_xml), status = 'old', iostat = ios)
      if ( ios /= 0 ) then
        call errore("phoebe", &
             "you must use the patch on pw.x calculation ahead of running ph.x", 1)
      end if
            
    end if

    !----------------------------------------------------------------------------
    ! In the first iteration, we write to file some general info
    ! (like k/q meshes, crystal, ...)
    
    if ( iq == 1 ) then

      ! Note: this unfolding doesn't work with spin
      allocate(energies_unfolded(nbnd,nk))
      energies_unfolded = 0.
      do i = 1,nkstot
        k_phonon = xk_collect(:,i) ! wavevector used in the phonon code
        CALL cryst_to_cart(1, k_phonon, at, -1) ! fold to crystal coords
        call get_point_index(iq_full, k_phonon, nk1, nk2, nk3, k1, k2, k3,&
             .false.)
        ! iq_full is the index of the vector in xk_collect inside my list
        do j = 1,nk
          if ( k_equiv(j) == iq_full ) then
            energies_unfolded(:,j) = et_collect(:,i)
          end if
        end do
      end do
        
      unit_phoebe = find_free_unit()
      file_phoebe = trim(prefix) // ".phoebe.0000.dat"
      open(unit = unit_phoebe, file = TRIM(file_phoebe), form = 'formatted', &
           access = 'sequential', status = 'replace', iostat = ios)
      write(unit_phoebe,*) "Phoebe-QE"
      write(unit_phoebe,*) nbnd, nelec, nspin
      write(unit_phoebe,*) nq1, nq2, nq3, nk1, nk2, nk3
      ! Crystal
      write(unit_phoebe,*)  alat, nat
      write(unit_phoebe,*) ((at(ii, jj)*alat, ii = 1, 3), jj = 1, 3)
      write(unit_phoebe,*) ((bg(ii, jj)*tPi/alat, ii = 1, 3), jj = 1, 3)
      write(unit_phoebe,*) (ityp(ii), ii = 1, nat)
      do jj = 1,nat
        write(unit_phoebe,*) (tau(ii, jj)*alat, ii = 1, 3)
      end do
      ! all q-points
      write(unit_phoebe,*) nq, nq_irr
      do jj = 1,nq
        write(unit_phoebe,*) (qgrid_full(ii, jj), ii = 1, 3)
      end do
      ! all k-points
      write(unit_phoebe,*) nk
      do jj = 1,nk
        write(unit_phoebe,*) (kgrid_full(ii, jj), ii = 1, 3)
      end do
      ! write KS energies at all k points
      do jj = 1,nk
        write(unit_phoebe,*) (energies_unfolded(ii, jj), ii = 1, nbnd)
      end do
      close(unit_phoebe)
    end if
    !
    !----------------------------------------------------------------------------
    ! Here we reconstruct the "star" of q-points equivalent
    ! to the current irreducible point under consideration
    !
    q_phonon = x_q(:,iq)    
    CALL cryst_to_cart(1, q_phonon, at, -1)
    
    ! iq_full is the index of q (from ph.x) in my grid of q points
    call get_point_index(iq_full, q_phonon, nq1, nq2, nq3, 0, 0, 0, .false.)
    
    ! check q point in grid
    if ( q_equiv(iq_full) /= iq_full ) then
      call errore("phoebe", "qpoint not ordered as expected 1", 1)
    end if

    qdiff(:) = q_phonon(:) -  qgrid_full(:,iq_full)
    qdiff = qdiff - nint(qdiff)

    if ( sum(qdiff**2) > 1.0d-5 ) then
      call errore("phoebe", "qpoint grid not ordered as expected 2", 1)
    end if
    
    ! number of equivalent q points in this star
    n_star = 0
    do iqq = 1,nq
      if ( q_equiv(iqq) == iq_full ) then
        n_star = n_star + 1
      end if
    end do
    allocate(q_star(3,n_star))
    j = 0
    do iqq = 1,nq
      if ( q_equiv(iqq) == iq_full ) then
        j = j + 1
        q_star(:,j) = qgrid_full(:,iqq)
      end if
    end do

    !-------------------------------------------------------------------
    ! We miss a couple of things to the el_ph_mat matrix
    ! 1) the perturbation is computed in the pattern basis
    !    we need to multiply by "u" to go in cartesian coordinates
    !    After this rotation, the 4th index represents an index
    !    over 3 cartesian coordinates and nat atoms
    ! 2) we need a further rotation to the phonon eigenvector basis
    !    so, we multiply the results by the phonon eigenvector matrix dyn.
    !    Note, I verified that dyn contains eigenvectors scaled by /sqrt(mass)
    !    and that has dimension dyn(3*nat,nPhModes)

    ! el_ph_mat_collect(nbnd, nbnd, nksqtot, nmodes)
    ! the el=ph coupling was computed in the pattern representation
    ! Now I want to rotate it to cartesian coordinates
    allocate(tmp_el_ph_mat(nbnd, nbnd, nksqtot, nmodes))
    tmp_el_ph_mat = cmplx(0.,0.,kind=dp)
    do i_pattern = 1,nmodes
      do i_mode = 1,nmodes
        tmp_el_ph_mat(:,:,:,i_mode) = tmp_el_ph_mat(:,:,:,i_mode) + &
             conjg(u(i_mode,i_pattern)) * el_ph_mat_collect(:,:,:,i_pattern)
      end do
    end do
    allocate(el_ph_mat_cartesian(nbnd, nbnd, nksqtot, nmodes))
    el_ph_mat_cartesian = cmplx(0.,0.,kind=dp)
    do i_pattern = 1,nmodes
      do i_mode = 1,nmodes
        el_ph_mat_cartesian(:,:,:,i_mode) = el_ph_mat_cartesian(:,:,:,i_mode) + &
             conjg(dyn(i_pattern,i_mode)) * tmp_el_ph_mat(:,:,:,i_pattern)
      end do
    end do
    deallocate(tmp_el_ph_mat)
    
    !----------------------------------------------------------------------------
    ! Here we start the coupling unfolding
    ! we use g(Sk,q) = g(k,S^{-1}q) to unfold the coupling
    ! which is equivalent to g(Sk,Sq) = g(k,q)
    
    allocate(gq_coupling(nbnd, nbnd, nk, nmodes, n_star))
    gq_coupling = cmplx(0., 0., kind=dp)

    ! here we actually apply g(Sk,Sq^irr) = g(k,q^irr)
    
    do iksq = 1,nksqtot ! loop over the k's computed by phonon.x
      ! (xk_collect(3, nkstot) k point coordinates
      ! ikks_collect(nksqtot) ! indeces of k vectors in the loop nksqtot
      ik_phonon = ikks_collect(iksq) ! index of k in xk_collect list
      k_phonon = xk_collect(:,ik_phonon) ! k wavevector in cartesian coordinates
      CALL cryst_to_cart(1, k_phonon, at, -1) ! in crystal coords
      !
      do isym = 1,nsym
        k_rotated = matmul(s(:,:,isym), k_phonon)
        q_rotated = matmul(s(:,:,isym), q_phonon) ! q_phonon is the irred q-point
        k_rotated(:) = k_rotated(:) - nint( k_rotated(:) )
        q_rotated(:) = q_rotated(:) - nint( q_rotated(:) )
        if ( t_rev(isym) == 1 ) then
          k_rotated = - k_rotated
          q_rotated = - q_rotated
        end if
        call is_point_in_grid(k_in_the_list, k_rotated, nk1, nk2, nk3, 0, 0, 0, .false.)
        call is_point_in_grid(q_in_the_list, q_rotated, nq1, nq2, nq3, 0, 0, 0, .false.)
        !
        if ( k_in_the_list .and. q_in_the_list ) then
          ! index of point in the full grid
          call get_point_index(ik_rotated, k_rotated, nk1, nk2, nk3, 0, 0, 0, .false.)
          call get_point_index(iq_rotated, q_rotated, nq1, nq2, nq3, 0, 0, 0, .false.)
          
          iq_star = 0
          do iqq = 1,n_star
            if ( sum((qgrid_full(:,iq_rotated)-q_star(:,iqq))**2) < 1.0d-5 ) then
              iq_star = iqq
              exit
            end if
          end do
          if ( iq_star == 0 ) call errore("phoebe", "q not found in star",1)
          
          gq_coupling(:,:,ik_rotated,:,iq_star) = el_ph_mat_cartesian(:,:,iksq,:)
        end if
        !
        if ( time_reversal ) then
          call is_point_in_grid(k_in_the_list, k_rotated, nk1, nk2, nk3, &
               k1, k2, k3, .true.)
          call is_point_in_grid(q_in_the_list, q_rotated, nq1, nq2, nq3, &
               0, 0, 0, .true.)
          IF ( k_in_the_list .and. q_in_the_list ) THEN
            call get_point_index(ik_rotated, k_rotated, nk1, nk2, nk3, &
                 k1, k2, k3, .true.)
            call get_point_index(iq_rotated, q_rotated, nq1, nq2, nq3, &
                 0, 0, 0, .true.)
            
            iq_star = 0
            do iqq = 1,n_star
              if (sum((qgrid_full(:,iq_rotated)-q_star(:,iqq))**2) < 1.0d-5) then
                iq_star = iqq
                exit
              end if
            end do
            if ( iq_star == 0 ) call errore("phoebe", "q not found in star",1)
            
            gq_coupling(:, :, ik_rotated, :, iq_star) = &
                 el_ph_mat_cartesian(:, :, iksq, :)
          end if
        end if
        !
      end do ! end nsym
    end do

    !------------ Rotate the phonon eigenvectors -------------------------------
    ! I must do this, so that we can keep their phase information,
    ! needed to Fourier transform them to real space
    ! Eq. 2.33, Maradudin & Vosko, Rev. Mod. Phys. (1968)
    ! https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.40.1

    ! Note: a1 is in alat units, cartesian cords
    ! tau is in alat units, cartesian coords
    
    allocate(q_star_cartesian(3,n_star))
    q_star_cartesian = q_star
    do iq_star = 1,n_star
      call cryst_to_cart(1, q_star_cartesian(:,iq_star), bg, 1)
    end do
    
    a1 = at(:,1)
    a2 = at(:,2)
    a3 = at(:,3)

    ! Note: I checked that dyn contains the eigenvectors
    ! normalized by sqrt of atomic mass
    ! to verify that, run:
    ! allocate(eigvec(3*nat,3*nat))
    ! eigvec = dyn * sqrt(amass(1) * amu_ry)
    ! dyn = cmplx(0.,0.,kind=dp)
    ! do ii = 1,3*nat
    !   dyn(ii,ii) = sqrt(w2(ii))
    ! end do
    ! dyn = matmul(eigvec, dyn)
    ! dyn = matmul(dyn,transpose(eigvec))
    ! allocate(freq(3*nat))
    ! CALL cdiagh (3*nat, dyn, 3*nat, freq, eigvec)
    ! print*, freq*ry_to_cmm1
    ! print*, w2*ry_to_cmm1
    ! print*, eigvec(:,1) ! this gives the same pattern of *.dyn* files
    ! deallocate(eigvec)

    allocate(ph_eigenvector(3,nat,3*nat))
    ph_eigenvector = cmplx(0.,0.,kind=dp)
    do k = 1,nat
      do i_cart = 1,3
        ii = 3 * (k-1) + i_cart
        do jj = 1,3*nat
          ph_eigenvector(i_cart,k,jj) = dyn(ii,jj)
        end do
      end do
    end do
    allocate(ph_star_eigenvector(3,nat,3*nat,n_star))
    ph_star_eigenvector = cmplx(0.,0.,kind=dp)
    ph_star_eigenvector(:,:,:,1) = ph_eigenvector
    
    !--------------------------------------------
    
    if ( n_star > 1 ) then ! if not, nothing to rotate
      allocate(kappa_to_k(nat))
      allocate(rotated_tau(3,nat))
      do iq_star = 2,n_star
        ! step 1: identify the symmetry S such that q^rotated = S q^irreducible

        ! iqq is the index of iq_star in the full grid
        call get_point_index(iqq, q_star(:,iq_star), nq1, nq2, nq3, 0, 0, 0, .false.)
        ! conveniently, I computed already what is the symmetry to use here
        isym = q_equiv_symmetry(iqq)        
        rotation = sr(:,:,isym)
        inv_rotation = transpose(rotation)
        translation(:) = ft(:,isym)
        ! transform in cartesian coordinates
        call cryst_to_cart(1,translation(:),at,+1)

        ! step 2: build the mapping kappa_to_k
        kappa_to_k = 0
        do kappa = 1,nat
          ! note how the translation must be applied before the rotation
          rotated_tau(:,kappa) = matmul(rotation,tau(:,kappa) + translation)
          k = 0
          do ii = 1,nat
            do ix = -3,3
              do iy = -3,3
                do iz = -3,3
                  vec(:) = rotated_tau(:,kappa) &
                       + ix * a1(:) + iy * a2(:) + iz * a3(:)
                  if ( sum( ( vec - tau(:,ii))**2 ) < 1.0e-5 ) then
                    k = ii
                  end if
                end do
              end do
            end do
          end do
          if ( k == 0 ) then
            call errore("phoebe", "tau not folded", 1)
          end if
          kappa_to_k(kappa) = k
        end do
        
        ! step 3: rotation
        
        do kappa = 1,nat
          k = kappa_to_k(kappa)
          vec(:) = matmul(inv_rotation,tau(:,kappa)-translation) - tau(:,k)
          arg = tpi * sum( q_star_cartesian(:,1) * vec(:) )
          phase = cmplx(cos(arg),sin(arg),kind=dp)          
          do i_eigenmode = 1,nmodes
            do i_cart = 1,3
              do j_cart = 1,3
                jj = 3 * (kappa-1) + j_cart
                ph_star_eigenvector(i_cart,k,i_eigenmode,iq_star) =        &
                     ph_star_eigenvector(i_cart,k,i_eigenmode,iq_star)     &
                     + rotation(i_cart,j_cart) * phase * dyn(jj,i_eigenmode)
              end do
            end do
          end do
        end do
      end do
      deallocate(kappa_to_k)
      deallocate(rotated_tau)
    end if

    ! ph frequencies
    allocate(ph_frequencies(3*nat))
    do ii = 1,nmodes
      if ( w2(ii) >= 0.d0 ) then
        ph_frequencies(ii) = sqrt(w2(ii))
      else
        ph_frequencies(ii) = - sqrt(-w2(ii))
      end if
    end do

    !---------------------------------------------
    !---------- Dump results to file -------------
    !---------------------------------------------

    unit_phoebe = find_free_unit()
    write(iq_char,"(I4.4)") iq
    file_phoebe = trim(prefix) // ".phoebe." // trim(iq_char) // ".dat"
    open(unit = unit_phoebe, file = TRIM(file_phoebe), form = 'formatted', &
         access = 'sequential', status = 'replace', iostat = ios)
    write(unit_phoebe,*) n_star
    do iqq = 1,n_star
      write(unit_phoebe,*) (q_star(ii,iqq), ii = 1,3) ! crystal coords
    end do
    do iqq = 1,n_star
      write(unit_phoebe,*) (q_star_cartesian(ii,iqq), ii = 1,3)
    end do
    write(unit_phoebe,*) (ph_frequencies(ii), ii = 1,nmodes) ! in rydbergs
    do iqq = 1,n_star
      do jj = 1,nmodes
        do k = 1,nat
          do i_cart = 1,3
            write(unit_phoebe,"(2ES24.16)") ph_star_eigenvector(i_cart,k,jj,iqq)
          end do
        end do
      end do
    end do
    write(unit_phoebe,*) ""
    write(unit_phoebe,"(2ES24.16)") (((((gq_coupling(ii, jj, kk, ll, iqq), &
         ii=1,nbnd), jj=1,nbnd), kk=1,nk), ll=1,nmodes), iqq=1,n_star)
    close(unit_phoebe)
    
    deallocate(gq_coupling, q_star_cartesian, ph_star_eigenvector)
    deallocate(ph_eigenvector, el_ph_mat_cartesian, ph_frequencies)
    
    return
  end subroutine elphfil_phoebe
  !
  !-----------------------------------------------------------------------
  !
SUBROUTINE elphfil_epa(iq)
  !-----------------------------------------------------------------------
  !
  !      Writes electron-phonon matrix elements to a file
  !      which is subsequently processed by the epa code
  !      Original routine written by Georgy Samsonidze
  !
  !-----------------------------------------------------------------------
  USE cell_base, ONLY : ibrav, alat, omega, tpiba, at, bg
  USE disp, ONLY : nq1, nq2, nq3, nqs, x_q, wq, lgamma_iq
  USE dynmat, ONLY : dyn, w2
  USE el_phon, ONLY : el_ph_mat, done_elph
  USE fft_base, ONLY : dfftp, dffts, dfftb
  USE gvect, ONLY : ngm_g, ecutrho
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, nsp, atm, ityp, tau
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, nelec, nks, nkstot, ngk
  USE lsda_mod, ONLY : nspin, isk
  USE modes, ONLY : nirr, nmodes, npert, npertx, u, t, tmq, &
       name_rap_mode, num_rap_mode
  USE lr_symm_base, ONLY : irgq, nsymq, irotmq, rtau, gi, gimq, &
       minus_q, invsymq
  USE mp, ONLY : mp_bcast, mp_sum
  USE mp_images, ONLY : intra_image_comm
  USE mp_pools, ONLY : npool, intra_pool_comm
  USE qpoint, ONLY : nksq, nksqtot, ikks, ikqs, eigqts
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, invs, ft, nrot, nsym, nsym_ns, &
       nsym_na, ft, sr, sname, t_rev, irt, time_reversal, &
       invsym, nofrac, allfrac, nosym, nosym_evc, no_t_rev
  USE wvfct, ONLY : nbnd, et, wg
  USE gvecw, ONLY : ecutwfc
  USE io_files, ONLY : prefix

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iq

  INTEGER :: iuelph, ios, irr, ii, jj, kk, ll
  character :: cdate*9, ctime*9, sdate*32, stime*32, &
       stitle*32, myaccess*10, mystatus*7
  CHARACTER(LEN=80) :: filelph

  REAL(DP), ALLOCATABLE :: xk_collect(:,:), wk_collect(:)
  REAL(DP), ALLOCATABLE :: et_collect(:,:), wg_collect(:,:)
  INTEGER, ALLOCATABLE :: ngk_collect(:)
  INTEGER, ALLOCATABLE :: ikks_collect(:), ikqs_collect(:)
  COMPLEX(DP), ALLOCATABLE :: el_ph_mat_collect(:,:,:,:)
  INTEGER :: ftau(3,48)
  INTEGER, EXTERNAL :: find_free_unit, atomic_number

  filelph = TRIM(prefix) // '.epa.k'

  DO irr = 1, nirr
     IF (.NOT. done_elph(irr)) RETURN
  ENDDO

  IF (iq .EQ. 1) THEN
     myaccess = 'sequential'
     mystatus = 'replace'
  ELSE
     myaccess = 'append'
     mystatus = 'old'
  ENDIF
  IF (ionode) THEN
     iuelph = find_free_unit()
     OPEN(unit = iuelph, file = TRIM(filelph), form = 'unformatted', &
          access = myaccess, status = mystatus, iostat = ios)
  ELSE
     iuelph = 0
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('elphfil_epa', 'opening file ' // filelph, ABS(ios))

  IF (iq .EQ. 1) THEN
     CALL date_and_tim(cdate, ctime)
     WRITE(sdate, '(A2,"-",A3,"-",A4,21X)') cdate(1:2), cdate(3:5), cdate(6:9)
     WRITE(stime, '(A8,24X)') ctime(1:8)
     WRITE(stitle, '("EPA-Complex",21X)')
     CALL cryst_to_cart(nqs, x_q, at, -1)
     ! write header
     IF (ionode) THEN
        WRITE(iuelph) stitle, sdate, stime
        WRITE(iuelph) ibrav, nat, nsp, nrot, nsym, nsym_ns, nsym_na, &
             ngm_g, nspin, nbnd, nmodes, nqs
        WRITE(iuelph) nq1, nq2, nq3, nk1, nk2, nk3, k1, k2, k3
        WRITE(iuelph) time_reversal, invsym, nofrac, allfrac, nosym, &
             nosym_evc, no_t_rev
        WRITE(iuelph) alat, omega, tpiba, nelec, ecutrho, ecutwfc
        WRITE(iuelph) dfftp%nr1, dfftp%nr2, dfftp%nr3
        WRITE(iuelph) dffts%nr1, dffts%nr2, dffts%nr3
        WRITE(iuelph) dfftb%nr1, dfftb%nr2, dfftb%nr3
        WRITE(iuelph) ((at(ii, jj), ii = 1, 3), jj = 1, 3)
        WRITE(iuelph) ((bg(ii, jj), ii = 1, 3), jj = 1, 3)
        WRITE(iuelph) (atomic_number(atm(ii)), ii = 1, nsp)
        WRITE(iuelph) (ityp(ii), ii = 1, nat)
        WRITE(iuelph) ((tau(ii, jj), ii = 1, 3), jj = 1, nat)
        WRITE(iuelph) ((x_q(ii, jj), ii = 1, 3), jj = 1, nqs)
        WRITE(iuelph) (wq(ii), ii = 1, nqs)
        WRITE(iuelph) (lgamma_iq(ii), ii = 1, nqs)
     ENDIF
     CALL cryst_to_cart(nqs, x_q, bg, 1)
  ENDIF

  ! collect data for current q-point
  ALLOCATE(xk_collect(3, nkstot))
  ALLOCATE(wk_collect(nkstot))
  ALLOCATE(et_collect(nbnd, nkstot))
  ALLOCATE(wg_collect(nbnd, nkstot))
  ALLOCATE(ngk_collect(nkstot))
  ALLOCATE(ikks_collect(nksqtot))
  ALLOCATE(ikqs_collect(nksqtot))
  ALLOCATE(el_ph_mat_collect(nbnd, nbnd, nksqtot, nmodes))
  IF (npool > 1) THEN
     CALL poolcollect(3, nks, xk, nkstot, xk_collect)
     CALL poolcollect(1, nks, wk, nkstot, wk_collect)
     CALL poolcollect(nbnd, nks, et, nkstot, et_collect)
     CALL poolcollect(nbnd, nks, wg, nkstot, wg_collect)
     CALL ipoolcollect(1, nks, ngk, nkstot, ngk_collect)
     CALL jpoolcollect(1, nksq, ikks, nksqtot, ikks_collect)
     CALL jpoolcollect(1, nksq, ikqs, nksqtot, ikqs_collect)
     CALL el_ph_collect(nmodes, el_ph_mat, el_ph_mat_collect, nksqtot, nksq)
  ELSE
     xk_collect(1:3, 1:nks) = xk(1:3, 1:nks)
     wk_collect(1:nks) = wk(1:nks)
     et_collect(1:nbnd, 1:nks) = et(1:nbnd, 1:nks)
     wg_collect(1:nbnd, 1:nks) = wg(1:nbnd, 1:nks)
     ngk_collect(1:nks) = ngk(1:nks)
     ikks_collect(1:nksq) = ikks(1:nksq)
     ikqs_collect(1:nksq) = ikqs(1:nksq)
     el_ph_mat_collect(1:nbnd, 1:nbnd, 1:nksq, 1:nmodes) = &
          el_ph_mat(1:nbnd, 1:nbnd, 1:nksq, 1:nmodes)
  ENDIF
  CALL cryst_to_cart(nkstot, xk_collect, at, -1)
  ! write data for current q-point
  IF (ionode) THEN
     WRITE(iuelph) nsymq, irotmq, nirr, npertx, nkstot, nksqtot
     WRITE(iuelph) minus_q, invsymq
     WRITE(iuelph) (irgq(ii), ii = 1, 48)
     WRITE(iuelph) (npert(ii), ii = 1, nmodes)
     WRITE(iuelph) (((rtau(ii, jj, kk), ii = 1, 3), jj = 1, 48), &
          kk = 1, nat)
     WRITE(iuelph) ((gi(ii, jj), ii = 1, 3), jj = 1, 48)
     WRITE(iuelph) (gimq(ii), ii = 1, 3)
     WRITE(iuelph) ((u(ii, jj), ii = 1, nmodes), jj = 1, nmodes)
     WRITE(iuelph) ((((t(ii, jj, kk, ll), ii = 1, npertx), &
          jj = 1, npertx), kk = 1, 48), ll = 1, nmodes)
     WRITE(iuelph) (((tmq(ii, jj, kk), ii = 1, npertx), &
          jj = 1, npertx), kk = 1, nmodes)
     WRITE(iuelph) (name_rap_mode(ii), ii = 1, nmodes)
     WRITE(iuelph) (num_rap_mode(ii), ii = 1, nmodes)
     WRITE(iuelph) (((s(ii, jj, kk), ii = 1, 3), jj = 1, 3), kk = 1, 48)
     WRITE(iuelph) (invs(ii), ii = 1, 48)
     ! FIXME: should disappear
     ftau(1,1:48) = NINT(ft(1,1:48)*dfftp%nr1)
     ftau(2,1:48) = NINT(ft(2,1:48)*dfftp%nr2)
     ftau(3,1:48) = NINT(ft(3,1:48)*dfftp%nr3)
     WRITE(iuelph) ((ftau(ii, jj), ii = 1, 3), jj = 1, 48)
     ! end FIXME
     WRITE(iuelph) ((ft(ii, jj), ii = 1, 3), jj = 1, 48)
     WRITE(iuelph) (((sr(ii, jj, kk), ii = 1, 3), jj = 1, 3), kk = 1, 48)
     WRITE(iuelph) (sname(ii), ii = 1, 48)
     WRITE(iuelph) (t_rev(ii), ii = 1, 48)
     WRITE(iuelph) ((irt(ii, jj), ii = 1, 48), jj = 1, nat)
     WRITE(iuelph) ((xk_collect(ii, jj), ii = 1, 3), jj = 1, nkstot)
     WRITE(iuelph) (wk_collect(ii), ii = 1, nkstot)
     WRITE(iuelph) ((et_collect(ii, jj), ii = 1, nbnd), jj = 1, nkstot)
     WRITE(iuelph) ((wg_collect(ii, jj), ii = 1, nbnd), jj = 1, nkstot)
     WRITE(iuelph) (isk(ii), ii = 1, nkstot)
     WRITE(iuelph) (ngk_collect(ii), ii = 1, nkstot)
     WRITE(iuelph) (ikks_collect(ii), ii = 1, nksqtot)
     WRITE(iuelph) (ikqs_collect(ii), ii = 1, nksqtot)
     WRITE(iuelph) (eigqts(ii), ii = 1, nat)
     WRITE(iuelph) (w2(ii), ii = 1, nmodes)
     WRITE(iuelph) ((dyn(ii, jj), ii = 1, nmodes), jj = 1, nmodes)
     WRITE(iuelph) ((((el_ph_mat_collect(ii, jj, kk, ll), ii = 1, nbnd), &
          jj = 1, nbnd), kk = 1, nksqtot), ll = 1, nmodes)
     CLOSE (unit = iuelph, status = 'keep')
  ENDIF
  CALL cryst_to_cart(nkstot, xk_collect, bg, 1)

  call elphfil_phoebe(iq, xk_collect, ikks_collect, et_collect, el_ph_mat_collect)

  DEALLOCATE(xk_collect)
  DEALLOCATE(wk_collect)
  DEALLOCATE(et_collect)
  DEALLOCATE(wg_collect)
  DEALLOCATE(ngk_collect)
  DEALLOCATE(ikks_collect)
  DEALLOCATE(ikqs_collect)
  DEALLOCATE(el_ph_mat_collect)

  RETURN

END SUBROUTINE elphfil_epa
   
!----------------------------------------------------------------------------
SUBROUTINE ipoolcollect( length, nks, f_in, nkstot, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... as poolcollect, for an integer vector
  !
  USE mp_pools,  ONLY : my_pool_id, npool, kunit, &
                        inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length, nks, nkstot
  ! first dimension of arrays
  ! number of k-points per pool
  ! total number of k-points
  INTEGER, INTENT(IN)  :: f_in (length,nks)
  ! pool-distributed function
  INTEGER, INTENT(OUT) :: f_out(length,nkstot)
  ! pool-collected function
  !
  INTEGER :: nbase, rest, nks1
  !
  nks1    = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1.ne.nks) &
     call errore('ipoolcollect','inconsistent number of k-points',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  ! copy the original points in the correct position of the list
  !
  f_out=0
  f_out(:,nbase+1:nbase+nks) = f_in(:,1:nks)
  !
  CALL mp_sum( f_out, inter_pool_comm )
  !
  RETURN
  !
END SUBROUTINE ipoolcollect

!----------------------------------------------------------------------------
SUBROUTINE jpoolcollect( length, nks, f_in, nkstot, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... as ipoolcollect, without kunit and with an index shift
  !
  USE mp_pools,  ONLY : my_pool_id, npool, kunit, &
                        inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length, nks, nkstot
  ! first dimension of arrays
  ! number of k-points per pool
  ! total number of k-points
  INTEGER, INTENT(IN)  :: f_in (length,nks)
  ! pool-distributed function
  INTEGER, INTENT(OUT) :: f_out(length,nkstot)
  ! pool-collected function
  !
  INTEGER :: nbase, rest, nks1
  !
  nks1    = ( nkstot / npool )
  !
  rest = ( nkstot - nks1 * npool )
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + 1
  !
  IF (nks1.ne.nks) &
     call errore('jpoolcollect','inconsistent number of k-points',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest
  !
  ! copy the original points in the correct position of the list
  !
  f_out=0
  f_out(:,nbase+1:nbase+nks) = f_in(:,1:nks) + nbase * kunit
  !
  CALL mp_sum( f_out, inter_pool_comm )
  !
  RETURN
  !
END SUBROUTINE jpoolcollect
   
!-----------------------------------------------------------------------
FUNCTION dos_ef (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE mp_pools, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_sum
  IMPLICIT NONE
  REAL(DP) :: dos_ef
  INTEGER :: ngauss, nbnd, nks
  REAL(DP) :: et (nbnd, nks), wk (nks), ef, degauss
  !
  INTEGER :: ik, ibnd
  REAL(DP), EXTERNAL :: w0gauss
  !
  !     Compute DOS at E_F (states per Ry per unit cell)
  !
  dos_ef = 0.0d0
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        dos_ef = dos_ef + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     ENDDO
  ENDDO
  !
  !    Collects partial sums on k-points from all pools
  !
  CALL mp_sum ( dos_ef, inter_pool_comm )
  !
  RETURN
END FUNCTION dos_ef

!a2F
subroutine lint ( nsym, s, minus_q, at, bg, npk, k1,k2,k3, &
     nk1,nk2,nk3, nks, xk, kunit, nkBZ, eqBZ, sBZ)
  !-----------------------------------------------------------------------
  !
  ! Find which k-points of a uniform grid are in the IBZ
  !
  use kinds, only : DP
  implicit none
  integer, intent (IN) :: nks, nsym, s(3,3,48), npk, k1, k2, k3, &
       nk1, nk2, nk3, kunit, nkBZ
  logical, intent (IN) :: minus_q
  real(kind=DP), intent(IN):: at(3,3), bg(3,3), xk(3,npk)
  integer, INTENT(OUT) :: eqBZ(nkBZ), sBZ(nkBZ)
  !
  real(kind=DP) :: xkr(3), deltap(3), deltam(3)
  real(kind=DP), parameter:: eps=1.0d-5
  real(kind=DP), allocatable :: xkg(:,:), xp(:,:)
  integer ::  i,j,k, ns, n, nk
  integer :: nkh
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  allocate (xkg( 3,nkBZ))
  !
  if(kunit < 1 .or. kunit > 2) call errore('lint','bad kunit value',kunit)
  !
  ! kunit=2: get only "true" k points, not k+q points, from the list
  !
  nkh = nks/kunit
  allocate (xp(3,nkh))
  if (kunit == 1) then
     xp(:,1:nkh) = xk(:,1:nkh)
  else
     do j=1,nkh
        xp(:,j) = xk(:,2*j-1)
     enddo
  end if
  do i=1,nk1
     do j=1,nk2
        do k=1,nk3
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        end do
     end do
  end do

  call cryst_to_cart (nkh,xp,at,-1)

  do nk=1,nkBZ
     do n=1,nkh
        do ns=1,nsym
           do i=1,3
              xkr(i) = s(i,1,ns) * xp(1,n) + &
                       s(i,2,ns) * xp(2,n) + &
                       s(i,3,ns) * xp(3,n)
           end do
           do i=1,3
              deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
              deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
           end do
           if ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) < eps .or. ( minus_q .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) < eps ) ) then
              eqBZ(nk) = n
              sBZ(nk) = ns
              go to 15
           end if
        end do
     end do
     call errore('lint','cannot locate  k point  xk',nk)
15   continue
  end do

  do n=1,nkh
     do nk=1,nkBZ
        if (eqBZ(nk) == n) go to 20
     end do
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     call errore('lint','cannot remap grid on k-point list',n)
20   continue
  end do

  deallocate(xkg)
  deallocate(xp)

  return
end subroutine lint
