! Phoebe

subroutine get_point_index(n,k_rotated,nk1,nk2,nk3,k1,k2,k3, time_reversal)
  ! given a point in crystal coordinates,
  ! finds the index of the point in the full grid
  use kinds, only: dp
  implicit none
  integer, intent(out) :: n
  integer, intent(in) :: nk1, nk2, nk3, k1, k2, k3
  integer :: i, j, k
  logical, intent(in) :: time_reversal
  real(dp), intent(in) :: k_rotated(3)
  real(dp) :: fac
  if ( time_reversal ) then
    fac = - 1.
  else
    fac = 1.
  end if
  i = mod ( nint ( fac * k_rotated(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
  j = mod ( nint ( fac * k_rotated(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
  k = mod ( nint ( fac * k_rotated(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
  n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
  return
end subroutine get_point_index
!
!----------------------------------------------------
!
subroutine is_point_in_grid(is_k_in_list, k_rotated, &
     nk1, nk2, nk3, k1, k2, k3, time_reversal)
  ! given a point in crystal coordinates,
  ! checks that it falls on the Monkhorst-Pack mesh
  use kinds, only: dp
  implicit none
  logical, intent(in) :: time_reversal
  logical, intent(out) :: is_k_in_list
  integer, intent(in) :: nk1,nk2,nk3,k1,k2,k3
  real(dp), intent(in) :: k_rotated(3)
  real(dp) :: xx, yy, zz, fac
  real(dp), parameter :: eps = 1.0d-5
  
  if ( time_reversal ) then
    fac = - 1.
  else
    fac = 1.
  end if
  xx = fac * k_rotated(1)*nk1 - 0.5d0*k1
  yy = fac * k_rotated(2)*nk2 - 0.5d0*k2
  zz = fac * k_rotated(3)*nk3 - 0.5d0*k3
  is_k_in_list = abs(xx-nint(xx))<=eps .and. &
       abs(yy-nint(yy))<=eps .and. &
       abs(zz-nint(zz))<=eps
  return
end subroutine is_point_in_grid
!
!-----------------------------------------------------------------------------
!
subroutine find_irreducible_grid(nk1, nk2, nk3, k1, k2, k3, xkg, equiv, &
     wkk, equiv_symmetry)
  ! Returns:
  ! 1) the full grid of k-points
  ! 2) the weight of k-points on the full grid
  !    although the reducible points should be discarded here
  ! 3) equiv_symmetry: the index of the rotation such that S*k^red = k^irr
  use kinds, only: dp
  use symm_base, only: s, nsym, t_rev, irt, time_reversal, sname
  implicit none
  integer, intent(in) :: nk1, nk2, nk3, k1, k2, k3
  integer, intent(out) :: equiv(nk1*nk2*nk3), wkk(nk1*nk2*nk3), &
       equiv_symmetry(nk1*nk2*nk3)
  real(dp), intent(out) :: xkg(3,nk1*nk2*nk3)
  !
  integer :: i, j, k, n, nk, nkr, ns, ik
  real(dp) :: xkr(3)
  logical :: in_the_list
  
  nkr = nk1*nk2*nk3
  do i=1,nk1
    do j=1,nk2
      do k=1,nk3
        !  this is nothing but consecutive ordering
        n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
        !  xkg are the components of the complete grid in crystal axis
        xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
        xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
        xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
      end do
    end do
  end do
  
  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
  
  do nk=1,nkr
    equiv(nk)=nk
    wkk(nk) = 1
    xkg(:,nk) = xkg(:,nk) - nint(xkg(:,nk))
  end do
    
  do nk=1,nkr
    !  check if this k-point has already been found equivalent to another
    if (equiv(nk) == nk) THEN
      wkk(nk) = 1
      !  check if there are equivalent k-point to this in the list
      !  (excepted those previously found to be equivalent to another)
      !  check both k and -k
      do ns=1,nsym
        do i=1,3
          xkr(i) = s(i,1,ns) * xkg(1,nk) &
               + s(i,2,ns) * xkg(2,nk) &
               + s(i,3,ns) * xkg(3,nk)
          xkr(i) = xkr(i) - nint( xkr(i) )
        end do
        if ( t_rev(ns) == 1 ) xkr = - xkr
        
        call is_point_in_grid(in_the_list, xkr, nk1, nk2, nk3, k1, k2, k3,&
             .false.)
        
        if (in_the_list) THEN
          
          call get_point_index(n, xkr, nk1, nk2, nk3, k1, k2, k3, .false.)
          
          if ( n>nk .and. equiv(n)==n ) then
            equiv(n) = nk
            wkk(nk) = wkk(nk) + 1
          else
            IF (equiv(n)/=nk .or. n<nk ) CALL errore('phoeve_kpoint_grid', &
                 'something wrong in the checking algorithm',1)
          end if
        end if
        if ( time_reversal ) then
          
          call is_point_in_grid(in_the_list, xkr, nk1, nk2, nk3, k1, k2, k3,&
               .true.)
          
          if (in_the_list) then
            
            call get_point_index(n, xkr, nk1, nk2, nk3, k1, k2, k3, .true.)
            
            if (n>nk .and. equiv(n)==n) then
              equiv(n) = nk
              wkk(nk) = wkk(nk) + 1
            else
              if (equiv(n)/=nk.or.n<nk) CALL errore('kpoint_grid', &
                   'something wrong in the checking algorithm',2)
            end if
          end if
        end if
      end do
    end if
  end do

  ! equiv now tells us if the point is irreducible
  ! now, find the rotation that maps the irreducible k-point to red point

  equiv_symmetry = 1 ! default at identity matrix

  do ik=1,nkr
    !  check if this k-point has already been found equivalent to another
    if (equiv(ik) == ik) then
      cycle ! identity matrix is the right one
    end if
    !
    !  check if there are equivalent k-point to this in the list
    !  (excepted those previously found to be equivalent to another)
    !  check both k and -k
    do ns=1,nsym
      xkr(:) = matmul(s(:,:,ns), xkg(:,ik))
      xkr(:) = xkr(:) - nint( xkr(:) )
      ! if ( t_rev(ns) == 1 ) xkr = - xkr ! this is for magnetisim
      
      call is_point_in_grid(in_the_list, xkr, nk1, nk2, nk3, k1, k2, k3,&
           .false.)
        
      if (in_the_list) THEN
        call get_point_index(n, xkr, nk1, nk2, nk3, k1, k2, k3, .false.)
        if ( equiv(ik) == n ) then
          equiv_symmetry(ik) = ns
          cycle
        end if
      end if
    end do
    
    if ( equiv_symmetry(ik) == 1 ) then
      call errore("phoebe", "Failed to find rotation",1)
    end if
    
  end do

  return 
end subroutine find_irreducible_grid
!
!---------------------------------------------------------
!
subroutine set_wavefunction_gauge(ik)
  ! Subroutine to fix gauge of the wavefunction and satisfy
  ! the symmetry properties of the wavefunction.
  !
  ! A few notes:
  ! The rotation that is performed is described in the Phoebe documentation page.
  ! As is, the code doesn't work with magnetism (spin-orbit, or LSDA)
  !
  ! Technical notes:
  ! evc(ig,ib): plane wave coefficients of G-vector ig and band ib.
  !             MPI-distributed over G-vectors
  ! ig_l2g: maps the local g vector to the global G-vector
  ! igk_k(ig,ik): the G-vectors of the wavefunction are ordered differently
  !               at every k-point. This array maps the local G-vector order to the
  !               global G-vector order
  ! g: list of G-vectors, ordered by magnitude of |G|^2, in cartesian coords
  ! ngm_g: number of G-vectors in the global list
  ! ngk(ik): MPI-local number of plane wave coefficients for evc
  !          note that, even in serial, ngm_g>>ngk(ik), for reasons
  
  use gvect, only: gstart, g, ig_l2g, ngm, ngm_g, mill_g
  use wavefunctions, only: evc ! plane wave coefficients, evc(npwx*npol,nbnd)
  use wvfct, only: et, & ! eigenvalues of the Hamiltonian et(nbnd,nkstot)
       nbnd, & ! number of bands
       npwx, & ! maximum number of PW for wavefunctions
       nbndx, & ! max number of bands used in iterative diag
       npw ! number of plane waves
  use constants, only: pi, ryToEv, tpi ! greek-pi
  use kinds, only: dp
  use parallel_include
  use input_parameters, only: calculation
  use mp_pools, only: intra_pool_comm, me_pool, root_pool, &
       nproc_pool, my_pool_id
  use mp, only: mp_bcast, mp_sum
  use mp_world, only: mpime
  use klist, only: igk_k, ngk, xk
  use cell_base, only: tpiba, bg, at
  use io_files, only: prefix, tmp_dir
  use start_k, only: nk1, nk2, nk3, k1, k2, k3
  USE symm_base, only: s, sr, nsym, t_rev, irt, time_reversal, ft
  use control_flags, only: restart
  use lsda_mod, only: nspin
  use noncollin_module, only: noncolin
  implicit none
  integer, intent(in) :: ik
  !
  integer :: g0_pool, ib, ib1, ib2, sizeSubspace, shap(2), nRows, nBands, &
       num_local_plane_Waves,  i, j, ik_global, ngm_g_, degeneracy_group_counter, &
       ios, ik_irr, isym, nbnd_, ig_rotated, ig1, ig2, ig, ib1_, degeneracy_group(nbnd)
  integer, save :: nk_full=0, nk1_, nk2_, nk3_
  integer, allocatable :: gmap(:)
  integer, allocatable, save :: xk_equiv(:), xk_equiv_symmetry(:), xk_weight(:)
  real(dp) :: theta, rotation(3,3), inv_rotation(3,3), translation(3), diff, &
       arg, this_rotated_g(3), this_g(3), xk_crys(3), xk_irr_from_file(3), &
       umklapp_Vector(3), xk_irr_from_file_cart(3), eigenvalues(nbnd)
  real(dp), allocatable, save :: g_global(:,:), xk_full_cart(:,:), &
       xk_full_cryst(:,:)
  real(dp), allocatable :: et_irr(:), rotated_g_global(:,:), g_global_(:,:)
  complex(dp) :: correction, xc, unitary_matrix(nbnd,nbnd), delta_matrix(nbnd,nbnd)
  complex(dp), allocatable :: gauge_coefficients(:), evc_collected(:), &
       phases(:), evc_irreducible(:), evc_rotated(:), evc_test(:)
  character(len=64) :: file_name
  character(len=4) :: ichar
  logical, save :: first = .true.
  logical, allocatable :: is_band_degenerate(:), is_band_degenerate_(:)
  logical :: any_prob, in_scf, in_nscf
  integer, parameter :: i_unit = 52
  
  in_scf = (trim(calculation) == 'scf') ! .and. (restart) 
  in_nscf = (trim(calculation) == 'nscf') .or. (trim(calculation) == 'bands') &
       .or. (trim(calculation)=="none")
  
  if ( (.not. in_scf) .and. (.not. in_nscf) ) then
    return
  end if

  if ( nspin /= 1 ) then
    call errore("phoebe", "Spin is not yet supported in phoebe", 1)
  end if
    if ( noncolin ) then
    call errore("phoebe", "Spin-orbit is not yet supported in phoebe", 1)
  end if
  ! Both can be supported, but we need to consider magnetic symmetries

  !
  ! Things only work without offset
  if ( k1 /= 0 .or. k2 /= 0 .or. k3 /= 0 ) then
    call errore("phoebe", "No k-point offset allowed",1)
    ! that's because the k+q mesh would not be commensurate
  end if
  !
  ! figure out which process holds the G=0 vector
  g0_pool = 0
  if ( gstart == 2 ) then
    g0_pool = me_pool
  end if
  call mp_sum(g0_pool, intra_pool_comm)
  if ( me_pool == g0_pool ) then
    if ( sum(g(:,1)**2)>1.0e-8 ) then
      call errore("set_wavefunction_gauge", "Unexpectedly G/=0", 1)
    end if
  end if
  !
  ! Sanity check, in case of particular band parallelizations
  shap = shape(evc)
  nRows = shap(1)
  nBands = shap(2)
    if ( nBands /= nbnd ) then
    call errore("set_wavefunction_gauge", "Unexpected",1)
  end if
  !
  !--------------------------------------------------------------
  ! Step 1
  ! Fix the G=0 plane wave coefficient to be positive and real
  ! for degenerate eigenstates, only the first band of the degenerate group
  ! can be set to have positive G=0 coefficient
  
  allocate(gauge_coefficients(nbnd))
  allocate(is_band_degenerate(nbnd))
  is_band_degenerate = .false.
  gauge_coefficients = cmplx(0.,0.,kind=dp)
  degeneracy_group = 0
  degeneracy_group_counter = 0
  if ( me_pool == g0_pool ) then
    ib = 0
    ! loop on bands
    do while ( ib < nbnd )
      ib = ib + 1      
      ! check if the band is degenerate, and get the degenerate subspace size
      sizeSubspace = 1
      do ib2 = ib+1,nbnd
        if (abs(et(ib,ik) - et(ib2,ik)) > 0.0001 / ryToeV) then  ! 0.1 meV
          exit
        end if
        sizeSubspace = sizeSubspace + 1;
      end do
      degeneracy_group_counter = degeneracy_group_counter + 1
      if (sizeSubspace == 1) then
        gauge_coefficients(ib) = evc(1,ib)
        degeneracy_group(ib) = degeneracy_group_counter
      else
        gauge_coefficients(ib:ib+sizeSubspace-1) = evc(1,ib)
        is_band_degenerate(ib:ib+sizeSubspace-1) = .true.
        degeneracy_group(ib:ib+sizeSubspace-1) = degeneracy_group_counter
      end if
      ib = ib + sizeSubspace - 1;      
    end do ! band loop    
  end if
 
  ! test: let's fix the gauge Re{c(G=0)}>0, Im{c(G=0)}=0
  call mp_bcast(gauge_coefficients, g0_pool, intra_pool_comm)
  ! for every band, find max
  do ib = 1,nbnd
    xc = gauge_coefficients(ib)
    ! Now, I compute and impose the gauge
    ! z = |z| e^(i theta)
    theta = atan( dimag(xc) / real(xc) )
    if ( real(xc) < 0. ) then ! rotation to make c(G) positive
      theta = theta + pi
    end if
    correction = cmplx(cos(-theta), sin(-theta), kind=dp)
    ! Impose gauge
    evc(:,ib) = evc(:,ib) * correction
  end do
  deallocate(gauge_coefficients)

  !----------------------------------------------------------------
  ! STEP 2:
  ! The first time this subroutine is called, we can setup some quantities
  ! in particular, the reducible list of kpoints, and the list of G-vectors
  
  ! we save the global list of g-vectors
  if ( first ) then
    first = .false.

    ! setup list of G-vectors
    ! Note: g is a local sublist of G-vectors, so we need an allreduce op
    allocate(g_global(3,ngm_g))
    g_global = 0.d0
    do i = 1,ngk(ik)
      g_global(:,ig_l2g(i)) = g(:,i)
    end do
    call mp_sum(g_global, intra_pool_comm)

    if ( trim(calculation) == 'scf' ) then
      ! setup kpoint grid parameters
      nk1_ = nk1
      nk2_ = nk2
      nk3_ = nk3

      ! Save info on the G grid, to check that runs are consistent        
      if ( my_pool_id == 0 .and. me_pool == root_pool ) then
        file_name = trim(tmp_dir) // trim(prefix) // ".phoebe.scf.0000.dat"
        open(unit = i_unit, file = TRIM(file_name), form = 'unformatted', &
             access = 'sequential', status = 'replace', iostat = ios)
        write(i_unit) ngm_g
        write(i_unit) nk1, nk2, nk3
        write(i_unit) g_global
        close(i_unit) 
      end if
      
    else
      ! here we read the file generated by the scf run, and check consistency
      ! check that the code restarts with the same grid
      if ( my_pool_id == 0 .and. me_pool == root_pool ) then
        if ( calculation == "none" ) then ! case of phonons
          file_name = trim(tmp_dir) // "/../../" // trim(prefix) // ".phoebe.scf.0000.dat"
        else
          file_name = trim(tmp_dir) // trim(prefix) // ".phoebe.scf.0000.dat"
        end if
        open(unit=i_unit, file=trim(file_name), form = 'unformatted', &
             access='sequential', status='old', iostat=ios)
        if ( ios /= 0 ) call errore("phoebe", "file not found", 1)
        read(i_unit) ngm_g_ ! # of global G vectors
        read(i_unit) nk1_, nk2_, nk3_ ! kgrid mesh

        allocate(g_global_(3,ngm_g_)) ! list of G-vectors
        
        read(i_unit) g_global_
        ! I checked that g vectors are the same
        
        close(i_unit) 
      end if
      call mp_bcast(nk1_, root_pool, intra_pool_comm)
      call mp_bcast(nk2_, root_pool, intra_pool_comm)
      call mp_bcast(nk3_, root_pool, intra_pool_comm)
      
      if ( ngm_g_ /= ngm_g ) then
        print*, ngm_g_, ngm_g, "!"
        call errore("phoebe", "Different number of Gvectors in restart", 1)
      end if
      
    end if
    
    if ( nk1_ <= 0 .or. nk2_ <= 0 .or. nk3_ <= 0 ) then
      call errore("phoebe","k-point grid not found. Using kpoints automatic?",1)
    end if
   
    !
    ! full grid and symmetry analysis
    nk_full = nk1_*nk2_*nk3_ ! total # of reducible kpoints
    allocate(xk_full_cryst(3,nk_full))
    allocate(xk_equiv(nk_full))
    allocate(xk_weight(nk_full))
    allocate(xk_equiv_symmetry(nk_full)) ! index of symmetry s.t. S(idx)*k^irr = k
    call find_irreducible_grid(nk1_, nk2_, nk3_, k1, k2, k3, xk_full_cryst, &
         xk_equiv, xk_weight, xk_equiv_symmetry)
    
    allocate(xk_full_cart(3,nk_full)) ! same as xk_full_cryst, but in cartesian coords
    xk_full_cart = xk_full_cryst
    do i = 1,nk_full
      call cryst_to_cart(1, xk_full_cart(:,i), at, -1)
    end do

  end if

  ! define ik_global as the index of this point
  ! in the full list of points that we use internally

  ! xk(:,ik) is the irreducible kpoint being computed now, in cartesian coords
  ! here we also get it in crystal coords
  DO i = 1, 3 
    xk_crys(i) = at(1,i)*xk(1,ik) + at(2,i)*xk(2,ik) + at(3,i)*xk(3,ik)
  end do
  ! find index of the irred. point in my global list of points
  ! this also folds point correctly in 1st BZ
  call get_point_index(ik_global, xk_crys, nk1_,nk2_,nk3_, k1,k2,k3, .false.)

  !---------------------------------------------
  ! Step 3:
  ! we split the code in two cases
  ! 1) If we are doing an scf calculation, we write the wavefunction to file
  ! 2) In a nscf calculation, we read the wavefunction and fix the gauge
    
  ! if we run the scf, we need to save info to file
  if ( in_scf ) then ! -----------------------

    ! print*, "Gauge saving"
    
    if ( me_pool == root_pool ) then
      write(ichar,"(I4.4)") ik_global
      file_name = trim(tmp_dir)//trim(prefix)//".phoebe.scf."//ichar//".dat"
      open(unit = i_unit, file = TRIM(file_name), form = 'unformatted', &
           access = 'sequential', status = 'replace', iostat = ios)
      write(i_unit) xk_crys
      write(i_unit) nbnd
      write(i_unit) et(:,ik)
      write(i_unit) is_band_degenerate(:)
    end if

    allocate(evc_collected(ngm_g))
    do ib1 = 1,nbnd
        ! I want to reorder evc to be aligned with global list of g_vectors
        !
        evc_collected = cmplx(0.,0.,kind=dp)
        do ig = 1,ngk(ik)
          evc_collected(ig_l2g(igk_k(ig,ik))) = evc(ig,ib1)
          ! Note: igk_k maps the ordering of g-vectors at k in
          ! the ordering of g(:,:) (|G|^2 ordering vs |G+k|^2 ordering)
          ! ig_l2g maps the G-index from MPI-local to MPI-global
        end do
        call mp_sum(evc_collected, intra_pool_comm)
        !
        if ( me_pool == root_pool ) then
          write(i_unit) ib1
          write(i_unit) evc_collected
        end if
        !
    end do ! band loop
    deallocate(evc_collected)
    
    if ( me_pool == root_pool ) then
      close(i_unit)      
    end if
    
  else ! nscf ! ----------------------------------------------------

    ! print*, "Gauge fixing"
    
    ! Here is the funny part
    ! We are doing this nscf calculation on a full grid of points
    
    ! First, we find the index of the current kpoint
    ! in the list of irreducible points
    
    ik_irr = xk_equiv(ik_global)
    isym = xk_equiv_symmetry(ik_global)
    
    rotation = sr(:,:,isym) ! such that R*k^red = k^irr, in cartesian space
    inv_rotation = transpose(rotation) ! Rotations are unitary

    ! Read from file the energies of the irreducible point
    allocate(et_irr(nbnd))
    allocate(is_band_degenerate_(nbnd))
    et_irr = 0.0d0
    if ( me_pool == root_pool ) then
      ! read info on g vectors and symmetries
      write(ichar,"(I4.4)") ik_irr
      if ( calculation == "none" ) then ! case of phonons
        file_name = trim(tmp_dir) // "/../../" //trim(prefix)// ".phoebe.scf."//ichar//".dat"
      else
        file_name = trim(tmp_dir) // trim(prefix) // ".phoebe.scf."//ichar//".dat"
      end if
      
      open(unit=i_unit, file=TRIM(file_name), form='unformatted', &
           access='sequential', status='old', iostat=ios)
      read(i_unit) xk_irr_from_file
      read(i_unit) nbnd_
      read(i_unit) et_irr(:)
      read(i_unit) is_band_degenerate_(:)
    end if
    call mp_bcast(et_irr, me_pool, intra_pool_comm)
    call mp_bcast(is_band_degenerate_, me_pool, intra_pool_comm)
    call mp_bcast(nbnd_, me_pool, intra_pool_comm)

    ! find the Umklapp vector between the current k and the k from file
    umklapp_vector = matmul(s(:,:,isym),xk_crys(:)) - xk_irr_from_file(:)
    ! make sure it has integer values
    if ( sum(abs(umklapp_vector)) - nint(sum(abs(umklapp_vector))) > 1.0e-5 ) then
      call errore("phoebe", "Umklapp with non-integer values. Wrong kPoints?", 1)
      ! this also makes us test whether we are reading the correct wavevector
    end if
    ! we use cartesian coordinates
    umklapp_vector = matmul(bg,umklapp_vector)

    ! Sanity check: et_irr should be roughly the same of the current energies et
    ! Note that this check should change when adding spin
    if ( nbnd_ /= nbnd ) call errore("phoebe","scf and nscf run with different bands",1)
    do ib = 1,nbnd
      if ( abs(et_irr(ib) - et(ib,ik)) > 1.0e-4 ) then
        call errore("phoebe","incorrect symmetry on energies",1)
      end if
      if ( is_band_degenerate_(ib) .neqv. is_band_degenerate(ib) ) then
        call errore("phoebe","Degeneracy info not consistent with the scf",1)
      end if        
    end do
    
    ! we reinforce the symmetry on energies, to reduce numerical noise
    ! Note: this causes some noise
    ! et(:,ik) = et_irr(:)

    ! deallocate(is_band_degenerate, et_irr)

    ! Note: we still need to fix gauge for  all bands
    ! except for those kpoints that coincide with the irreducible ones
    ! if no degeneracy, no need to rotate wavefunction, and return
!    if ( isym == 1 ) then ! this is the index of the identity matrix
!      if ( me_pool == root_pool ) close(i_unit)
!      return
!    end if

    !-----------------------------------------------------
    ! Actual gauge fixing of degenerate states starts here

    ! build the list of rotated G vectors
    allocate(rotated_g_global(3,ngm_g))
    do ig1 = 1,ngm_g
      ! note that G vectors are in cartesian coordinates, in units of 2Pi/alat
      this_g(:) = g_global(:,ig1)
      this_rotated_g(:) = matmul(inv_rotation, this_g(:)) + umklapp_vector
      rotated_g_global(:,ig1) = this_rotated_g(:)
    end do

    !----------------------------------------
    ! find the index mapping between G -> R G
    allocate(gmap(ngm_g))
    gmap = 0
    do ig1 = 1,ngm_g
      ! this search is expensive, so we go parallel within pool
      if ( mod(ig1-1,nproc_pool) /= me_pool ) cycle
      this_g(:) = g_global(:,ig1)
      ! it seems that the first g-vector is 0.
      ! then there are non-zero vectors, then it's again a lot of zero vectors
      ! here we make sure that the first gmap refers to 0 g-vector
      if ( (sum(this_g**2) < 1.0e-6) .and. (ig1>10) ) cycle
      ig_rotated = 0
      do ig2 = 1,ngm_g
        diff = sum(( g_global(:,ig1) - rotated_g_global(:,ig2) )**2)
        if ( diff < 1.0e-6 ) then
          gmap(ig2) = ig1
          exit
        end if
      end do
    end do
    call mp_sum(gmap, intra_pool_comm)
    deallocate(rotated_g_global)

    !-------------------------------------------------------
    ! compute phases due to translation
    ! these are the same for all bands
    allocate(phases(ngm_g))
    phases = cmplx(0.,0.,kind=dp)
    ! ft contains  fractional translations in crystal axis
    ! transform in cartesian coordinates
    translation = matmul(at,ft(:,isym))

    ! These are phases that the pw coefficients gain on rotation
    xk_irr_from_file_cart = matmul(bg,xk_irr_from_file)
    do ig1 = 1,ngm_g
      if ( mod(ig1-1,nproc_pool) /= me_pool ) cycle
      if ( (sum(g_global(:,ig1)**2) < 1.0e-6) .and. (ig1>10) ) cycle
      ! arg = - tpi * dot_product( matmul(rotation,xk_irr_from_file_cart) + g_global(:,ig1) , translation )
      ! arg = - tpi * dot_product( xk(:,ik) + g_global(:,ig1) , translation )
      arg = tpi * dot_product( matmul(rotation,xk_irr_from_file_cart) + g_global(:,ig1) , translation )
      phases(ig1) = cmplx(cos(arg), sin(arg), kind=dp)
    end do
    call mp_sum(phases, intra_pool_comm)
    ! print*, translation, isym
    ! write(*,"(3F12.5)") matmul(s(:,:,isym),xk_crys(:)) - xk_irr_from_file(:)
    ! write(*,"(6F12.5)") xk_crys(:), xk_irr_from_file(:)

    ! time_invariance = .false.
    
    !-------------------------------------------------------
    ! Rotate wavefunction plane wave coefficients
    
    allocate(evc_irreducible(ngm_g))
    allocate(evc_rotated(ngm_g))

    ! create a matrix that rotates the wavefunctions according to symmetry
    unitary_matrix = cmplx(0.,0.,kind=dp)
   
    do ib1 = 1,nbnd
      evc_irreducible(:) = cmplx(0.,0.,kind=dp)
      evc_rotated(:) = cmplx(0.,0.,kind=dp)
        
      ! Read the wavefunction at the irreducible point
      ! remember that the array is not distributed, and ordered like g_global

      if ( me_pool == root_pool ) then
        ! read info on g vectors and symmetries
        read(i_unit) ib1_
        read(i_unit) evc_irreducible          
      end if
      call mp_bcast(ib1_, me_pool, intra_pool_comm)
      if ( ib1 /= ib1_ ) then
        call errore("phoebe", "Unexpected degenerate band", 1)
      end if
      call mp_bcast(evc_irreducible, me_pool, intra_pool_comm)

      !---------------------------------------------------------------------
      ! Rotate the wavefunction

      ! Giannozzi's quote:
      ! no, no!  evc(i,n) = i-th component of the n-th band; the i-th component
      ! corrsponds to (k+G)(i) = k(ik)+G(igk_k(i,ik))   where ik is the index of
      ! k-points. G is the array of G-vectors

      ! create evc_rotated: a wfc, not parallel distributed, with G-vectors
      ! ordered like g_global, which satisfies the symmetries
      do ig = 1,ngm_g
        if ( gmap(ig) > 0 ) then
          evc_rotated(ig) = evc_irreducible(gmap(ig)) * phases(ig)
        end if
      end do

      allocate(evc_collected(ngm_g))
      evc_collected = cmplx(0.,0.,kind=dp)
      do ig = 1,ngk(ik)
        evc_collected(ig_l2g(igk_k(ig,ik))) = evc(ig,ib1)
      end do
      ! print*, ib1, degeneracy_group(ib1), et(ib1,ik), abs(sum(conjg(evc_rotated)*evc_collected))
      deallocate(evc_collected)
      
      ! substitute back in QE array
      ! In doing so, we must reorder the G-vectors and distribute the wavefunction
      !do ig = 1,ngk(ik)
      !  evc(ig,ib1) = evc_rotated(ig_l2g(igk_k(ig,ik)))
      !end do
      
     do ib2 = 1,nbnd
       if ( degeneracy_group(ib1) == degeneracy_group(ib2) ) then
         do ig = 1,ngk(ik)
           unitary_matrix(ib1,ib2) = unitary_matrix(ib1,ib2) + &
                conjg(evc(ig,ib2)) * evc_rotated(ig_l2g(igk_k(ig,ik)))
         end do
       end if
     end do

      ! print*, abs(sum(conjg(evc(:,ib1)) * evc(:,ib1)))
      ! print*, abs(sum(conjg(evc(:,ib1)) * evc_rotated(ig_l2g(igk_k(:,ik)))))
      ! print*, abs(sum(conjg(evc(:,ib1)) * evc_irreducible(:)))
      ! print*, abs(sum(conjg(evc_irreducible(:)) * evc_rotated(ig_l2g(igk_k(:,ik)))))
      ! do ig = 1,30
      !   print*,ig,ig_l2g(ig),igk_k(ig,ik),gmap(ig),evc(ig,ib1),evc_rotated(ig_l2g(igk_k(ig,ik)))
      ! end do
      
    end do
    
    call mp_sum(unitary_matrix, intra_pool_comm)

    ! This is actually wrong!
    ! unitary_matrix = ( unitary_matrix + conjg(transpose(unitary_matrix)) ) * 0.5d0

!    do ib1 = 1,nbnd
!      write(*,"(24F11.7)") unitary_matrix(ib1,:)
!    end do
        
!    do ib1 = 1,nbnd
!      print*, ib1, evc(1,ib1), sum(unitary_matrix(ib1,:)*evc(1,:))
!    end do
!    print*, "!!"

    !------------------------------
    ! reinforce matrix is unitary
    delta_matrix = cmplx(0.,0.,kind=dp)
    do ib1 = 1,nbnd
      delta_matrix(ib1,ib1) = cmplx(1.,0.,kind=dp)
    end do
    delta_matrix = delta_matrix - matmul(unitary_matrix,conjg(transpose(unitary_matrix)))
    delta_matrix = matmul(matmul(conjg(transpose(unitary_matrix)),delta_matrix),unitary_matrix)
    do ib1 = 1,nbnd
      delta_matrix(ib1,ib1) = delta_matrix(ib1,ib1) + cmplx(1.,0.,kind=dp)
    end do
    call zpotrf("L",nbnd,delta_matrix,nbnd,ib2)
    if ( ib2 /= 0 ) call errore("phoebe","Cholesky failed",1)
    unitary_matrix = matmul(unitary_matrix,delta_matrix)
    
!    call mp_bcast(unitary_matrix, root_pool, intra_pool_comm)

!    do ib1 = 1,nbnd
!      print*, ib1, evc(1,ib1), sum(unitary_matrix(ib1,:)*evc(1,:))
!    end do

!    do ib1 = 1,nbnd
!      print*, ib1, sum(conjg(evc(:,ib1))*evc(:,ib1))
!    end do
    
!    do ig = 1,ngk(ik)
!      evc(ig,:) = matmul(unitary_matrix,evc(ig,:))
!    end do    

    unitary_matrix = matmul( unitary_matrix, conjg(transpose(unitary_matrix)) ) 
!    do ib1 = 1,nbnd
!      write(*,"(24F11.7)") unitary_matrix(ib1,:)
!    end do
    
    deallocate(is_band_degenerate, et_irr)
    deallocate(evc_rotated, evc_irreducible, gmap, phases)
    if ( me_pool == root_pool ) close(i_unit)

  end if
  
  return
end subroutine set_wavefunction_gauge





!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE c_bands( iter )
  !----------------------------------------------------------------------------
  !! Driver routine for the Hamiltonian diagonalization ones.
  !! It reads the Hamiltonian and an initial guess of the wavefunctions
  !! from a file and computes initialization quantities for the
  !! diagonalization routines.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, isolve, restart
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions,        ONLY : evc
  USE bp,                   ONLY : lelfield
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE check_stop,           ONLY : check_stop_now
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  !! iteration index
  !
  ! ... local variables
  !
  REAL(DP) :: avg_iter
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios
  ! ik : counter on k points
  ! ik_: k-point already done in a previous run
  LOGICAL :: exst
  !
  !
  CALL start_clock( 'c_bands' ); !write (*,*) 'start c_bands' ; FLUSH(6)
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands( ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  ! ... (not needed for a single k-point: this is done in wfcinit, 
  ! ...  directly from file, in order to avoid wasting memory)
  !
  DO ik = 1, ik_
     IF ( nks > 1 .OR. lelfield ) &
        CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
  ENDDO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSEIF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")')
  ELSEIF ( isolve == 2 ) THEN
     WRITE( stdout, '(5X,"PPCG style diagonalization")')
  ELSEIF ( isolve == 3 ) THEN
     WRITE( stdout, '(5X,"ParO style diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
  ENDIF
  !
  ! ... For each k point diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) CALL phase_factor(ik)
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     CALL g2_kin( ik )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF ( nks > 1 .OR. lelfield ) &
          CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... diagonalization of bands for k-point ik
     !
     call diag_bands ( iter, ik, avg_iter )
     !
     ! ... save wave-functions to be used as input for the
     ! ... iterative diagonalization of the next scf iteration
     ! ... and for rho calculation
     !
     IF ( nks > 1 .OR. lelfield ) &
          CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     !
     IF (ik <= nkdum) THEN
        IF (check_stop_now()) THEN
           CALL save_in_cbands( ik, ethr, avg_iter, et )
           RETURN
        ENDIF
     ENDIF
     !
  ENDDO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, &
       '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
       ethr, avg_iter
  !
  CALL stop_clock( 'c_bands' ); !write (*,*) 'stop c_bands' ; FLUSH(6)
  !
  RETURN
  !
END SUBROUTINE c_bands
!
!----------------------------------------------------------------------------
SUBROUTINE diag_bands( iter, ik, avg_iter )
  !----------------------------------------------------------------------------
  !! Driver routine for diagonalization at each k-point. Types of iterative
  !! diagonalizations currently in use:
  !
  !! * Davidson algorithm (all-band);
  !! * Conjugate Gradient (band-by-band);
  !! * Projected Preconditioned Conjugate Gradient (block);
  !! * Parallel Orbital update (all-band).
  !
  !! Internal procedures:
  !
  !! * \(\textrm{diag_bands_gamma}\)(): optimized algorithms for gamma sampling
  !!                                    of the BZ (real Hamiltonian);
  !! * \(\textrm{diag_bands_k}\)(): general algorithm for arbitrary BZ sampling
  !!                                (complex Hamiltonian);
  !! * \(\textrm{test_exit_cond}\)(): the test on the iterative diagonalization.
  !
  USE kinds,                ONLY : DP
  USE buffers,              ONLY : get_buffer
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfc, iunefieldp, iunefieldm
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE gvect,                ONLY : gstart
  USE wvfct,                ONLY : g2kin, nbndx, et, nbnd, npwx, btype
  USE control_flags,        ONLY : ethr, lscf, max_cg_iter, max_ppcg_iter, isolve, &
                                   gamma_only, use_para_diag
  USE noncollin_module,     ONLY : npol
  USE wavefunctions,        ONLY : evc
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE scf,                  ONLY : v_of_0
  USE bp,                   ONLY : lelfield, evcel, evcelp, evcelm, bec_evcel, &
                                   gdir, l3dstring, efield, efield_cry
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE klist,                ONLY : nks, ngk
  USE mp_bands,             ONLY : nproc_bgrp, intra_bgrp_comm, inter_bgrp_comm, &
                                   my_bgrp_id, nbgrp
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  !! iteration index
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  REAL(KIND=DP), INTENT(INOUT) :: avg_iter
  !! average number of H*psi products
  !
  ! ... local variables
  !
  REAL(KIND=DP) :: cg_iter, ppcg_iter
  ! (weighted) number of iterations in Conjugate-Gradient
  INTEGER :: npw, ig, dav_iter, ntry, notconv, nhpsi
  ! number of iterations in Davidson
  ! number or repeated call to diagonalization in case of non convergence
  ! number of notconverged elements
  INTEGER :: ierr, ipw, ibnd, ibnd_start, ibnd_end
  !
  LOGICAL :: lrot
  ! .TRUE. if the wfc have already be rotated
  !
  INTEGER, PARAMETER :: sbsize = 5, rrstep = 7
  ! block dimensions used in PPCG 
  !
  ! Davidson diagonalization uses these external routines on groups of nvec bands
  EXTERNAL h_psi, s_psi, g_psi
  ! subroutine h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
  ! subroutine s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
  ! subroutine g_psi(npwx,npw,nvec,psi,eig)   computes G*psi -> psi
  !------------------------------------------------------------------------
  ! CG diagonalization uses these external routines on a single band
  EXTERNAL hs_1psi, s_1psi, hs_psi
  ! subroutine hs_1psi(npwx,npw,psi,hpsi,spsi)  computes H*psi and S*psi
  ! subroutine s_1psi(npwx,npw,psi,spsi)        computes S*psi (if needed)
  ! In addition to the above the initial wfc rotation uses h_psi, and s_psi
  !------------------------------------------------------------------------
  ! PPCG diagonalization uses these external routines on groups of bands
  ! subroutine h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
  ! subroutine s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
  !------------------------------------------------------------------------
  ! ParO diagonalization uses these external routines on a single band
  ! subroutine hs_1psi(npwx,npw,psi,hpsi,spsi)  computes H*psi and S*psi
  ! subroutine g_1psi(npwx,npw,psi,eig)         computes G*psi -> psi
  ! In addition to the above the initial wfc rotation uses h_psi, and s_psi
  external g_1psi

  ALLOCATE( h_diag( npwx, npol ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' diag_bands ', ' cannot allocate h_diag ', ABS(ierr) )
  !
  ALLOCATE( s_diag( npwx, npol ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' diag_bands ', ' cannot allocate s_diag ', ABS(ierr) )
  !
  ipw=npwx
  CALL mp_sum(ipw, intra_bgrp_comm)
  IF ( nbndx > ipw ) &
     CALL errore ( 'diag_bands', 'too many bands, or too few plane waves',1)
  !
  ! ... allocate space for <beta_i|psi_j> - used in h_psi and s_psi
  !
  CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm )
  !
  npw = ngk(ik)
  IF ( gamma_only ) THEN
     !
     CALL diag_bands_gamma()
     !
  ELSE
     !
     CALL diag_bands_k()
     !
  ENDIF

  ! PHOEBE:
  ! Apply gauge of wavefunction
  call set_wavefunction_gauge(ik)
   
  !
  ! ... deallocate work space
  !
  CALL deallocate_bec_type( becp )
  DEALLOCATE( s_diag )
  DEALLOCATE( h_diag )
  !
  IF ( notconv > MAX( 5, nbnd / 4 ) ) THEN
     !
     CALL errore( 'c_bands', 'too many bands are not converged', 1 )
     !
  ELSEIF ( notconv > 0 ) THEN
     !
     WRITE( stdout, '(5X,"c_bands: ",I2, " eigenvalues not converged")' ) notconv
     !
  ENDIF
  !
  RETURN
  !
 CONTAINS
  !
  ! ... internal procedures
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_bands_gamma()
    !-----------------------------------------------------------------------
    !
    ! ... Diagonalization of a real Hamiltonian
    !
    IMPLICIT NONE
    !
    IF ( isolve == 1 .OR. isolve == 2 .OR. isolve == 3) THEN
       !
       ! ... (Projected Preconditioned) Conjugate-Gradient diagonalization
       !
       ! ... h_diag is the precondition matrix
       !
       IF ( isolve == 1 .OR. isolve == 2 ) THEN
          FORALL( ig = 1 : npw )
             h_diag(ig,1) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
          END FORALL
       ELSE
          FORALL( ig = 1 : npw )
             h_diag(ig, 1) = g2kin(ig) + v_of_0
          END FORALL
          CALL usnldiag( npw, h_diag, s_diag )
       END IF
       !
       ntry = 0
       !
       CG_loop : DO
          !
          IF ( isolve == 1 .OR. isolve == 2 ) THEN
             lrot = ( iter == 1 .AND. ntry == 0 )
             !
             IF ( .NOT. lrot ) THEN
                !
                CALL rotate_wfc( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
                !
                avg_iter = avg_iter + 1.D0
                !
             ENDIF
          ENDIF
          !
          IF ( isolve == 1 ) THEN
             CALL rcgdiagg( hs_1psi, s_1psi, h_diag, &
                         npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
                         ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
             !
             avg_iter = avg_iter + cg_iter
             !
          ELSE IF ( isolve == 2 ) THEN
             CALL ppcg_gamma( h_psi, s_psi, okvan, h_diag, &
                         npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
                         0.1d0*ethr, max_ppcg_iter, notconv, ppcg_iter, sbsize , rrstep, iter )
             !
             avg_iter = avg_iter + ppcg_iter
             !
          ELSE
             !
             CALL paro_gamma_new( h_psi, s_psi, hs_psi, g_1psi, okvan, &
                        npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
             !
             avg_iter = avg_iter + nhpsi/float(nbnd) 
             ! write (6,*) ntry, avg_iter, nhpsi
             !
          ENDIF
          !
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  CG_loop
          !
       ENDDO CG_loop
       !
    ELSE
       !
       ! ... Davidson diagonalization
       !
       ! ... h_diag are the diagonal matrix elements of the
       ! ... hamiltonian used in g_psi to evaluate the correction
       ! ... to the trial eigenvectors
       !
       h_diag(1:npw, 1) = g2kin(1:npw) + v_of_0
       !
       CALL usnldiag( npw, h_diag, s_diag )
       !
       ntry = 0
       !
       david_loop: DO
          !
          lrot = ( iter == 1 )
          !
          IF ( use_para_diag ) then
             !
!             ! make sure that all processors have the same wfc
             CALL pregterg( h_psi, s_psi, okvan, g_psi, &
                         npw, npwx, nbnd, nbndx, evc, ethr, &
                         et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                                   !    BEWARE gstart has been removed from call 
             !
          ELSE
             !
             CALL regterg (  h_psi, s_psi, okvan, g_psi, &
                         npw, npwx, nbnd, nbndx, evc, ethr, &
                         et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
                                   !    BEWARE gstart has been removed from call
          ENDIF
          !
          avg_iter = avg_iter + dav_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  david_loop
          !
       ENDDO david_loop
       !
    ENDIF
    !
    !
    RETURN
    !
  END SUBROUTINE diag_bands_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_bands_k()
    !-----------------------------------------------------------------------
    !! Complex Hamiltonian diagonalization.
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ipol
    REAL(DP) :: eps=0.000001d0
    ! --- Define a small number ---
    !
    !write (*,*) ' enter diag_bands_k'; FLUSH(6)
    IF ( lelfield ) THEN
       !
       ! ... save wave functions from previous iteration for electric field
       !
       evcel = evc
       !
       !... read projectors from disk
       !
       IF (.NOT.l3dstring .AND. ABS(efield)>eps ) THEN
          CALL get_buffer (evcelm(:,:,gdir), nwordwfc, iunefieldm, ik+(gdir-1)*nks)
          CALL get_buffer (evcelp(:,:,gdir), nwordwfc, iunefieldp, ik+(gdir-1)*nks)
       ELSE
          DO ipol = 1, 3
             IF ( ABS(efield_cry(ipol))>eps ) THEN
                CALL get_buffer( evcelm(:,:,ipol), nwordwfc, iunefieldm, ik+(ipol-1)*nks )
                CALL get_buffer( evcelp(:,:,ipol), nwordwfc, iunefieldp, ik+(ipol-1)*nks )
             ENDIF
          ENDDO
       ENDIF
       !
       IF ( okvan ) THEN
          !
          CALL allocate_bec_type( nkb, nbnd, bec_evcel )
          !
          CALL calbec( npw, vkb, evcel, bec_evcel )
          !
       ENDIF
       !
    ENDIF
    !
    !write (*,*) ' current isolve value ( 0 Davidson, 1 CG, 2 PPCG)', isolve; FLUSH(6)
    IF ( isolve == 1 .OR. isolve == 2 .OR. isolve == 3 ) THEN
       !
       ! ... (Projected Preconditioned) Conjugate-Gradient diagonalization
       !
       ! ... h_diag is the precondition matrix
       !
       !write (*,*) ' inside CG solver branch '
       !
       h_diag = 1.D0
       IF ( isolve == 1 .OR. isolve == 2) THEN
          FORALL( ig = 1 : npwx )
             h_diag(ig,:) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
          END FORALL
       ELSE
          FORALL( ig = 1 : npwx )
             h_diag(ig, :) = g2kin(ig) + v_of_0
          END FORALL
          CALL usnldiag( npw, h_diag, s_diag )
       ENDIF
       !
       ntry = 0
       !
       CG_loop : DO
          !
          IF ( isolve == 1 .OR. isolve == 2 ) THEN
             lrot = ( iter == 1 .AND. ntry == 0 )
             !
             IF ( .NOT. lrot ) THEN
                !
                CALL rotate_wfc( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
                !
                avg_iter = avg_iter + 1.D0
                !
             ENDIF
          ENDIF
          !
          IF ( isolve == 1) then
             CALL ccgdiagg( hs_1psi, s_1psi, h_diag, &
                         npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), &
                         ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
             !
             avg_iter = avg_iter + cg_iter
             !
          ELSE IF ( isolve == 2) then
             ! BEWARE npol should be added to the arguments
             CALL ppcg_k( h_psi, s_psi, okvan, h_diag, &
                         npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), &
                         0.1d0*ethr, max_ppcg_iter, notconv, ppcg_iter, sbsize , rrstep, iter )
             !
             avg_iter = avg_iter + ppcg_iter
             !
          ELSE 
             !
             CALL paro_k_new( h_psi, s_psi, hs_psi, g_1psi, okvan, &
                        npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), ethr, notconv, nhpsi )
             !
             avg_iter = avg_iter + nhpsi/float(nbnd) 
             ! write (6,*) ntry, avg_iter, nhpsi
             !
          ENDIF
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  CG_loop
          !
       ENDDO CG_loop
       !
    ELSE
       !
       ! ... Davidson diagonalization
       !
       ! ... h_diag are the diagonal matrix elements of the
       ! ... hamiltonian used in g_psi to evaluate the correction
       ! ... to the trial eigenvectors
       !
       DO ipol = 1, npol
          !
          h_diag(1:npw, ipol) = g2kin(1:npw) + v_of_0
          !
       ENDDO
       !
       CALL usnldiag( npw, h_diag, s_diag )
       !
       ntry = 0
       !
       david_loop: DO
          !
          lrot = ( iter == 1 )
          !
          IF ( use_para_diag ) THEN
             !
             CALL pcegterg( h_psi, s_psi, okvan, g_psi, &
                            npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                            et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
             !
          ELSE
             !
             CALL cegterg( h_psi, s_psi, okvan, g_psi, &
                           npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                           et(1,ik), btype(1,ik), notconv, lrot, dav_iter, nhpsi )
          ENDIF
          !
          avg_iter = avg_iter + dav_iter
          !
          ! ... save wave-functions to be used as input for the
          ! ... iterative diagonalization of the next scf iteration
          ! ... and for rho calculation
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT david_loop
          !
       ENDDO david_loop
       !
    ENDIF
    !
    IF ( lelfield .AND. okvan ) CALL deallocate_bec_type( bec_evcel )
    !
    RETURN
    !
  END SUBROUTINE diag_bands_k
  !
  !-----------------------------------------------------------------------
  FUNCTION test_exit_cond()
    !-----------------------------------------------------------------------
    !! This logical function is .TRUE. when iterative diagonalization
    !! is converged.
    !
    IMPLICIT NONE
    !
    LOGICAL :: test_exit_cond
    !
    !
    test_exit_cond = .NOT. ( ( ntry <= 5 ) .AND. &
         ( ( .NOT. lscf .AND. ( notconv > 0 ) ) .OR. &
         (       lscf .AND. ( notconv > 5 ) ) ) )
    !
  END FUNCTION test_exit_cond
  !
END SUBROUTINE diag_bands
!
!----------------------------------------------------------------------------
SUBROUTINE c_bands_efield( iter )
  !----------------------------------------------------------------------------
  !! Driver routine for Hamiltonian diagonalization under an electric field.
  !
  USE noncollin_module,     ONLY : npol
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : nberrycyc, fact_hepsi, &
                                   evcel, evcelp, evcelm, gdir, l3dstring,&
                                   efield, efield_cry
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : nbnd, npwx
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter
  !! iteration index
  !
  ! ... local variables
  !
  INTEGER :: inberry, ipol, ierr
  !
  !
  ALLOCATE( evcel ( npol*npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcel ', ABS( ierr ) )
  ALLOCATE( evcelm( npol*npwx, nbnd, 3  ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcelm ', ABS( ierr ) )
  ALLOCATE( evcelp( npol*npwx, nbnd, 3 ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcelp ', ABS( ierr ) )
  ALLOCATE( fact_hepsi(nks, 3), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate fact_hepsi ', ABS( ierr ) )
  !
  DO inberry = 1, nberrycyc
     !
     !...set up electric field hermitean operator
     !
     FLUSH(stdout)
     IF (.NOT.l3dstring) THEN
        CALL h_epsi_her_set (gdir, efield)
     ELSE
        DO ipol=1,3
           CALL h_epsi_her_set(ipol, efield_cry(ipol))
        ENDDO
     ENDIF
     FLUSH(stdout)
     !
     CALL c_bands( iter )
     !
  ENDDO
  !
  DEALLOCATE( fact_hepsi )
  DEALLOCATE( evcelp )
  DEALLOCATE( evcelm )
  DEALLOCATE( evcel  )
  !
  RETURN
  !
END SUBROUTINE c_bands_efield
!
!------------------------------------------------------------------------------
SUBROUTINE c_bands_nscf( )
  !----------------------------------------------------------------------------
  !! Driver routine for Hamiltonian diagonalization routines
  !! specialized to non-self-consistent calculations (no electric field).
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE basis,                ONLY : starting_wfc
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, restart, isolve, io_level, iverbosity
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions,        ONLY : evc
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE check_stop,           ONLY : check_stop_now
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(DP) :: avg_iter, ethr_
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios
  ! ik_: k-point already done in a previous run
  ! ik : counter on k points
  LOGICAL :: exst
  !
  REAL(DP), EXTERNAL :: get_clock
  !
  !
  CALL start_clock( 'c_bands' )
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands( ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  !
  DO ik = 1, ik_
     CALL get_buffer( evc, nwordwfc, iunwfc, ik )
  ENDDO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSEIF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")' )
  ELSEIF ( isolve == 2 ) THEN
     WRITE( stdout, '(5X,"PPCG style diagonalization")' )
  ELSEIF ( isolve == 3 ) THEN
     WRITE( stdout, '(5X,"ParO style diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve )
  ENDIF
  !
  ! ... For each k point (except those already calculated if restarting)
  ! ... diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) CALL phase_factor(ik)
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     CALL g2_kin( ik )
     ! 
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
          CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... calculate starting  wavefunctions
     !
     IF ( iverbosity > 0 .AND. npool == 1 ) THEN
        WRITE( stdout, 9001 ) ik, nks
     ELSEIF ( iverbosity > 0 .AND. npool > 1 ) THEN
        WRITE( stdout, 9002 ) ik, nks
     ENDIF
     !
     IF ( TRIM(starting_wfc) == 'file' ) THEN
        !
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        !
     ELSE
        !
        CALL init_wfc( ik )
        !
     ENDIF
     !
     ! ... diagonalization of bands for k-point ik
     !
     CALL diag_bands( 1, ik, avg_iter )
     !
     ! ... save wave-functions (unless disabled in input)
     !
     IF ( io_level > -1 ) CALL save_buffer( evc, nwordwfc, iunwfc, ik )
     !
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     IF (ik <= nkdum) THEN
        !
        ! ... stop requested by user: save restart information,
        ! ... save wavefunctions to file
        !
        IF ( check_stop_now() ) THEN
           CALL save_in_cbands( ik, ethr, avg_iter, et )
           RETURN
        ENDIF
        !
     ENDIF
     !
     ! report about timing
     !
     IF ( iverbosity > 0 ) THEN
        WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
        FLUSH( stdout )
     ENDIF
     !
  ENDDO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, '(/,5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)' ) &
       ethr, avg_iter
  !
  CALL stop_clock( 'c_bands' )
  !
  RETURN
  !
  ! formats
  !
9002 FORMAT(/'     Computing kpt #: ',I5, '  of ',I5,' on this pool' )
9001 FORMAT(/'     Computing kpt #: ',I5, '  of ',I5 )
9000 FORMAT( '     total cpu time spent up to now is ',F10.1,' secs' )
  !
END SUBROUTINE c_bands_nscf
