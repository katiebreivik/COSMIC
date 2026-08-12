    subroutine assign_commons()
        use track_support
        implicit none
        
        !to assign common variables when METISSE is used with COSMIC
          
        REAL(dp) :: ecsn,ecsn_mlow
        COMMON /SNVARS1/ ecsn,ecsn_mlow
         
        real(dp) :: d
        
        if(front_end == COSMIC) then
        ! use inputs from COSMIC
        
            if (Mec_core > 0.d0) ecsn = Mec_core
            d = (Mec_core-Mup_core)
            if (Mup_core > 0.d0 .and. d>tiny ) ecsn_mlow = Mup_core
            
        else
            print*,'Error: Front end mismtach in assign commons'
            print*,'expected 2 (COSMIC); got ', front_end
        endif

    end subroutine

    subroutine get_bhspin(bhspin,id)
        use track_support, only: tarr,dp
        implicit none
        integer, intent(in) :: id
        real(dp), intent(out) :: bhspin

        bhspin = tarr(id)% pars% bhspin
    end subroutine
    
    subroutine check_error(err)
        use track_support, only: code_error
        integer, intent(out) :: err
        err = 0
        if(code_error) err = 1
    end subroutine
    
    
    subroutine assign_error()
        use track_support, only: code_error
        code_error = .true.
    end subroutine
      
    subroutine initialize_metisse_front_cmc()
    ! passing strings with c/cmc is not very realiable
    ! we set front end like this avoid possible seg faults
        call initialize_front_end('cosmic')
    end subroutine

    subroutine get_COSMIC_input()
        use track_support
        use z_support, only: Z_accuracy_limit, get_csafe_string
        
        ! takes inputs from cosmic and assigns them
        ! to appropiate variables in METISSE
        
        character(len=strlen) :: path_to_tracks, path_to_he_tracks
        real(dp) :: z_match_limit
        LOGICAL METISSE_verbose
        COMMON/ METISSEVARS/ path_to_tracks,path_to_he_tracks,&
                     z_match_limit, METISSE_verbose
        
        ! remove the null charcater if any
        call get_csafe_string(path_to_tracks,METALLICITY_DIR)
        call get_csafe_string(path_to_he_tracks,METALLICITY_DIR_HE)
        Z_accuracy_limit = z_match_limit
        verbose = METISSE_verbose
    
    end subroutine

    subroutine load_metisse_tracks_cmc(z, zpars, ierr)
    ! CMC has no Python in its process, so it can't populate METISSE's
    ! track arrays the way COSMIC's own evolve.py does (read_MIST_track +
    ! c_m_interface.set_tracks_from_python, before zcnsts is ever called).
    ! front_end==COSMIC's zcnsts path assumes exactly that and segfaults
    ! otherwise. Every other front_end reads tracks from disk directly in
    ! Fortran (already-working code, see METISSE_zcnsts.f90's "case default"),
    ! so: borrow that path just long enough to load the tracks, using the
    ! same METALLICITY_DIR/METALLICITY_DIR_HE/Z_accuracy_limit CMC already
    ! passes in via get_COSMIC_input(), then switch back to COSMIC so
    ! METISSE_hrdiag uses COSMIC's own remnant/kick physics as before.
        use track_support, only: dp, strlen, amuse_metallicity_dir, amuse_metallicity_dir_he, verbose
        use z_support, only: Z_accuracy_limit, get_csafe_string
        implicit none
        real(dp), intent(in) :: z
        real(dp), intent(out) :: zpars(20)
        integer, intent(out) :: ierr

        character(len=strlen) :: path_to_tracks, path_to_he_tracks
        real(dp) :: z_match_limit
        logical :: METISSE_verbose
        COMMON/ METISSEVARS/ path_to_tracks,path_to_he_tracks, z_match_limit, METISSE_verbose

        ! NOTE: get_COSMIC_input() (below in this file) is NOT usable here --
        ! it writes straight into METALLICITY_DIR/Z_accuracy_limit/verbose,
        ! but those get reset by metisse_defaults.inc, included inside
        ! zcnsts's own first-time setup, which runs AFTER this call returns.
        ! amuse_metallicity_dir[_he] are read by METISSE_zcnsts's case(AMUSE)
        ! *after* that reset, which is exactly the ordering we need.
        call get_csafe_string(path_to_tracks, amuse_metallicity_dir)
        call get_csafe_string(path_to_he_tracks, amuse_metallicity_dir_he)

        call initialize_front_end('amuse')
        call zcnsts(z, zpars)
        call initialize_front_end('cosmic')

        ! Z_accuracy_limit/verbose get the same metisse_defaults.inc reset,
        ! but case(AMUSE) has no equivalent restore for them (only the two
        ! path variables above) -- so the very first metallicity-file match
        ! just above used METISSE's own default tolerance/verbosity, not
        ! CMC's configured values. Re-apply CMC's values now so every
        ! subsequent zcnsts call (reload-decision checks) uses them correctly.
        Z_accuracy_limit = z_match_limit
        verbose = METISSE_verbose

        call assign_commons()
        call check_error(ierr)
    end subroutine load_metisse_tracks_cmc

    logical function check_path_change() result (load_tracks)
        use track_support, only: strlen,METALLICITY_DIR, METALLICITY_DIR_HE
        use z_support, only: get_csafe_string

        character(len=strlen) :: path_to_tracks, path_to_he_tracks
        COMMON/ METISSEVARS/ path_to_tracks,path_to_he_tracks

        character(len=strlen) :: string1,string2
        load_tracks = .false.

        ! remove the null charcater if any
        call get_csafe_string(path_to_tracks,string1)
        call get_csafe_string(path_to_he_tracks, string2)

        if((trim(path_to_tracks)/=trim(METALLICITY_DIR)) .or. &
            (trim(path_to_he_tracks)/=trim(METALLICITY_DIR_HE))) load_tracks = .true.
    end function

    
