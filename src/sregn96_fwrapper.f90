! Fortran bind(C) wrappers for the RayleighWaveKernel subroutines.
! These allow the C++ test to call the original Fortran implementation
! for comparison.

subroutine sregn96_c(thk, vp, vs, rhom, nlayer, t, cp, cg, &
                     dispu, dispw, stressu, stressw, &
                     dc2da, dc2db, dc2dh, dc2dr, iflsph) &
    bind(C, name="sregn96_")
  use iso_c_binding
  use RayleighWaveKernel
  implicit none
  integer(c_int), value       :: nlayer, iflsph
  real(c_float),  intent(in)  :: thk(nlayer), vp(nlayer), vs(nlayer), rhom(nlayer)
  real(c_double), intent(inout) :: t, cp, cg
  real(c_double), intent(inout) :: dispu(nlayer), dispw(nlayer)
  real(c_double), intent(inout) :: stressu(nlayer), stressw(nlayer)
  real(c_double), intent(inout) :: dc2da(nlayer), dc2db(nlayer)
  real(c_double), intent(inout) :: dc2dh(nlayer), dc2dr(nlayer)

  call sregn96(thk, vp, vs, rhom, nlayer, t, cp, cg, &
               dispu, dispw, stressu, stressw, &
               dc2da, dc2db, dc2dh, dc2dr, iflsph)
end subroutine sregn96_c

subroutine sregn96_hti_c(thk, vp, vs, rhom, nlayer, t, cp, cg, &
                         dispu, dispw, stressu, stressw, &
                         dc2da, dc2db, dc2dh, dc2dr, dc2dgc, dc2dgs, iflsph) &
    bind(C, name="sregn96_hti_")
  use iso_c_binding
  use RayleighWaveKernel
  implicit none
  integer(c_int), value       :: nlayer, iflsph
  real(c_float),  intent(in)  :: thk(nlayer), vp(nlayer), vs(nlayer), rhom(nlayer)
  real(c_double), intent(inout) :: t, cp, cg
  real(c_double), intent(inout) :: dispu(nlayer), dispw(nlayer)
  real(c_double), intent(inout) :: stressu(nlayer), stressw(nlayer)
  real(c_double), intent(inout) :: dc2da(nlayer), dc2db(nlayer)
  real(c_double), intent(inout) :: dc2dh(nlayer), dc2dr(nlayer)
  real(c_double), intent(inout) :: dc2dgc(nlayer), dc2dgs(nlayer)

  call sregn96_hti(thk, vp, vs, rhom, nlayer, t, cp, cg, &
                   dispu, dispw, stressu, stressw, &
                   dc2da, dc2db, dc2dh, dc2dr, dc2dgc, dc2dgs, iflsph)
end subroutine sregn96_hti_c

subroutine sregnpu_c(thk, vp, vs, rhom, nlayer, t, cp, cg, &
                     dispu, dispw, stressu, stressw, &
                     t1, cp1, t2, cp2, &
                     dc2da, dc2db, dc2dh, dc2dr, &
                     du2da, du2db, du2dh, du2dr, &
                     iflsph) &
    bind(C, name="sregnpu_")
  use iso_c_binding
  use RayleighWaveKernel
  implicit none
  integer(c_int), value       :: nlayer, iflsph
  real(c_float),  intent(in)  :: thk(nlayer), vp(nlayer), vs(nlayer), rhom(nlayer)
  real(c_double), intent(inout) :: t, cp, cg
  real(c_double), intent(inout) :: dispu(nlayer), dispw(nlayer)
  real(c_double), intent(inout) :: stressu(nlayer), stressw(nlayer)
  real(c_double), intent(inout) :: t1, cp1, t2, cp2
  real(c_double), intent(inout) :: dc2da(nlayer), dc2db(nlayer)
  real(c_double), intent(inout) :: dc2dh(nlayer), dc2dr(nlayer)
  real(c_double), intent(inout) :: du2da(nlayer), du2db(nlayer)
  real(c_double), intent(inout) :: du2dh(nlayer), du2dr(nlayer)

  call sregnpu(thk, vp, vs, rhom, nlayer, t, cp, cg, &
               dispu, dispw, stressu, stressw, &
               t1, cp1, t2, cp2, &
               dc2da, dc2db, dc2dh, dc2dr, &
               du2da, du2db, du2dh, du2dr, &
               iflsph)
end subroutine sregnpu_c
