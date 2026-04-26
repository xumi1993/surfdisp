! C-callable wrapper for the original Fortran surfdisp96 subroutine.
! Uses Fortran 2003 bind(C) and value attribute so integers are passed by value,
! matching the C function signature in surfdisp.h.
subroutine surfdisp96_f(thkm, vpm, vsm, rhom, nlayer, iflsph, iwave, &
                         mode, igr, kmax, t, cg) bind(C, name="surfdisp96_f")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int), value, intent(in)  :: nlayer, iflsph, iwave, mode, igr, kmax
  real(c_float),         intent(in)  :: thkm(nlayer), vpm(nlayer)
  real(c_float),         intent(in)  :: vsm(nlayer),  rhom(nlayer)
  real(c_double),        intent(in)  :: t(kmax)
  real(c_double),        intent(inout) :: cg(kmax)

  ! Local copies because surfdisp96 declares intent(in) arrays separately
  real(c_float)  :: thkm_l(nlayer), vpm_l(nlayer), vsm_l(nlayer), rhom_l(nlayer)
  real(c_double) :: t_l(kmax), cg_l(kmax)
  integer(c_int) :: nl, ifl, iw, md, ig, km

  nl  = nlayer
  ifl = iflsph
  iw  = iwave
  md  = mode
  ig  = igr
  km  = kmax
  thkm_l = thkm
  vpm_l  = vpm
  vsm_l  = vsm
  rhom_l = rhom
  t_l    = t
  cg_l   = 0.0d0

  call surfdisp96(thkm_l, vpm_l, vsm_l, rhom_l, nl, ifl, iw, md, ig, km, t_l, cg_l)

  cg = cg_l
end subroutine surfdisp96_f
