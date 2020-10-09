!     integer, parameter  :: dp = kind(1.0d0)
      real(dp) :: PI

      integer  :: el1_nbeg, el1_nend
      integer  :: el2_nbeg, el2_nend
      integer  :: el3_nbeg, el3_nend
      integer  :: out_nbeg, out_nend

      common/defpar/pi,el1_nbeg,el1_nend,el2_nbeg,el2_nend,el3_nbeg, &
                    el3_nend,out_nbeg, out_nend
