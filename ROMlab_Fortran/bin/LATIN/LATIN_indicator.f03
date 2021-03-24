module LATIN_indicator


  use Structure
  use Matrix_op
  use global_stage
  use local_stage

  implicit none

  contains


  function indicator1(G,L)  result(ind)

    !type(Input), intent(inout) :: I
    type(global), intent(inout) :: G
    type(local), intent(inout) :: L


    double precision :: ind
    ind =  100*( norm2( L%d_Ep - G%d_Ep )/norm2( L%d_Ep ) + norm2( L%sigma - G%sigma )/norm2( L%sigma ) )


  end function indicator1





end module LATIN_indicator
