program flreg
  use precision
  use string
  use fileio
  use simulation
  use parser
  use partition
  implicit none
  character(len=str_medium) :: input_name, tmp_string

  ! Initialize parallel environment
  call parallel_init

  ! Initialize the random number generator
  call random_init

  ! Parse the command line
  call parallel_get_inputname(input_name)

  ! Parse the input file
  call parser_init
  call parser_parsefile(input_name)

  ! ! Geometry initialization
  ! call geometry_init

  ! ! Data initialization
  ! call data_init
  ! call optdata_init




end program flreg
