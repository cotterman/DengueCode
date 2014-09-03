#############################################################################
############### test code to work with python's rpy package  ################
#############################################################################


run_test_wrap = function(x_from_python){
  cubed = x_from_python^3
  print("let's test this")
  return(cubed)
}
#run_test_wrap(2)
