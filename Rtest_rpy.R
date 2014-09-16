#############################################################################
############### test code to work with python's rpy package  ################
#############################################################################

convert_factor_to_numeric = function(data, vars_to_convert){
  #print(vars_to_convert)
  for (name in vars_to_convert) {
    if(name=="Sexo"){
      data[,name] = as.numeric(data$Sexo)-1
    }
    data[,name] = as.numeric(as.character(data[,name]))
  } 
  return(data)
}


run_test_wrap = function(x_from_python){
  
  print("Variable types in beginning of R program")
  print( sapply(x_from_python , class) )
  
  #obtain list of intensity variables
  Xindices = grep("X",colnames(x_from_python)) #column indices in dataset
  intensity_vars = colnames(x_from_python[Xindices]) #names
  #convert these variables to numerics
  x_from_python = convert_factor_to_numeric(x_from_python, intensity_vars)
  
  #make sure variable types are as expected
  print("Variable types right before run_predictions function")
  print( sapply(x_from_python , class) )
  
  
}
#run_test_wrap(2)
