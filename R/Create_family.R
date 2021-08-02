
#' @title  createShell
#' @author Dieter Henrik Heiland
#' @description createShell
#' @param string The string contains a function that need to be run in a external shell
#' @param object.to.load The object that need to be processed (the function will save the object) ... name.obj ....
#' @param wd The folder in which the object will be saved
#' @param nameFromObject The name of the object in your environment as character ... "name.obj"...
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export


createShell <- function(string, object.to.load, nameFromObject, wd=getwd()){
  string.start <- ' #!/usr/bin/Rscript \n \n \n '
  
  # Save temp env.
  base::setwd(wd)
  base::saveRDS(object.to.load, "temp.obj.RDS")
  
  #Add used environment
  add <- base::paste0(nameFromObject,' <- readRDS("temp.obj.RDS") \n')
  
  #Add wd
  wd <- base::paste0("setwd(", " '", wd ,"'", " ) \n " )
  
  save <- base::paste0(' saveRDS(', nameFromObject , ', ', nameFromObject,'.RDS ) \n')
  
  string.out <- base::paste0(string.start,
                       wd,
                       add, 
                       string,
                       save,
                       '\n \n print("Done") \n \n'
  )
  
  
  base::writeLines(string.out, "shell.R")
}





