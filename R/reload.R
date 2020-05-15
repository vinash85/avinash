reload <- function(lib ){
   require(lib)
   try( 
      eval(parse(text=paste( "detach(package:", lib, ", unload = TRUE)", sep="")))
   )
   libr <- paste("package:",lib,sep="")
   try(library.dynam.unload( lib, system.file(package = lib)))
   require(lib, character.only = TRUE)
   
}
recompile <- function( path, lib, reload=T, Ccompile=F ){
   try( 
      eval(parse(text=paste( "detach(package:", lib, ", unload = TRUE)", sep="")))
   )
   libr <- paste("package:",lib,sep="")
   try(library.dynam.unload( lib, system.file(package = lib)))
   if(Ccompile) path <- paste( "--vanilla   CMD INSTALL -c --no-lock ", path )
   else path <- paste( "--vanilla   CMD INSTALL --no-lock ", path )
   system2( 'R', path  )
   if(reload) require(lib, character.only = TRUE)
}
recompile_r = function(path, lib.name, reload=T, ...){
   # tryCatch( eval(parse(text=paste( "detach(package:", lib.name, ", unload = TRUE)", sep=""))), error = function(e) print("not loaded... "))
   # lib.namer <- paste("package:",lib.name,sep="")
   # tryCatch(library.dynam.unload( lib.name, system.file(package = lib.name)), error = function(e) print("not loaded... "))
   remove.packages(lib.name,...)
   install.opts <- "--no-lock "
   install.packages(path, repos=NULL, clean = TRUE,INSTALL_opts = install.opts, ...) 
   if(reload) require(lib.name, character.only = TRUE)
}
# recompile_r_avinash = function()
   # recompile_r(lib.name ="avinash", path="~/shortcuts/avinash", lib = "/homes6/asahu/R/x86_64-pc-linux-gnu-library/3.6")
