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
recompile_r = function(path, lib, reload=T, Ccompile=F){
   eval(parse(text=paste( "detach(package:", lib, ", unload = TRUE)", sep="")))
   libr <- paste("package:",lib,sep="")
   try(library.dynam.unload( lib, system.file(package = lib)))
   if(Ccompile) install.opts <- "--no-lock " 
   else install.opts <- "--vanilla  --no-lock "
   # system2( 'R', path  )
   install.packages(path, repos=NULL, clean = TRUE) 
   if(reload) require(lib, character.only = TRUE)
}
