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
