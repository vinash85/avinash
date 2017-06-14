intersectBed <- function(functionstring="intersectBed",  a, b, outfile=NULL, opt.string="")
{
  options(scipen =99) # not to use scientific notation when writing out
  #create temp files
  #Sys.umask(mode="0") #write permission for file
  save.out <- F
  if(!is.null(outfile)){
    out = outfile
    file.remove(out)
    save.out <- T
  }else{
    out   =tempfile()
  }
  file.create(out)
  Sys.chmod(out)
  if(!is.character(a)){
  a.file=tempfile()
  write.table(a ,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  }else{
    a.file <- a
  }
  if(!is.character(b)){
    b.file=tempfile()
    write.table(b ,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  }else{
    b.file <- b
  }
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))

  res=fread(out)
  if(!is.character(a)) unlink(a.file)
  if(!is.character(b))   unlink(b.file);
  
  if(!save.out) unlink(out)
  return(res)
}

