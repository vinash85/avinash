intersectBed <- function(functionstring="intersectBed", a, b,opt.string="")
{
  #options(scipen =99) # not to use scientific notation when writing out
  #create temp files
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
  out   =tempfile()
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))

  res=fread(out)
  if(!is.character(a)) unlink(a.file)
  if(!is.character(b))   unlink(b.file);
  
  unlink(out)
  return(res)
}

