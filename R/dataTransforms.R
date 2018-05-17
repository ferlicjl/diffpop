#' eventsToCSV
#'
#' convert .events file to CSV
#'
#' @param directory directory of events file
#' @param prefixes a list of prefixes to convert
#' @param size initial population size
#' @param label initial probabilty that a cell receives a unique barcode
#'
#' @export
eventsToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".events$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".events$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = T)
    write.csv(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F)
  }
}

#' labelToCSV
#'
#' convert .label file to CSV
#'
#' @param directory directory of label file
#' @param prefixes a list of prefixes to convert
#'
#' @export
labelToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".label$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".label$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = T)
    write.csv(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F)
  }
}

#' popToCSV
#'
#' convert .pop file to CSV
#'
#' @param directory directory of pop file
#' @param prefixes a list of prefixes to convert
#'
#' @export
popToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".pop$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".pop$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = T)
    write.csv(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F)
  }
}

#' sdiToCSV
#'
#' convert sdi file to CSV
#'
#' @param directory directory of sdi file
#' @param prefixes a list of prefixes to convert
#'
#' @export
sdiToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".sdi$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".sdi$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = T)
    write.csv(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F)
  }
}

#' typesToCSV
#'
#' convert .types file to CSV
#'
#' @param directory directory of types file
#' @param prefixes a list of prefixes to convert
#'
#' @export
typesToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".types$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".types$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = T)
    write.csv(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F)
  }
}

#' allToCSV
#'
#' convert .all file to CSV
#'
#' @param directory directory of types file
#' @param prefixes a list of prefixes to convert
#'
#' @export
allToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".all$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".*.all$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    tdat = read.delim(f, header = F, sep = " ")
    tdat = tdat[1:(nrow(tdat)-4),c(2, 4, 6, 8)]
    names(tdat) = c("barcode", "mutation", "fitness", "count")
    write.csv(tdat, paste(name, ".csv", sep = ""), row.names = F, quote = F)
  }
}

#' paramsToCSV
#'
#' convert params file to CSV
#'
#' @param directory directory of params file
#' @param prefixes a list of prefixes to convert
#'
#' @export
paramsToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".params$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".params$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = F)
    write.table(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F, col.names = F, sep = ",")
  }
}

#' params2ToCSV
#'
#' convert params2 file to CSV
#'
#' @param directory directory of params2 file
#' @param prefixes a list of prefixes to convert
#'
#' @export
params2ToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".params2$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".params2$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = F)
    write.table(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F, col.names = F, sep = ",")
  }
}

#' mutToCSV
#'
#' convert mut file to CSV
#'
#' @param directory directory of mut file
#' @param prefixes a list of prefixes to convert
#'
#' @export
mutToCSV = function(directory = ".", prefixes = NULL){
  dir = R.utils::getAbsolutePath(directory)

  files = c()
  if(length(prefixes) == 0)
    files = grep(paste(".mut$", sep = ""), list.files(dir, full.names = T),value=TRUE)
  else{
    for(p in prefixes){
      files = c(files, grep(paste(p, ".mut$", sep = ""), list.files(dir, full.names = T),value=TRUE) )
    }
  }

  print(files)
  for(f in files){
    name = gsub("\\.", "_", f)
    t_file = read.table(f, header = F)
    write.table(t_file, paste(name, ".csv", sep = ""), row.names = F, quote = F, col.names = F, sep = ",")
  }
}

#' outputToCSV
#'
#' convert all output files to CSV
#'
#' @param directory directory of output files
#' @param prefixes a list of prefixes to convert
#'
#' @export
outputToCSV = function(directory = ".", prefixes = NULL){
  eventsToCSV(directory, prefixes)
  labelToCSV(directory, prefixes)
  mutToCSV(directory, prefixes)
  popToCSV(directory, prefixes)
  sdiToCSV(directory, prefixes)
  typesToCSV(directory, prefixes)
  paramsToCSV(directory, prefixes)
  params2ToCSV(directory, prefixes)
  allToCSV(directory, prefixes)
}


