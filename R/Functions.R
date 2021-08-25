#' This function takes the summary results from RBPMap and converts it to a list of lists
#'
#' @param path.to.file The path of the RBPMap text file
#' @return the resulting list of lists summarizing RBPMap results
#'
inputRBP <- function(path.to.file){

  inputed.data = readLines(path.to.file)

  asterik.lines = grep(inputed.data, pattern = "[*]")

  asterik.split = list()

  chr.names = c()

  for(i in 2:length(asterik.lines)){
    asterik.split[[i-1]] = inputed.data[asterik.lines[i-1]:asterik.lines[i]]
    chr.names[i-1] = asterik.split[[i-1]][ grep(asterik.split[[i-1]],pattern = "[=]")-1]
  }

  names(asterik.split) = chr.names


  RBP.summary = list()
  for(j in 1:length(asterik.split)){

    if(length(grep(asterik.split[[j]],pattern = "No motifs found.")) == 0 & length(grep(asterik.split[[j]], pattern = "Error")) == 0){
      select.asterik.split = asterik.split[[j]]
      protein.heading.index = grep(select.asterik.split, pattern = "Protein")
      blank.space = which(select.asterik.split == "")
      start.end.df = data.frame('start'=protein.heading.index,
                                'end'=blank.space[3:length(blank.space)])
      start.end.df$start = start.end.df$start + 1
      start.end.df$end = start.end.df$end - 1

      protein.list = list()
      protein.df.list = list()
      for(i in 1:nrow(start.end.df)){
        protein.list[[i]] = select.asterik.split[start.end.df$start[i]:start.end.df$end[i]]
        writeLines(protein.list[[i]],con = "temp.txt")
        protein.df.list[[i]] = read.delim("temp.txt")
      }

      names(protein.list) = select.asterik.split[protein.heading.index]
      names(protein.df.list) = select.asterik.split[protein.heading.index]
      RBP.summary[[j]] = protein.df.list

    }else{
      RBP.summary[[j]] = NULL
      print(paste(j, "is null"))
    }
  }


  names(RBP.summary) = chr.names

  inputRBP = RBP.summary

  return(inputRBP)
}

#' This function takes the information on particular gene locations and creates a file to search for RBPs on RBPMap
#'
#' @param chr.vector Chromosome locations
#' @param start.vector Start positions
#' @param end.vector End positions
#' @param strand.vector Strand positions
#' @param bprange BP range
#' @path.to.output output file
#'
convert_to_RBP <- function(chr.vector, start.vector, end.vector, strand.vector, bprange=100,
                           path.to.output = "output.txt"){

  rbp.vector = c()
  for(i in 1:length(chr.vector)){
    rbp.vector[i] = paste(chr.vector[i],
                          ":",
                          start.vector[i]-bprange/2,
                          "-",
                          end.vector[i]+bprange/2,
                          ":", strand.vector[i],sep="")
  }
  writeLines(rbp.vector, con=path.to.output)
}


#' This function samples a list from an RBPMap results file
#'
#' @param path.to.file The path of the RBPMap text file
#' @param samplecount Number of samples needed
#' @return the resulting list of lists summarizing RBPMap results
#'
sampleRBP <- function(path.to.file, samplecount){

  inputed.data = readLines(path.to.file)

  asterik.lines = grep(inputed.data, pattern = "[*]")

  asterik.split = list()

  chr.names = c()

  for(i in 2:length(asterik.lines)){
    asterik.split[[i-1]] = inputed.data[asterik.lines[i-1]:asterik.lines[i]]
    chr.names[i-1] = asterik.split[[i-1]][ grep(asterik.split[[i-1]],pattern = "[=]")-1]
  }

  names(asterik.split) = chr.names

  sample.x = sample(x=1:length(asterik.split),size = samplecount,
         replace = F)

  print(paste("Sampling", samplecount, "regions"))
  asterik.split = asterik.split[sample.x]
  chr.names = chr.names[sample.x]

  RBP.summary = list()
  for(j in 1:length(asterik.split)){

    if(length(grep(asterik.split[[j]],pattern = "No motifs found.")) == 0 & length(grep(asterik.split[[j]], pattern = "Error")) == 0){
      select.asterik.split = asterik.split[[j]]
      protein.heading.index = grep(select.asterik.split, pattern = "Protein")
      blank.space = which(select.asterik.split == "")
      start.end.df = data.frame('start'=protein.heading.index,
                                'end'=blank.space[3:length(blank.space)])
      start.end.df$start = start.end.df$start + 1
      start.end.df$end = start.end.df$end - 1

      protein.list = list()
      protein.df.list = list()
      for(i in 1:nrow(start.end.df)){
        protein.list[[i]] = select.asterik.split[start.end.df$start[i]:start.end.df$end[i]]
        writeLines(protein.list[[i]],con = "temp.txt")
        protein.df.list[[i]] = read.delim("temp.txt")
      }

      names(protein.list) = select.asterik.split[protein.heading.index]
      names(protein.df.list) = select.asterik.split[protein.heading.index]
      RBP.summary[[j]] = protein.df.list

    }else{
      RBP.summary[[j]] = NULL
      print(paste(j, "is null"))
    }
  }


  names(RBP.summary) = chr.names

  inputRBP = RBP.summary

  return(inputRBP)
}
