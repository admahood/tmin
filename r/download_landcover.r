# make yearly landcover mosaics from modis
library(httr)
library(tidyverse)
library(RCurl)


dir.create("data/MCD12Q1", recursive=TRUE)
# input your username and password here as character strings
uu <- "admahood" 
# pp <-

u_p <- paste0(uu,":",pp)
# define tiles
url1 <-paste0("https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2001.01.01/")
u_p_url <- paste0("https://", u_p, "@https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/2001.01.01/")

filenames <- RCurl::getURL(url1, userpwd = u_p)
# write to a temporary file
cat(filenames, file = 'tmp.txt')

# read the temporary file as a fixed width text file
dir_listing1 <- read_fwf('tmp.txt', fwf_empty('tmp.txt'),skip = 37) %>%
  dplyr::select(X1) %>%
  filter(str_detect(X1, "\\S{41}.hdf.xml"))%>%
  mutate(X1 = str_extract(X1, "\\S{41}.hdf"),
         tiles = substr(X1,18,23))
#world
tiles <- dir_listing1$tiles %>% unique
#coterminous US
# tiles <- c("h08v04","h09v04","h10v04","h11v04","h12v04","h13v04",
#           "h08v05","h09v05","h10v05","h11v05","h12v05",
#           "h08v06","h09v06","h10v06","h11v06")

tiles <- c("h13v14" ,"h14v14" ,"h12v13" ,"h13v13" ,"h11v12" ,"h12v12", "h13v12", "h08v11", "h11v11", "h12v11",
            "h13v11", "h14v11", "h10v10", "h11v10", "h12v10", "h13v10", "h14v10", "h08v09", "h09v09", "h10v09",
            "h11v09", "h12v09", "h13v09", "h14v09", "h08v08", "h09v08", "h10v08", "h11v08", "h12v08", "h13v08",
            "h03v07", "h07v07", "h08v07", "h09v07", "h10v07", "h11v07", "h12v07", "h03v06", "h07v06", "h08v06",
            "h09v06", "h10v06", "h11v06", "h07v05", "h08v05", "h09v05", "h10v05", "h11v05", "h12v05", "h08v04",
            "h09v04", "h10v04", "h11v04", "h12v04", "h13v04", "h14v04", "h06v03", "h07v03", "h08v03", "h09v03",
            "h10v03", "h11v03", "h12v03", "h13v03", "h14v03", "h15v03", "h28v03", "h29v03", "h09v02", "h10v02",
            "h11v02", "h12v02", "h13v02", "h14v02", "h15v02", "h16v02", "h17v02", "h12v01", "h13v01", "h14v01",
            "h15v01", "h16v01", "h17v01", "h16v00", "h17v00")

years = 2016:2019 #2017 disappeared!!

local_data <- "data/MCD12Q1"
s3_path <- "s3://earthlab-natem/modis-burned-area/input/landcover"

# corz <- detectCores()-1
# registerDoParallel(corz)
for (y in 1:length(years)) {
  
  
  url1 <-paste0("https://e4ftl01.cr.usgs.gov/MOTA/MCD12Q1.006/", years[y], ".01.01/")
  
  existing_files <- system(paste0("aws s3 ls ",
                                  file.path(s3_path, years[y]),"/"), intern = TRUE) %>%
    map(.,function(d) {
      xx <- data.frame(text = unlist(str_split(d, pattern = " ")))
      yy <- data.frame(file_size_bytes = xx[5,1], file_name = xx[6,1])
    }) %>%
    bind_rows()
  
  # "MCD12Q1.A2001001.h15v07.006.2018053201902.hdf"-> example
  # nchar(example)
  
  fl <- RCurl::getURL(url1, userpwd = u_p) %>%
    strsplit("> <")
  ff <- fl[[1]][str_detect(fl[[1]], "\\S{41}.hdf<")]
  fdf <- data.frame(file_name = str_extract(ff, "\\S{41}.hdf"),
                    file_size_megabytes = str_extract(ff, "\\d{2}M")) %>%
    mutate(file_size_bytes = substr(file_size_megabytes,1,2)%>%as.numeric()*1000000) %>%
    filter(str_extract(file_name, "h\\d{2}v\\d{2}") %in% tiles,
           !file_name %in% existing_files$file_name) %>%
    mutate(file_name = as.character(file_name))
  
  if(nrow(existing_files) != length(tiles)){
    dir.create(file.path(local_data, years[y]), showWarnings = T, recursive = T)
    for(i in fdf$file_name){
      t0<-Sys.time()
      output_file_name <- file.path(local_data,years[y],i)
      
      if(!file.exists(output_file_name)){
        httr::GET(paste0(url1,i),
                  authenticate(uu,pp),
                  write_disk(path = output_file_name, overwrite = TRUE))
        # download.file(paste0(u_p_url, "/",years[y], "/", output_file_name),
        #               output_file_name, method = "wget", quiet = TRUE)
      }
      
      if(file.info(output_file_name)$size <= 6000000){
        httr::GET(paste0(url1,i),
                  authenticate(uu,pp),
                  write_disk(path = output_file_name, overwrite=TRUE))
      }
    }
    
    system(paste(
      "aws s3 sync",
      file.path(local_data,years[y]),
      file.path("s3://earthlab-natem/modis-burned-area/input/landcover",years[y])
    ))
    unlink(file.path(local_data, years[y]), recursive = TRUE)
    print(Sys.time() - t0)
    print(years[y])
  }
}

