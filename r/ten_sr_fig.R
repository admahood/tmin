# figure for ten simple rules paper
library(tidyverse)
library(sf)
library(raster)

# system("aws s3 cp s3://earthlab-amahood/na_files/tables/fired_events_s1_t5_2020153.csv data/na_fired.csv")

# setting up the data ================================
template_path <- "template/template.tif"
template <- raster(template_path)

# system("aws s3 cp s3://earthlab-amahood/na_files/tables/fired_events_s1_t5_2020153.csv data/na_fired.csv")
raw_events_file <- "data/na_fired.csv"


# loading in fire event data frame
df <- read_csv(raw_events_file) %>%
  dplyr::select(id,date,x,y) %>%
  #centering the pixels on the raster cells
  mutate(x = x + (res(template)[1]/2),
         y = y - (res(template)[1]/2),
         year = as.numeric(substr(date, 1,4))) %>%
  # removing about 8000 repeat pixels (i.e. adjacent month detections)
  distinct(x,y,id, .keep_all = T)


# first, serial =================================================
set.seed(2020)
n_events<- c(100, 300, 1000, 3000, 10000)
ids <- df$id %>% unique()

result_df<- data.frame(rep=NA, n_events = NA, time=NA)
counter<-1
for(i in 1:length(n_events)){
  x<- microbenchmark(
    df_poly <- df %>%
      filter(id %in% sample(x = ids, replace = FALSE, size = n_events[i]))%>%
      st_as_sf(coords = c("x","y"), crs = crs(template, asText=TRUE)) %>%
      group_by(id)%>%
      st_buffer(dist = 1+(res(template)[1]/2), endCapStyle = "SQUARE")%>%
      dplyr::summarize(),
    times = 10, unit = "s"
   )
    
  for(j in 1:length(x$time)){
    print(paste(n_events[i], x$time[j]))
    result_df[counter, 1] <- j
    result_df[counter, 2] <- n_events[i]
    result_df[counter, 3] <- x$time[j]
    counter<- counter + 1}
}
result_df$cat <- "serial"

# parallel optimized =====================================
library(foreach);library(doParallel)
doParallel::registerDoParallel(parallel::detectCores()-2)
rep<-1:10
# ids <- unique(df$id)
result_dfp<- data.frame(rep=NA, n_events = NA, time=NA, cat=NA)
counter<-1
for(i in 1:length(n_events)){
  for(j in rep){
  df_grps <- df %>%
    filter(id %in% sample(x = ids, replace = FALSE, size = n_events[i]))%>%
    mutate(grp = cut(id, breaks=6))
  
  grps <- unique(df_grps$grp)
  
  x<- microbenchmark(
          df_poly <- foreach(g = grps) %dopar% {

          zzz <- df_grps %>%
            filter(grp == g) %>%
            group_by(id) %>%
            st_as_sf(coords = c("x","y"), crs = crs(template, asText=TRUE)) %>%
            st_buffer(dist = 1+(res(template)[1]/2), endCapStyle = "SQUARE")%>%
            dplyr::summarize() 
          return(zzz)
          
        },
        bind_rows(df_poly), times=1
        )

    print(paste(n_events[i], x$time))
    result_dfp[counter, 1] <- j
    result_dfp[counter, 2] <- n_events[i]
    result_dfp[counter, 3] <- sum(x$time)
    result_dfp[counter, 4] <- "par_op"
    counter<- counter + 1
}}
rbind(result_df, result_dfp) %>%
  mutate(cat = ifelse(cat == "serial", 
                      "Minimally Optimized", 
                      "Fully Optimized")) %>%
ggplot(aes(n_events, time/10000000, color = cat)) +
  # geom_smooth(se = F) +
  geom_jitter(width = .1, size=1, alpha=0.7) +
  scale_x_log10(breaks = c(100,300,1000,3000,10000)) +
  theme_classic() +
  xlab("Number of Groups") +
  ylab("Computation Time")+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        axis.text.y = element_blank(), 
        axis.line = element_line(arrow = arrow(type="closed", angle = 20,
                                               length=unit(0.15, "inches"))),
        axis.ticks = element_blank()) +
  ggsave("tsr_fig1b.png", width =3.5, height=5)
