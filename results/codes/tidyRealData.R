library(readxl)
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Read individual years
listFiles <- list.files("Real data", pattern = "*.RData", full.names = TRUE)
readFiles <- lapply(listFiles, loadRData)
str(readFiles)

# Read data along time

ireland <- read_excel("Real data/ireland.xlsx",
                             sheet = "Optimal")
str(ireland)

# ireland <-  ireland %>%
#   select("Year", "Genotype", "Location", "Bloc", "yld_ton_ha") %>%
#   rename(Yield = yld_ton_ha, Environment = Location) %>%
#   mutate_at("Yield", as.numeric) %>%
#   mutate_if(is.character, as.factor)  %>%
#   group_by(Genotype, Environment, Year) %>%
#   mutate(Mean = mean(Yield))  %>%
#   summarise(Mean = sort(unique(Mean)))


ireland <-  ireland |>
  select("Year", "Genotype", "Location", "Bloc", "yld_ton_ha") |>
  rename(Yield = yld_ton_ha, Environment = Location) |>
  mutate_at("Yield", as.numeric) |>
  mutate_at(c("Year", "Genotype", "Environment"), as.factor)
ireland$Bloc <-  substr(ireland$Bloc, start = -7)

# Rename factor levels  to be anonymous
levels(ireland$Genotype) <- paste0("g", 1:length(levels(ireland$Genotype)))
levels(ireland$Environment) <- paste0("e", 1:length(levels(ireland$Environment)))

# add env to block name and turns it into factor
ireland$Bloc <- as.factor(paste0(ireland$Environment, ireland$Bloc))

str(ireland)

save(ireland, file = "Real data/ireland.RData")
#saveRDS(ireland, file = "Real data/ireland.RDS")

abc <- ireland |> filter(Year == "2010")

# separate the data into train and test

train <- subset(ireland,  grepl('1$|2$', Bloc))
test <- subset(ireland,  grepl('3$|4$', Bloc))


save(train, file = "Real data/train_ireland.RData")
save(test, file = "Real data/test_ireland.RData")

# Data set 2 --------------------------------------------------------------

data2 <- read_csv("Real data/data2.csv")
data2$YEAR <- as.factor(data2$YEAR)
data2$Entry_Number <- as.factor(data2$Entry_Number)
data2$LOCATION <- as.factor(data2$LOCATION)
data2$YIELD <- as.numeric(data2$YIELD)
data2$BLOCK <- as.factor(data2$BLOCK)

data2 <-  data2 %>%
  rename(GENOTYPE = Entry_Number) %>%
  group_by(GENOTYPE, LOCATION, YEAR) %>%
  mutate(MEAN = mean(YIELD))  %>%
  summarise(MEAN = sort(unique(MEAN)))

str(data2)
# Rename factor levels  to be anonymous
levels(data2$GENOTYPE) <- paste0("g", 1:length(levels(data2$GENOTYPE)))
levels(data2$LOCATION) <- paste0("e", 1:length(levels(data2$LOCATION)))

save(data2, file = "Real data/data2.RData")


