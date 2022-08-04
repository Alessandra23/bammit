library(dplyr)
complete_vcu <- readRDS("~/Documents/GitHub/bammit/Real data/vcus/complete_vcu.rds")
ws <- readRDS("~/Documents/GitHub/bammit/Real data/vcus/ws.rds")
ws2020 <- ws |> dplyr:: filter(YEAR == "2020") #|>
  #select(Location, winter_spring, ACZ)

vcu <- complete_vcu |> select(`Plot fresh yield`, drought_core,
                              ACZ, ENV, Experiment, Species)

ws2020 <- ws2020 |> filter(Location %in% complete_vcu$ENV) |>
  select(Location, ALT, LON, LAT, ACZ, Country, VCU_CODE)
vcu$Location <- vcu$ENV

names(vcu)
names(ws2020)
df11 <- merge(vcu, ws2020, by = c("Location", "ACZ"))

df12 <- unique(df11)




# Create a small data set -------------------------------------------------

set.seed(1983)
posws2020 <- sample(1:nrow(ws2020), 5, replace = TRUE)
posws2020 <- ws2020[posws2020,]
posws2020 <- posws2020[-1,]

posvcu <- complete_vcu |> filter(ENV %in% posws2020$Location)
set.seed(1983)
posvcurows <- sample(1:nrow(posvcu), 5, replace = TRUE)
posvcu <- posvcu[posvcurows,] |> rename(Location = ENV)




dfm <- merge(posws2020,posvcu, by = c("Location", "ACZ") )


View(dfm[c(1,2),])

# Tests -------------------------------------------------------------------






df <- left_join(complete_vcu, ws2020) |> head()

complete_vcu |> names()
ws |> names()

which(!(complete_vcu$ENV %in% as.factor(ws$Location)))

ws |> group_by(Location) |> summarise(c = length(Location))
complete_vcu |> group_by(ENV) |> summarise(c = length(ENV))



ws

aa <- complete_vcu |> filter(Species == "Bread")
aa$ENV |> levels()
ws2020 |> filter()

names(ws2020)


cvcu2020 <- complete_vcu |> select(Yield, Experiment, Species, ACZ, GID, ENV)

cvcu2020$ENV |> levels()
ws2020$Location |> as.factor() |> levels()

complete_vcu$Species

merge(cvcu2020, ws2020)

dim(complete_vcu)
dim(ws)

complete_vcu$Experiment
complete_vcu$Species
complete_vcu$Treatment
complete_vcu$Replicate
complete_vcu$Plot |> str()
complete_vcu$`Growth Habit` |> str()

complete_vcu$Emergence
complete_vcu$Ripening
complete_vcu$`Harvest Date`
complete_vcu$YL_N_diseases |> quantile()
complete_vcu$Pow_AUDPC
complete_vcu$Pow_N_diseases
complete_vcu$drought_core


complete_vcu$ENV
ws$env
ws$Location
ws$YYYYMMDD
ws$env


complete_vcu$`Year of Sowing`
ws$YEAR |> as.factor()

ws$Location
ws$local
ws$daysFromStart
ws$Altitiude |> str()
all.equal(ws$ALT, ws$Altitiude)

ws$Country |> as.factor() |> levels()
ws$

  data.frame(complete_vcu$Experiment, complete_vcu$drought_core)


merge(complete_vcu, ws)

colnames(complete_vcu)
colnames(ws)
head(ws)
head(complete_vcu)


# fake df -----------------------------------------------------------------

df1 <- data.frame(acz = as.vector(replicate(3, sample(c("A", "B", "C", "D")))))
