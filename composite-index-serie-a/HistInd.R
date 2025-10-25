# clearing memory
rm(list=ls())
gc()

# setting the seed
#set.seed(666)

# loading packages
library(dplyr)
library(ggplot2)
library(Compind)
library(reshape2)
library(ggimage)
library(png)


sPath = "C:/Users/frank/OneDrive/Desktop/Composite indicators/"

Ata <- paste0(sPath, "Ata.png")
Bol <- paste0(sPath, "Bol.png")
Fio <- paste0(sPath, "Fio.png")
Int <- paste0(sPath, "Int.png")
Juv <- paste0(sPath, "Juv.png")
Laz <- paste0(sPath, "Laz.png")
Mil <- paste0(sPath, "Mil.png")
Nap <- paste0(sPath, "Nap.png")
Rom <- paste0(sPath, "Rom.png")
Sas <- paste0(sPath, "Sas.png")
Tor <- paste0(sPath, "Tor.png")
Udi <- paste0(sPath, "Udi.png")

Images <- c(rep(Ata,8),rep(Bol,8),rep(Fio,8),rep(Int,8),rep(Juv,8),rep(Laz,8),rep(Mil,8),rep(Nap,8),rep(Rom,8),rep(Sas,8),rep(Tor,8),rep(Udi,8))
Images

img <- png::readPNG(Images)

df <- cbind(df, Images)

# loading data
df <- read.csv("C:/Users/frank/OneDrive/Desktop/CompInd.csv")
summary(df)
names(df)

table(df$Year,df$Team)

#---------------------------------------------------------------
# HISTORICAL INDEX
#---------------------------------------------------------------

df1 <- df %>%
  select(Year, Team, Position, PPG, DR, CupWins, EuropePart, Trophies, Ncoach)

# compute goalposts
goals1 <- apply(df1[df1$Year=="2016/17",c(3:9)], 2, median)

#df1_norm <- normalise_ci(df1, c(3:9), c("NEG","POS","POS","POS","POS","POS","NEG"), method=1, z.mean=100, z.std=10)
#df1_norm

# AMPI
CI1 <- ci_ampi(df1, 
              indic_col=c(3:9),
              gp=goals1,
              time=df1$Year, 
              polarity= c("NEG","POS","POS","POS","POS","POS","NEG"), 
              penalty="POS")

# score
CI1$ci_ampi_est
# penalties
CI1$ci_penalty
# normalized values
CI1$ci_norm

# attaching scores
df1_finale <- cbind(unique(df1$Team),CI1$ci_ampi_est)
names(df1_finale)[1] <- "Team"

# plotting trends
df1_finale <- melt(df1_finale, id="Team")
ggplot(df1_finale, aes(x=variable, y=value, color=as.factor(Team), group=as.factor(Team)))+
  geom_line()+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=15)) +
  labs(x="Years",
       y="Historical Score",
       color="Team")

#---------------------------------------------------------------
# DEFENSIVE INDEX
#---------------------------------------------------------------

df2 <- df %>%
  select(Year, Team, XGS, GS, PercPar, PortInv, Amm, Fouls, Int, TackW)

# compute goalposts
goals2 <- apply(df2[df2$Year=="2016/17",c(3:10)], 2, median)

# AMPI
CI2 <- ci_ampi(df2, 
              indic_col=c(3:10),
              gp=goals2,
              time=df2$Year, 
              polarity= c("POS","NEG","POS","POS","NEG","NEG","POS","POS"), 
              penalty="NEG")

# score
CI2$ci_ampi_est
# penalties
CI2$ci_penalty
# normalized values
CI2$ci_norm

# attaching scores
df2_finale <- cbind(unique(df2$Team),CI2$ci_ampi_est)
names(df2_finale)[1] <- "Team"

# plotting trends
df2_finale <- melt(df2_finale, id="Team")
ggplot(df2_finale, aes(x=variable, y=value, color=as.factor(Team), group=as.factor(Team)))+
  geom_line()+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=15)) +
  labs(x="Years",
       y="Defensive Score",
       color="Team")

#---------------------------------------------------------------
# OFFENSIVE INDEX
#---------------------------------------------------------------

df3 <- df %>%
  select(Year, Team, Bposs, XG, GF, Shoot, PassCompl, Dribbling)

# compute goalposts
goals3 <- apply(df3[df3$Year=="2016/17",c(3:8)], 2, median)

# AMPI
CI3 <- ci_ampi(df3, 
               indic_col=c(3:8),
               gp=goals3,
               time=df3$Year, 
               polarity= c("POS","NEG","POS","POS","POS","POS"), 
               penalty="POS")

# score
CI3$ci_ampi_est
# penalties
CI3$ci_penalty
# normalized values
CI3$ci_norm

# attaching scores
df3_finale <- cbind(unique(df3$Team),CI3$ci_ampi_est)
names(df3_finale)[1] <- "Team"

# plotting trends
df3_finale <- melt(df3_finale, id="Team")
ggplot(df3_finale, aes(x=variable, y=value, color=as.factor(Team), group=as.factor(Team)))+
  geom_line()+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=15)) +
  labs(x="Years",
       y="Offensive Score",
       color="Team")


#---------------------------------------------------------------
# MARKET INDEX
#---------------------------------------------------------------

#df4 <- dfm %>%
 # select(YEAR, SQUAD, VALTEAM, VALKEYPLAYERS, DEBIT, TRANSF, STADIUM, STPROP, OLD, MANAGER) %>%  #age is neg, wage neg, over pos (penality)
  #filter(YEAR!=2024)

df4 <- dfc %>%
  select(YEAR, SQUAD, DIRITTI, COSTOROSA, VALKEYPLAYERS, DEBIT, TRANSF, STADIUM, OLD, MANAGER)
  #filter(YEAR != 2024)

# compute goalposts
goals4 <- apply(df4[df4$YEAR=="2016/17",c(3:10)], 2, median)

# AMPI
CI4 <- ci_ampi(df4, 
               indic_col=c(3:10),
               gp=goals4,
               time=df4$YEAR, 
               polarity= c("POS","POS","POS","NEG","POS","POS","NEG","POS"), 
               penalty="NEG")

# score
CI4$ci_ampi_est
# penalties
CI4$ci_penalty
# normalized values
CI4$ci_norm

# attaching scores
df4_finale <- cbind(unique(df4$SQUAD),CI4$ci_ampi_est)
names(df4_finale)[1] <- "Team"

# plotting trends
df4_finale <- melt(df4_finale, id="Team")
ggplot(df4_finale, aes(x=variable, y=value, color=as.factor(Team), group=as.factor(Team)))+
  geom_line()+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=15)) +
  labs(x="Years",
       y="Market Score",
       color="Team")

#---------------------------------------

CI <- CI1$ci_ampi_est*0.3 + CI2$ci_ampi_est*0.3 + CI3$ci_ampi_est*0.3 + CI4$ci_ampi_est*0.2
df_finale <- cbind(unique(df$Team),CI)
names(df_finale)[1] <- "Team"
names(df_finale)[2] <- "Year"

df_finale <- melt(df_finale, id="Team")
ggplot(df_finale, aes(x=variable, y=value, color=as.factor(Team), group=as.factor(Team)))+
  #geom_image(aes(image="https://banner2.cleanpng.com/20180628/bss/kisspng-atalanta-b-c-201718-serie-a-juventus-stadium-s-5b34e43537bfd1.7463974615301929492284.jpg"), size=.05)
  geom_line(size=1.1)+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=15)) +
  labs(x="Years",
       y="TCCI",
       color="Team")

names(df4_finale)[2] <- "Year"
dfg <- df4_finale %>%
  filter(Year=="2023/24")

head(dfg[order(dfg,decreasing = T),],10)

df_finale <- cbind(df_finale, Images)

ggplot(df_finale, aes(x=variable, y=value, image = png::readPNG(Ata), color=as.factor(Team), group=as.factor(Team)))+
  geom_image(size=0.05)+
  geom_line(size=1.1)+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=15)) +
  labs(x="Years",
       y="TCCI",
       color="Team")
