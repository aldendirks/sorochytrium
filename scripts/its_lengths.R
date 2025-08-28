# START

# Load libraries
library(tidyverse)
library(ggrepel)
library(hrbrthemes)
library(ggridges)
library(ggExtra)
library(geomtextpath)

################################################################################

# PREPARE DATA

# ITSoneDB
# Load data
itsonedb.ena <- read_tsv("ITSoneDB_lengths_ENA.txt", col_names = TRUE)
itsonedb.hmm <- read_tsv("ITSoneDB_lengths_HMM.txt", col_names = TRUE)
itsonedb.both <- read_tsv("ITSoneDB_lengths_both.txt", col_names = TRUE) 
itsonedb.rep <- read_tsv("ITSoneDB_lengths_rep.txt", col_names = TRUE)
itsonedb.infernal <- read_tsv("ITSoneDB_lengths_infernal.tsv", col_names = c("Accession", "ITS1", "Five", "ITS2", "Full", "Notes"), skip=1) # Only processed rep with Infernal
itsonedb.infernal <- itsonedb.infernal[itsonedb.infernal$ITS1 > 0,]
# Correct rep df with Infernal
itsonedb.rep <- itsonedb.rep %>%
  left_join(itsonedb.infernal %>% select(Accession, ITS1), by = "Accession") %>%
  mutate(Length = if_else(!is.na(ITS1), ITS1, Length)) %>%
  select(-ITS1)
# Add Sorochytrium
soro.df <- data.frame(NA,NA,"Sorochytrium milnesiophthora",4158)
names(soro.df) <- names(itsonedb.rep)
itsonedb.rep <- rbind(itsonedb.rep, soro.df)

# ROD
# Load data
rod.its1 <- read_tsv("ROD_lengths_ITS1.txt", col_names = TRUE)
rod.5_8S <- read_tsv("ROD_lengths_5_8S.txt", col_names = TRUE)
rod.its2 <- read_tsv("ROD_lengths_ITS2.txt", col_names = TRUE) 
rod.full <- read_tsv("ROD_lengths_full.txt", col_names = TRUE)
rod.infernal <- read_tsv("ROD_lengths_infernal.tsv", col_names = c("File", "ITS1", "Five", "ITS2", "Full", "Notes"), skip=1)
rod.infernal[which(rod.infernal$Notes == "5.8S_dup"),]
#rod.infernal <- rod.infernal[-which(rod.infernal$Notes == "5.8S_dup"),]
rod.infernal <- separate_wider_delim(rod.infernal, cols = File, delim = "-", names = c("Region", "Accession", "Sequence", "Sequence_end")) 
rod.infernal$Sequence <- paste(rod.infernal$Sequence, rod.infernal$Sequence_end, sep = "-")
rod.infernal$Sequence <- gsub("_", "/", rod.infernal$Sequence)
rod.infernal <- rod.infernal[!duplicated(rod.infernal[,-c(1,4)]),][,-c(1,4)]
# Combine data
rod.df <- full_join(rod.its1, rod.its2, by = c("Accession", "Sequence", "Taxonomy")) 
colnames(rod.df) <- c("Accession", "Sequence", "Taxonomy", "ITS1", "ITS2")
rod.df <- full_join(rod.df, rod.full, by = c("Accession", "Sequence", "Taxonomy")) 
colnames(rod.df) <- c("Accession", "Sequence", "Taxonomy", "ITS1", "ITS2", "Full")
rod.df$Full <- as.numeric(rod.df$Full)
rod.df <- full_join(rod.df, rod.5_8S, by = c("Accession", "Sequence", "Taxonomy")) 
colnames(rod.df) <- c("Accession", "Sequence", "Taxonomy", "ITS1", "ITS2", "Full", "Five")
# Replace with Infernal values - I originally only replaced if not NA, but actually NA from Infernal is also informative
# rod.df <- rod.df %>%
#   left_join(rod.infernal %>% select(Accession, Sequence, ITS1, ITS2, Full, Five), by = c("Accession", "Sequence"), suffix = c("", ".new")) %>%
#   mutate(
#     ITS1 = if_else(!is.na(ITS1.new), ITS1.new, ITS1),
#     ITS2 = if_else(!is.na(ITS2.new), ITS2.new, ITS2),
#     Full = if_else(!is.na(Full.new), Full.new, Full),
#     Five = if_else(!is.na(Five.new), Five.new, Five)
#   ) %>%
#   select(-ITS1.new, -ITS2.new, -Full.new, -Five.new)
# Replace with Infernal values - added OR statement so will replace even NAs if also with match
rod.df <- rod.df %>%
  left_join(rod.infernal %>% select(Accession, Sequence, ITS1, ITS2, Full, Five), by = c("Accession", "Sequence"), suffix = c("", ".new")) %>%
  mutate(
    ITS1 = if_else(!is.na(ITS1.new) | (Accession %in% rod.infernal$Accession & Sequence %in% rod.infernal$Sequence), ITS1.new, ITS1),
    ITS2 = if_else(!is.na(ITS2.new) | (Accession %in% rod.infernal$Accession & Sequence %in% rod.infernal$Sequence), ITS2.new, ITS2),
    Full = if_else(!is.na(Full.new) | (Accession %in% rod.infernal$Accession & Sequence %in% rod.infernal$Sequence), Full.new, Full),
    Five = if_else(!is.na(Five.new) | (Accession %in% rod.infernal$Accession & Sequence %in% rod.infernal$Sequence), Five.new, Five)
  ) %>%
  select(-ITS1.new, -ITS2.new, -Full.new, -Five.new)
# Remove bad rows
rod.df.na <- rod.df[which(rowSums(is.na(rod.df[,4:7])) > 2),]
rod.df.bad <- rod.df[!is.na(rod.df$Five) & (rod.df$Five >= 200 | rod.df$Five < 0),] # I don't want to remove the NA values
rod.df <- rod.df[!(rod.df$Sequence %in% rod.df.bad$Sequence),]
rod.df.bad <- rod.df[!is.na(rod.df$ITS2) & (rod.df$ITS2 < 0),] # I don't want to remove the NA values
rod.df <- rod.df[!(rod.df$Sequence %in% rod.df.bad$Sequence),]
accessions.bad <- c("CATKWG010002106.1/226048-210526", "CATKWG010001933.1/26808-11404") # Removed because Infernal was wrong and sequences from same gneome in db
rod.df <- rod.df[!rod.df$Sequence %in% accessions.bad,]
# Separate taxonomy and simplify kingdom classifications
rod.df <- separate_wider_delim(rod.df, cols = Taxonomy, delim = ";", names = c("Domain", "Supergroup", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) 
protists <- unique(rod.df$Kingdom)[c(3:5,7:19)]
rod.df$Kingdom[rod.df$Kingdom %in% protists] <- "Other"
# Add Sorochytrium
soro.df <- data.frame(NA,NA,"Eukaryota","Opisthokonta","Fungi",
             "Blastocladiomycota","Blastocladiomycetes","Blastocladiales",
             "Sorochytriaceae","Sorochytrium","Sorochytrium_milnesiophthora",
             4158,1217,5529,154)
names(soro.df) <- names(rod.df)
rod.df <- rbind(rod.df, soro.df)
# Subset data
rod.df.full <- rod.df[!is.na(rod.df$Full),]
rod.df.full <- rod.df.full[!is.na(rod.df.full$ITS1) & !is.na(rod.df.full$ITS2),] # 6013 rows with NA values in ITS1 or ITS2, which do not plot
# Check data
# summary(rod.infernal)
# length(which(rod.infernal$ITS1 + rod.infernal$Five + rod.infernal$ITS2 - 1 == rod.infernal$Full)) # 2927 rows have ITS subregions that sum to full ITS length, 2 NAs
# summary(rod.df)
# length(which(rod.df$ITS1 + rod.df$Five + rod.df$ITS2 + 1 != rod.df$Full)) # All ITS1 + 5.8S + ITS2 are same as full length (plus 1 for some reason), except Infernal sequences
# length(which(rod.df$ITS1 + rod.df$Five + rod.df$ITS2 - 1 == rod.df$Full)) # The Infernal sequences are minus 1 for some reason
# rod.infernal %>%
#   count(Sequence) %>%
#   filter(n > 1)
# rod.infernal %>%
#   count(Accession, Sequence) %>%
#   filter(n > 1)

# UNITE
# Load data
unite.its1 <- read_tsv("UNITE_lengths_ITS1.txt", col_names = TRUE)
unite.5_8S <- read_tsv("UNITE_lengths_5_8S.txt", col_names = TRUE)
unite.its2 <- read_tsv("UNITE_lengths_ITS2.txt", col_names = TRUE) 
unite.full <- read_tsv("UNITE_lengths_full.txt", col_names = TRUE)
unite.full$Length <- as.numeric(unite.full$Length)
unite.infernal <- read_tsv("UNITE_lengths_infernal.tsv", col_names = c("Accession", "ITS1", "Five", "ITS2", "Full", "Notes"), skip=1)
unite.infernal.short <- read_tsv("UNITE_lengths_infernal_short.tsv", col_names = c("Accession", "ITS1", "Five", "ITS2", "Full", "Notes"), skip=1)
# Combine data
unite.df <- full_join(unite.its1, unite.its2, by = c("Accession", "Taxonomy")) 
colnames(unite.df) <- c("Accession", "Taxonomy", "ITS1", "ITS2")
unite.df <- full_join(unite.df, unite.full, by = c("Accession", "Taxonomy")) 
colnames(unite.df) <- c("Accession", "Taxonomy", "ITS1", "ITS2", "Full")
unite.df <- full_join(unite.df, unite.5_8S, by = c("Accession", "Taxonomy")) 
colnames(unite.df) <- c("Accession", "Taxonomy", "ITS1", "ITS2", "Full", "Five")
# Replace with Infernal values - added OR statement so will replace even NAs if also with match
unite.df <- unite.df %>%
  left_join(unite.infernal %>% select(Accession, ITS1, ITS2, Full, Five), by = c("Accession"), suffix = c("", ".new")) %>%
  mutate(
    ITS1 = if_else(!is.na(ITS1.new) | (Accession %in% unite.infernal$Accession), ITS1.new, ITS1),
    ITS2 = if_else(!is.na(ITS2.new) | (Accession %in% unite.infernal$Accession), ITS2.new, ITS2),
    Full = if_else(!is.na(Full.new) | (Accession %in% unite.infernal$Accession), Full.new, Full),
    Five = if_else(!is.na(Five.new) | (Accession %in% unite.infernal$Accession), Five.new, Five)
  ) %>%
  select(-ITS1.new, -ITS2.new, -Full.new, -Five.new)
unite.df <- unite.df %>%
  left_join(unite.infernal.short %>% select(Accession, ITS1, ITS2, Full, Five), by = c("Accession"), suffix = c("", ".new")) %>%
  mutate(
    ITS1 = if_else(!is.na(ITS1.new) | (Accession %in% unite.infernal.short$Accession), ITS1.new, ITS1),
    ITS2 = if_else(!is.na(ITS2.new) | (Accession %in% unite.infernal.short$Accession), ITS2.new, ITS2),
    Full = if_else(!is.na(Full.new) | (Accession %in% unite.infernal.short$Accession), Full.new, Full),
    Five = if_else(!is.na(Five.new) | (Accession %in% unite.infernal.short$Accession), Five.new, Five)
  ) %>%
  select(-ITS1.new, -ITS2.new, -Full.new, -Five.new)
# Remove bad rows
unite.df.na <- unite.df[which(rowSums(is.na(unite.df[,3:6])) > 2),]
unite.df.bad <- unite.df[!is.na(unite.df$Five) & (unite.df$Five >= 200 | unite.df$Five < 0),] # I don't want to remove the NA values
unite.df <- unite.df[!(unite.df$Accession %in% unite.df.bad$Accession),]
unite.df.bad <- unite.df[!is.na(unite.df$ITS1) & (unite.df$ITS1 < 0),] # I don't want to remove the NA values
unite.df <- unite.df[!(unite.df$Accession %in% unite.df.bad$Accession),]
unite.df.bad <- unite.df[!is.na(unite.df$ITS2) & (unite.df$ITS2 < 0),] # I don't want to remove the NA values
unite.df <- unite.df[!(unite.df$Accession %in% unite.df.bad$Accession),]
accessions.bad <- c("MF148321")
unite.df <- unite.df[!unite.df$Accession %in% accessions.bad,]
# Separate taxonomy and simplify kingdom
unite.df <- separate_wider_delim(unite.df, cols = Taxonomy, delim = ";", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) 
unite.df <- unite.df %>% mutate(Kingdom = str_sub(Kingdom, 4, -1),
         Phylum = str_sub(Phylum, 4, -1),
         Class = str_sub(Class, 4, -1),
         Order = str_sub(Order, 4, -1),
         Family = str_sub(Family, 4, -1),
         Genus = str_sub(Genus, 4, -1),
         Species = str_sub(Species, 4, -1))
protists <- unique(unite.df$Kingdom)[c(2,4:8,10:22)]
unite.df$Kingdom[unite.df$Kingdom %in% protists] <- "Other"
# Add Sorochytrium
soro.df <- data.frame(NA,"Fungi",
                      "Blastocladiomycota","Blastocladiomycetes","Blastocladiales",
                      "Sorochytriaceae","Sorochytrium","Sorochytrium_milnesiophthora",
                      4158,1217,5529,154)
names(soro.df) <- names(unite.df)
unite.df <- rbind(unite.df, soro.df)
# Subset data
unite.df.full <- unite.df[!is.na(unite.df$Full),]
unite.df.full <- unite.df.full[!is.na(unite.df.full$ITS1) & !is.na(unite.df.full$ITS2),] # 4021 rows with NA values in ITS1 or ITS2, which do not plot
# Check data
# which(unite.df$ITS1 + unite.df$ITS2 > unite.df$Full) # No ITS1 + ITS2 are greater than full length, good
# which(unite.df$ITS1 + unite.df$Five + unite.df$ITS2 + 1 != unite.df$Full) # All ITS1 + 5.8S + ITS2 are same as full length (plus 1 for some reason), good

# Combine ROD and UNITE
combined.df <- rbind(rod.df.full[,-c(2:4)], unite.df.full)

################################################################################

# EXPLORE DATA

# ITSoneDB - 2000 long complete
# ENA-identified sequences (617,004 sequences)
sort(itsonedb.ena$Length, decreasing=T)
max(itsonedb.ena$Length, na.rm = T) # 6007 bp
mean(itsonedb.ena$Length, na.rm = T) # 241 bp
median(itsonedb.ena$Length, na.rm = T) # 231 bp
# HMM-identified sequences (548,081 sequences)
sort(itsonedb.hmm$Length, decreasing=T)
max(itsonedb.hmm$Length, na.rm = T) # 15416 bp
mean(itsonedb.hmm$Length, na.rm = T) # 197 bp
median(itsonedb.hmm$Length, na.rm = T) # 180 bp
# Both methods (130,705 sequences)
sort(itsonedb.both$Length, decreasing=T)
max(itsonedb.both$Length, na.rm = T) # 4294 bp
mean(itsonedb.both$Length, na.rm = T) # 239 bp
median(itsonedb.both$Length, na.rm = T) # 222 bp
# Rep sequences, Infernal corrected (186,438 sequences)
sort(itsonedb.rep$Length, decreasing=T)
max(itsonedb.rep$Length, na.rm = T) # 4609 bp
mean(itsonedb.rep$Length, na.rm = T) # 421 bp
median(itsonedb.rep$Length, na.rm = T) # 405 bp
sum(itsonedb.rep$Length > 1250) # 287 sequences greater than 1250 bp
# Infernal - values are larger because I only analyzed large sequences
sort(itsonedb.infernal$ITS1, decreasing=T)
max(itsonedb.infernal$ITS1, na.rm = T) # 4609 bp
mean(itsonedb.infernal$ITS1, na.rm = T) # 943 bp
median(itsonedb.infernal$ITS1, na.rm = T) # 809 bp
sort(itsonedb.infernal$Five, decreasing=T)
sort(itsonedb.infernal$ITS2, decreasing=T)
sort(itsonedb.infernal$Full, decreasing=T) # No full sequences because FASTA analyzed just flanking regions to ITS1
# ITS1 sequences larger than 400 bp
sum(itsonedb.rep$Length >= 400)

# ROD - 2000 long complete
# Full ITS length
sort(rod.df$Full, decreasing=T)
sort(rod.df$Full, decreasing=F)
max(rod.df$Full, na.rm = T) # 1788 bp, 1787 bp
mean(rod.df$Full, na.rm = T) # 672 bp
median(rod.df$Full, na.rm = T) # 610 bp
sd(rod.df$Full, na.rm = T) # 192 bp
rod.df[which(rod.df$Full == 1788),]
sort(rod.df.full$Full, decreasing=F)
sort(rod.df.full$ITS1, decreasing=F)
rod.df.full[which(rod.df.full$ITS1 == 61),]
sort(rod.df.full$ITS2, decreasing=F)
rod.df.full[which(rod.df.full$ITS2 == 76),]
rod.df.full[which(rod.df.full$ITS1 < 100 & rod.df.full$ITS2 > 400),]
# ITS1
sort(rod.df$ITS1, decreasing=T)
sort(rod.df$ITS1, decreasing=F)
max(rod.df$ITS1, na.rm = T) # 849 bp
mean(rod.df$ITS1, na.rm = T) # 272 bp
median(rod.df$ITS1, na.rm = T) # 236 bp
rod.df[which(rod.df$ITS1 == 849),]
sort(rod.its1$Length, decreasing=T)
sort(rod.infernal$ITS1, decreasing=T)
rod.df.full[rod.df.full$Kingdom == "Fungi",] # 9390 fungi 
large <- rod.df.full[rod.df.full$Kingdom == "Fungi" & rod.df.full$ITS1 >= 400,] # 90 fungi over 400
# Portion fo phyla (0.2% over 400 ITS1)
per1 <- rod.df[rod.df$Kingdom == "Fungi",] %>%
  count(Phylum) %>%
  mutate(percentage = (n / nrow(rod.df[rod.df$Kingdom == "Fungi",])) * 100)
per2 <- large  %>%
  count(Phylum)%>%
  mutate(percentage = (n / nrow(large)) * 100)
# ITS2
sort(rod.df$ITS2, decreasing=T)
sort(rod.df$ITS2, decreasing=F)
max(rod.df$ITS2, na.rm = T) # 907 bp
mean(rod.df$ITS2, na.rm = T) # 239 bp
median(rod.df$ITS2, na.rm = T) # 217 bp
rod.df[which(rod.df$ITS2 == 907),] 
sort(rod.its2$Length, decreasing=T)
sort(rod.infernal$ITS2, decreasing=T)
# 5.8S
sort(rod.df$Five, decreasing=T)
max(rod.df$Five, na.rm = T) # 199 bp because remove sequences with longer 5.8S
mean(rod.df$Five, na.rm = T) # 160 bp
median(rod.df$Five, na.rm = T) # 160 bp

# UNITE - 2000 long complete, 2000 short complete
# Full ITS length
sort(unite.df$Full, decreasing = T)
max(unite.df$Full, na.rm = T) # 1883 bp, 1428 bp
mean(unite.df$Full, na.rm = T) # 553 bp
median(unite.df$Full, na.rm = T) # 535 bp
sd(unite.df$Full, na.rm = T) # 106 bp
unite.df[which(unite.df$Full == 1883),] # Succinea_putris, mollusc
unite.df[which(unite.df$Full == 1428),] # Cyathodium_spruceanum, marchantiophyta
sort(unite.df.full$Full, decreasing = F)
unite.df[which(unite.df$Full == 184),] # Rhizocarpon_rubescens, Lecanoromycetes
# ITS1
sort(unite.df$ITS1, decreasing=T)
sort(unite.df$ITS1, decreasing=F)
max(unite.df$ITS1, na.rm = T) # 1890 bp, 1696 bp
mean(unite.df$ITS1, na.rm = T) # 197 bp
median(unite.df$ITS1, na.rm = T) # 187 bp
unite.df[which(unite.df$ITS1 ==1890),]
unite.infernal[which(unite.infernal$Accession=="KC602370"),]
unite.df[which(unite.df$ITS1 ==1696),]
unite.infernal[which(unite.infernal$Accession=="KP400558"),]
unite.df.full[unite.df.full$Kingdom == "Fungi",] # 293187 fungi 
large <- unite.df.full[unite.df.full$Kingdom == "Fungi" & unite.df.full$ITS1 >= 400,]
# Portion fo phyla (0.2% over 400 ITS1)
per1 <- unite.df.full[unite.df.full$Kingdom == "Fungi",] %>%
  count(Phylum) %>%
  mutate(percentage = (n / nrow(unite.df.full[unite.df.full$Kingdom == "Fungi",])) * 100)
per2 <- large  %>%
  count(Phylum)%>%
  mutate(percentage = (n / nrow(large)) * 100)
# ITS2
sort(unite.df$ITS2, decreasing=T)
sort(unite.df$ITS2, decreasing=F)
max(unite.df$ITS2, na.rm = T) # 1431 bp, 1341 bp
mean(unite.df$ITS2, na.rm = T) # 198 bp
median(unite.df$ITS2, na.rm = T) # 193 bp
unite.df[which(unite.df$ITS2==1431),]
unite.infernal[which(unite.infernal$Accession=="LC535111"),]
unite.df[which(unite.df$ITS2==1341),]
unite.infernal[which(unite.infernal$Accession=="OU942837"),]
unite.df[which(unite.df$ITS2==1164),]
unite.infernal[which(unite.infernal$Accession=="JQ619173"),]
# 5.8S
sort(unite.df$Five, decreasing=T)
max(unite.df$Five, na.rm = T) # 199 bp, removed >= 200
mean(unite.df$Five, na.rm = T) # 158 bp
median(unite.df$Five, na.rm = T) # 158 bp

################################################################################

# VISUALIZE DATA
# Understand normal data ellipses: https://ggplot2.tidyverse.org/reference/stat_ellipse.html
# Discussion on meaning of confidence ellipses: https://stats.stackexchange.com/questions/217374/real-meaning-of-confidence-ellipse

# Assign colors
# Default purple: #C77CFF
group.colors <- c(Fungi = "#C77CFF", Metazoa = "#00BFC4", Viridiplantae = "#7CAE00", Other ="#F8766D")

#### ITSoneDB

# Plot A - histogram of rep dataset
itsonedb.hist <- ggplot(itsonedb.rep, aes(x=Length)) +
  geom_vline(xintercept=itsonedb.rep[itsonedb.rep$Species=="Sorochytrium milnesiophthora",]$Length, color="#C77CFF", linewidth=1, linetype="solid", alpha=0.8) + 
    geom_histogram(alpha=0.8, binwidth=20, color = NA, fill = "black") +
  geom_vline(xintercept=mean(itsonedb.rep$Length), color="grey", linewidth=1, linetype="dashed") +  
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +  # Remove x-axis padding
  scale_y_continuous(expand = c(0,0)) +  # Remove y-axis padding
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITSoneDB ITS1 length histogram",x="ITS1 length (bp)", y = "Count") 
itsonedb.hist
ggsave("figures/ITS1_histogram.pdf", width=10, height=5)

# Plot B - histogram of rep dataset, log transformed count (y axis)
# pseudo_log_trans(sigma = 1) adds a small constant (sigma) before applying the log, so log10(count + 1) instead of log10(count).
itsonedb.hist.log <- ggplot(itsonedb.rep, aes(x=Length)) +
  geom_vline(xintercept=itsonedb.rep[itsonedb.rep$Species=="Sorochytrium milnesiophthora",]$Length, color="#C77CFF", linewidth=1, linetype="solid", alpha=0.8) + 
    geom_histogram(alpha=0.8, binwidth=20, color = NA, fill = "black") +
  geom_vline(xintercept=mean(itsonedb.rep$Length), color="grey", linewidth=1, linetype="dashed") + 
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +  # Remove x-axis padding
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = c(0, 1, 10, 100, 1000, 10000, 20000), expand = c(0,0)) +
  guides(y = guide_axis_logticks(negative_small = 1)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITSoneDB ITS1 length histogram",x="ITS1 length (bp)", y = "Count (log10 scale)")
itsonedb.hist.log
ggsave("figures/ITS1_histogram-log.pdf", width=10, height=5)

# Plot C - density plot of rep dataset
itsonedb.dp <- ggplot(itsonedb.rep, aes(x=Length)) +
  geom_vline(xintercept=itsonedb.rep[itsonedb.rep$Species=="Sorochytrium milnesiophthora",]$Length, color="#C77CFF", linewidth=1, linetype="solid", alpha=0.8) + 
  geom_density(fill="black", color="#e9ecef", alpha=0.8) +
  geom_vline(xintercept=mean(itsonedb.rep$Length), color="grey", linewidth=1, linetype="dashed") + 
  theme_bw() +
  scale_x_continuous(expand = c(0,0)) +  # Remove x-axis padding
  scale_y_continuous(expand = c(0,0)) +  # Remove y-axis padding
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITSoneDB ITS1 length density plot", x="ITS1 length (bp)", y = "Density")
itsonedb.dp
ggsave("figures/ITS1_density.pdf", width=10, height=5)

# Plot D - density plot of rep dataset, log10 transformation of counts
itsonedb.rep.top <- itsonedb.rep %>% arrange(desc(Length)) %>% top_n(6)
itsonedb.dp.log <- ggplot(itsonedb.rep, aes(x=Length)) +
  geom_textvline(data = itsonedb.rep.top, aes(xintercept=Length, label=Species), color = "#00BFC4", linewidth = 1) +
  geom_textvline(label = "Sorochytrium milnesiophthora", xintercept=4158, color = "#C77CFF", linewidth = 1) +
  geom_density(fill="#E69F00", color="#E69F00", alpha=0.5) +
  geom_vline(xintercept=mean(itsonedb.rep$Length), color="#E69F00", linewidth=1, linetype="dashed") + 
  theme_bw() +
  scale_x_continuous(expand = c(0,0), limits=c(0,4700)) +  # Remove x-axis padding
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-7), breaks = c(0.000, 0.0001, 0.001), expand = c(0,0)) +
  guides(y = guide_axis_logticks(negative_small = 0.00001)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITSoneDB ITS1 length density plot", x="ITS1 length (bp)", y = "Density (log10 scale)")
itsonedb.dp.log
ggsave("figures/ITS1_density-log.pdf", width=10, height=5)

# Plot E - ridgeline plot showing different datasets
# See https://r-graph-gallery.com/294-basic-ridgeline-plot.html
# Not nice looking, too right skewed
itsonedb.df <- rbind(itsonedb.ena, itsonedb.hmm, itsonedb.both, itsonedb.rep)
itsonedb.rp <- ggplot(itsonedb.df, aes(x = Length, y = File, fill = File)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
itsonedb.rp
ggsave("figures/ITS1_ridgeplot.pdf", width=5, height=5)

# Plot F - ridgeline plot showing different datasets, reduced length
# Better, interesting, but not really showing much
itsonedb.rp.reduced <- ggplot(subset(itsonedb.df, Length < 1000), aes(x = Length, y = File, fill = File)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
itsonedb.rp.reduced
ggsave("figures/ITS1_ridgeplot-reduced.pdf", width=5, height=5)

# Plot G - overlapping density plots for different datasets, reduced length
itsonedb.dp.overlap <- ggplot(subset(itsonedb.df, Length<1000), aes(x=Length, color=File, fill=File)) +
  geom_density(alpha=0.3) +
  geom_vline(xintercept=mean(itsonedb.df[itsonedb.df$File=="rep",]$Length, na.rm=T), color="#999999", linewidth=1, linetype="dashed") + 
  geom_vline(xintercept=mean(itsonedb.df[itsonedb.df$File=="both",]$Length, na.rm=T), color="#E69F00", linewidth=1, linetype="dashed") + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9", "red")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "red")) +
  labs(title="ITSoneDB ITS1 length density plot",x="ITS1 length (bp)", y = "Density") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0)) +  # Remove x-axis padding
  scale_y_continuous(expand = c(0,0)) +  # Remove y-axis padding
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22))
itsonedb.dp.overlap 
ggsave("figures/ITS1_density_overlap.pdf", width=10, height=5)

#### ROD

# Plot A - ROD ITS1 and ITS2 scatter plot for full length ITS sequences
rod.sp <- ggplot(rod.df.full, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  stat_ellipse() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
rod.sp
ggsave("figures/ITS_scatterplot_ROD.pdf", rod.sp, width=6, height=5)

# Plot B - ROD ITS1 and ITS2 scatter plot for full length ITS sequences with sum of ITS1 and ITS2 lines (remember, missing 5.8S length, not full ITS)
k_values <- c(250, 500, 750, 1000, 1250)
line_data <- data.frame(
  intercept = k_values,
  slope = -1 # A line of the form y = k - x has a slope of -1
)
rod.sp <- ggplot(rod.df.full, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_abline(data = line_data, aes(intercept = intercept, slope = slope),
              color = "grey", linetype = "dashed", alpha=0.5) +
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  stat_ellipse() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
rod.sp
ggsave("figures/ITS_scatterplot_ROD-lines.pdf", rod.sp, width=6, height=5)

#### UNITE

# Plot A - UNITE ITS1 and ITS2 scatter plot for full length ITS sequences
unite.sp <- ggplot(unite.df.full, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_abline(data = line_data, aes(intercept = intercept, slope = slope),
              color = "grey", linetype = "dashed", alpha=0.5) +
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  stat_ellipse() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
unite.sp
ggsave("figures/ITS_scatterplot_UNITE-full.pdf", unite.sp, width=6, height=5)

# Plot B - UNITE ITS1 and ITS2 scatter plot for full length ITS sequences, sampled for reduced dataset size
sum(unite.df.full$Full > 1000) # 3684
sum(unite.df.full$Full < 1000) # 422090
sum(unite.df.full$Full < 500) # 159756
unite.df.full.sampled <- unite.df.full[-sample(which(unite.df.full$Full<500), 100000),]
unite.df.full.sampled <- unite.df.full.sampled[-sample(which(unite.df.full.sampled$Full<750), 75000),]
unite.df.full.sampled <- unite.df.full.sampled[-sample(which(unite.df.full.sampled$Full<800), 75000),]
unite.df.full.sampled <- unite.df.full.sampled[-sample(which(unite.df.full.sampled$Full<1000), 75000),]
summary(unite.df.full.sampled)
unite.sp.sampled <- ggplot(unite.df.full.sampled, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_abline(data = line_data, aes(intercept = intercept, slope = slope),
              color = "grey", linetype = "dashed", alpha=0.5) +
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  stat_ellipse() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
unite.sp.sampled
ggsave("figures/ITS_scatterplot_UNITE-sampled.pdf", unite.sp.sampled, width=6, height=5)

#### Combined

# Plot A - combined ROD and UNITE ITS datsets scatterplot
combined.sp <- ggplot(combined.df, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=2,stroke=NA) +
  stat_ellipse() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp
ggsave("figures/ITS_scatterplot_combined-full.pdf", combined.sp, width=6, height=5)
# ggsave("figures/ITS_scatterplot_combined-full-ding.pdf", combined.sp, width=6, height=5, useDingbats = TRUE) # Doesn't reduce file size

# Plot B - combined ROD and UNITE ITS datsets scatterplot, sampled to reduce size
sum(combined.df$Full > 1000)
sum(combined.df$Full < 1000)
sum(combined.df$Full < 500)
combined.df.sampled <- combined.df[-sample(which(combined.df$Full<500), 150000),]
combined.df.sampled <- combined.df.sampled[-sample(which(combined.df.sampled$Full<750), 150000),]
combined.df.sampled <- combined.df.sampled[-sample(which(combined.df.sampled$Full<800), 50000),]
combined.df.sampled <- combined.df.sampled[-sample(which(combined.df.sampled$Full<1000), 50000),]
summary(combined.df.sampled)
combined.sp.sampled <- ggplot(combined.df.sampled, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  stat_ellipse() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.sampled
ggsave("figures/ITS_scatterplot_combined-sampled.pdf", combined.sp.sampled, width=6, height=5)

# Plot C - combined ROD and UNITE hexbin graph
# https://r-charts.com/correlation/hexbin-chart-ggplot2/
combined.hx <- ggplot(combined.df, aes(x=ITS1, y=ITS2)) + 
  geom_hex(binwidth = 20) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="ITS1 and ITS2 length hexplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.hx
ggsave("figures/ITS_hexplot_combined-full.pdf", combined.hx, width=6, height=5)

# Plot D - combined ROD and UNITE smoothed scatterplot
# Information here: https://r-charts.com/correlation/smooth-scatter-plot/
# I could try to replicate with ggplot... maybe... https://stackoverflow.com/questions/13094827/how-to-reproduce-smoothscatters-outlier-plotting-in-ggplot
palette <- hcl.colors(30, palette = "inferno")
combined.smooth <- smoothScatter(combined.df$ITS1, 
                                 combined.df$ITS2, 
                                 nrpoints = 0,
                                 xaxs="i", yaxs="i",
                                 xlab = "ITS1", ylab = "ITS2",)
dev.print(pdf, 'ITS_smooth_combined-full1.pdf')
combined.smooth <- smoothScatter(combined.df$ITS1, 
                                 combined.df$ITS2, 
                                 nrpoints = 0,
                                 bandwidth = 1,
                                 nbin = 1000,
                                 xaxs="i", yaxs="i",
                                 xlab = "ITS1", ylab = "ITS2",)
dev.print(pdf, 'ITS_smooth_combined-full2.pdf')
combined.smooth <- smoothScatter(combined.df$ITS1, 
                                 combined.df$ITS2, 
                                 nrpoints = 0,
                                 xaxs="i", yaxs="i",
                                 xlab = "ITS1", ylab = "ITS2",
                                 colramp = colorRampPalette(palette))
dev.print(pdf, 'ITS_smooth_combined-full3.pdf')
combined.smooth <- smoothScatter(combined.df$ITS1, 
                                 combined.df$ITS2, 
                                 nrpoints = 0,
                                 bandwidth = 1,
                                 nbin = 200,
                                 xaxs="i", yaxs="i",
                                 xlab = "ITS1", ylab = "ITS2",
                                 colramp = colorRampPalette(palette))
dev.print(pdf, 'ITS_smooth_combined-full4.pdf')
combined.smooth <- smoothScatter(combined.df$ITS1, 
                                 combined.df$ITS2, 
                                 nrpoints = 0,
                                 nbin = 1000,
                                 xaxs="i", yaxs="i",
                                 xlab = "ITS1", ylab = "ITS2",
                                 colramp = colorRampPalette(palette))
dev.print(pdf, 'ITS_smooth_combined-full5.pdf')
combined.smooth <- smoothScatter(combined.df$ITS1, 
                                 combined.df$ITS2, 
                                 nrpoints = 0,
                                 bandwidth = 1,
                                 nbin = 1000,
                                 xaxs="i", yaxs="i", 
                                 xlab = "ITS1", ylab = "ITS2",
                                 colramp = colorRampPalette(palette))
dev.print(pdf, 'ITS_smooth_combined-full6.pdf')

# Plot E - combined ROD and UNITE ITS datsets scatterplot, fungi only
combined.sp.fungi <- ggplot(subset(combined.df, Kingdom=="Fungi"), aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="Fungi ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.fungi
ggsave("figures/ITS_scatterplot_combined-fungi.pdf", combined.sp.fungi, width=6, height=5)

# Plot E - combined ROD and UNITE ITS datsets scatterplot, fungi only, phylum labelled
combined.sp.fungi.phylum <- ggplot(subset(combined.df, Kingdom=="Fungi"), aes(x=ITS1, y=ITS2, color=Phylum)) + 
  geom_point(alpha = 1,stroke=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="Fungi ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.fungi.phylum
ggsave("figures/ITS_scatterplot_combined-fungi-phyla.pdf", combined.sp.fungi.phylum, width=6, height=5)

# Plot F - combined ROD and UNITE ITS datsets scatterplot, fungi only, label points
combined.sp.fungi.labelled <- ggplot(subset(combined.df, Kingdom=="Fungi"), aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=2,stroke=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  geom_text_repel(data=subset(combined.df, (ITS1>900 | ITS2>900) & (Kingdom=="Fungi")), aes(ITS1,ITS2,label=Species), max.overlaps=20) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) +
  labs(title="Fungi ITS1 and ITS2 length scatterplot", x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.fungi.labelled
ggsave("figures/ITS_scatterplot_combined-fungi.pdf", combined.sp.fungi, width=6, height=5)

# Plot G - combined ROD and UNITE ITS datasets scatterplot, marginal density plots colored by kingdom
# See https://r-graph-gallery.com/277-marginal-histogram-for-ggplot2.html
combined.sp <- ggplot(combined.df, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=2,stroke=NA) +
  stat_ellipse(
    aes(group = Kingdom),
    type = "norm",
    color="white", alpha=0.5, 
    linewidth = 2
  ) +
  # Foreground colored ellipses for each group
  stat_ellipse(
    aes(group = Kingdom),
    type = "norm",
    linewidth = 0.8,
    linetype = 1
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) + 
  theme(legend.position = "none") +
  labs(x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.dp <- ggMarginal(combined.sp, type="density", groupColour=T, groupFill=T, alpha=0.5,size = 3) 
combined.sp.dp
ggsave("figures/ITS_scatterplot-density_combined.pdf", combined.sp.dp, width=5, height=5)

# Plot H - combined ROD and UNITE ITS datasets scatterplot, marginal density plots colored by kingdom, no ellipses
combined.sp <- ggplot(subset(combined.df, Full<5000), aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_color_manual(values=group.colors) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) + 
  theme(legend.position = "none") +
  labs(x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.dp <- ggMarginal(combined.sp, type="density", groupColour=T, groupFill=T, alpha=0.5,size = 3) 
combined.sp.dp
ggsave("figures/ITS_scatterplot-density_combined_noellipse.pdf", combined.sp.dp, width=5, height=5)

# Plot I - combined ROD and UNITE ITS datasets scatterplot, marginal density plots colored by kingdom, no ellipses
# The marginal plots are off, because of limits expanding the plot, the density plots are not lined up
combined.sp <- ggplot(subset(combined.df, Full<5000), aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  scale_color_manual(values=group.colors) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) + 
  theme(legend.position = "none") +
  labs(x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.dp <- ggMarginal(combined.sp, type="density", groupColour=T, groupFill=T, alpha=0.5,size = 3) 
combined.sp.dp
ggsave("figures/ITS_scatterplot-density_combined_noellipse_padded.pdf", combined.sp.dp, width=5, height=5)

# Plot J - full dataset
combined.sp <- ggplot(combined.df, aes(x=ITS1, y=ITS2, color=Kingdom)) + 
  geom_point(alpha = 0.2,size=1,stroke=NA) +
  scale_color_manual(values=group.colors) +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(hjust=0.5, size=22)) + 
  theme(legend.position = "none") +
  labs(x="ITS1 length (bp)", y = "ITS2 length (bp)")
combined.sp.dp <- ggMarginal(combined.sp, type="density", groupColour=T, groupFill=T, alpha=0.5,size = 3) 
combined.sp.dp
ggsave("figures/ITS_scatterplot-density_combined_noellipse_padded_full.pdf", combined.sp.dp, width=5, height=5)

# Plot K - try removing padding, chatGPT solution using patchwork
# Refactor to show fungi density plot on top
library(patchwork)
combined.df$Kingdom <- fct_relevel(combined.df$Kingdom, "Other", "Viridiplantae", "Metazoa","Fungi") 
p_main <- ggplot(combined.df, aes(x = ITS1, y = ITS2, color = Kingdom)) +
  geom_point(alpha = 0.2, size = 1, stroke = NA) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250), breaks = c(300, 600, 900, 1200)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250), breaks = c(300, 600, 900, 1200)) +
  scale_color_manual(values = group.colors) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 22),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(x = "ITS1 length (bp)", y = "ITS2 length (bp)")
# Top marginal density
p_top <- ggplot(combined.df, aes(x = ITS1, fill = Kingdom, color = Kingdom)) +
  geom_density(alpha = 0.5, size = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_fill_manual(values = group.colors) +
  scale_color_manual(values = group.colors) +
  theme_void() +
  theme(legend.position = "none")
# Right marginal density
p_right <- ggplot(combined.df, aes(x = ITS2, fill = Kingdom, color = Kingdom)) +
  geom_density(alpha = 0.5, size = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1250)) +
  scale_fill_manual(values = group.colors) +
  scale_color_manual(values = group.colors) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none")
# Layout with patchwork
pathplot <- p_top + plot_spacer() + p_main + p_right + plot_layout(ncol=2, nrow=2, heights= c(1,4), widths = c(4, 1))
pathplot
ggsave("figures/ITS_scatterplot-density_combined_noellipse_full_patchwork.pdf", pathplot, width=7.5, height=5)

# Plot L - add outlier panel
# For information on spacing: https://patchwork.data-imaginist.com/articles/guides/layout.html
p_outlier <- ggplot(subset(combined.df, ITS1>3000), aes(x = ITS1, y = ITS2, color = Kingdom)) +
  geom_point(size = 1.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(4100, 4225), breaks = c(4150)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1250), breaks = c(300, 600, 900, 1200)) +
  scale_color_manual(values = group.colors) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),   # Remove redundant y-axis labels
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = "none"
  ) +
  labs(x = "ITS1 length (bp)")
layout <- "
AAAAAAAAAAA###
BBBBBBBBBBCDDD
"
pathplot <- p_top + p_main + p_outlier + p_right +
  plot_layout(design = layout, heights= c(1,3.5))
pathplot
ggsave("figures/ITS_scatterplot-density_combined_noellipse_full_patchwork_outlier.pdf", pathplot, width=6, height=5)
