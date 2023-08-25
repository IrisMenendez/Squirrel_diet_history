
# INFO ----
# Project: GM p4 squirrels - dietary evolution
# R version 4.3.0 (2023-04-21) -- "Already Tomorrow"

# PACKAGES ----
library(StereoMorph)
library(geomorph)
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(dplyr)
library(scales)
library(strap)
library(Momocs)
library(MASS)
library(deeptime)
library(ggplot2)
library(cowplot)

# SEED ----
set.seed(15)

# DATA ----

# + Diet ----
# Read ecological and taxonomical database
squi_table <- read.csv2("data/DATABASE_Ardillas.csv")

# extract desired variables in a dataframe
trait_dat <- data.frame(Sp_red = squi_table$SP_reducido, 
                        Sp = squi_table$HandbookSp, 
                        Genus = squi_table$GENUS.HANDBOOK, 
                        Tribe = squi_table$TRIBE, 
                        Subfamily = squi_table$SUBFAMILY, 
                        Diet = squi_table$diet_Iris, 
                        DietCVA = squi_table$diet_CVA, 
                        meanWeight = squi_table$meanWeight_grams, 
                        stringsAsFactors = T)
# Here we have the info for extant species

# + Morphology ----
# ++ 1. Prepare filenames
# We have morphological data for extant and fossil specimens
# all filenames: 
files_all <- list.files("data/All_shapes")
files_all <- files_all[grep('.txt', files_all)]

# separate filenames in extant and fossil:
files_foss <- files_all[substr(files_all,1,4)=="FOSS"]
files_extant <- files_all[substr(files_all,1,4)!="FOSS"]

# get species names (reduced) (spred) from each extant specimen filename: 
files_extant_spred <-substr(files_extant,6,12)

# dataframe with diet data for all extant sepcimens:
trait_extant <- trait_dat[match(files_extant_spred, trait_dat$Sp_red), ]

# Get the filenames of the specimens belonging to species for which there's trait info:
files_extant_INF <- files_extant[trait_extant$Diet != "NO_INF"]

# Get a vector of the filenames that we will actually read: 
to_read_files <- c(files_foss, files_extant_INF)

# + Phylogeny ----
# Read tree (tip labels named as in the shape data already)
tree <- read.tree("data/phy7_reduced_name.tree")
tree <- ladderize(tree, right = FALSE)
# We have a phylogenetic tree of extant species

# Extract terrestrial (Xerinae) and arboreal (the rest) squirrels in different phylogenies
terr_sp <- as.vector(trait_extant_inf[trait_extant_inf$Subfamily == "Xerinae",]$Sp_red)
arb_sp <-  as.vector(trait_extant_inf[trait_extant_inf$Subfamily != "Xerinae",]$Sp_red)

# drop terrestrial species to have a phylogeny of only arboreal species:
tree_arb <- drop.tip(tree, terr_sp)

# drop arboreal species to have a phylogeny of only terrestrial species:
tree_terr <- drop.tip(tree, arb_sp)

# ++ 2. Import morpho data
# Opening shapes data with StereoMorph (fossil and extant)
Shapes_folder <- "data/All_shapes/"
dodo <- readShapes(paste0(Shapes_folder, to_read_files))
head(dodo$curves.pixel)

# ++ Table of specimens included (TABLE S1) ----
# Create a table with specimen information: 
# - institution
# - specimen name 
# - species reduced name 
# - species full name

# Get names of specimen shapes
shape_names <- row.names(as.matrix(dodo$curves.pixel))

# + 1. species reduced names
# Get species names (reduced) from each specimen shape
shape_spred <-substr(shape_names, 6, 12)

# a) fossil specimens
shape_names_fossil <- shape_names[substr(shape_names, 1, 4) == "FOSS"]
split1 <- unlist(strsplit(shape_names_fossil, "[(]"))
shape_spred_fossil <- substr(split1[seq(1, length(split1), 2)], 6, 12)  

# b) extant specimens
shape_names_extant <- shape_names[substr(shape_names, 1, 4) != "FOSS"]
shape_spred_extant <- substr(shape_names_extant, 6, 12)  

shape_spred_extant %in% shape_spred_fossil #There are species both in extant and fossil dataset!

# + 2. species full names
# Correspondence sp_red and species names (extant):
correspondence_extant <- trait_dat[, c('Sp_red', 'Sp')]
correspondence_extant$Sp_red <- as.character(correspondence_extant$Sp_red)
correspondence_extant$Sp <- as.character(correspondence_extant$Sp)

# Correspondence sp_red and species names (fossil):
ages_fossil <- read.table('data/extinct_sp_names_v2.csv', 
                                    sep = ';', header = TRUE, 
                                    stringsAsFactors = FALSE)

correspondence_fossil <- ages_fossil[, c('sp_short', 'sp_full')]

# + 3. Table extant specimens
table_specimens_extant <- data.frame(Museum = substr(shape_names_extant, 1, 4), 
                        Specimen = regmatches(shape_names_extant, 
                                              regexpr('[0-9]+', 
                                                      shape_names_extant)), 
                        Species_red = substr(shape_names_extant, 6, 12), 
                        Species = NA,
                        stringsAsFactors = FALSE)

table_specimens_extant$Species <- correspondence_extant$Sp[match(table_specimens_extant$Species_red, correspondence_extant$Sp_red)]

# + 4. Table fossil specimens
table_specimens_fossil <- data.frame(Museum = substr(shape_names_fossil, 1, 4), 
                      Specimen = regmatches(shape_names_fossil, 
                                            regexpr('[0-9]+', 
                                                    shape_names_fossil)), 
                      Species_red = substr(shape_names_fossil, 6, 12), 
                      Species = NA, stringsAsFactors = FALSE)
table_specimens_fossil$Species <- correspondence_fossil$sp_full[match(table_specimens_fossil$Species_red, correspondence_fossil$sp_short)]

# + 5. Save table
table_specimens <- rbind(table_specimens_extant, table_specimens_fossil)
write.csv2(table_specimens, "Tables/specimens_table_v4.csv", row.names = FALSE)


# GM DATASET ----
# + Procrustes superimposition ----
p4 <- readland.shapes(dodo,  nCurvePts = c(17, 17, 17, 17), scaled = T)
GPA <- gpagen(p4, curves = p4$curves , ProcD = FALSE, PrinAxes = FALSE) # GPA perfomed with minimized bending energy
plotAllSpecimens(GPA$coords)

# extract coordinates and centroid size (Csize) for fossils
coords_fossil <- GPA$coords[,,shape_names_fossil]
plotAllSpecimens(coords_fossil)
#dimnames(coords_fossil)[[3]] <- shape_spred_fossil

Csize_fossil <- GPA$Csize[shape_names_fossil]
#names(Csize_fossil) <- shape_spred_fossil

# extract coordinates and centroid size (Csize) for extant species
coords_extant <- GPA$coords[,,shape_names_extant]
plotAllSpecimens(coords_extant)
dimnames(coords_extant)[[3]] <- shape_spred_extant

Csize_extant <- GPA$Csize[shape_names_extant]
names(Csize_extant) <- shape_spred_extant

## Get mean centroid size per species (extant species)
group <- as.factor(rownames(two.d.array(coords_extant)))
mCsize_extant <- rowsum(Csize_extant, group)/as.vector(table(group))# rowsum()
max(table(group))
mean(table(group))

# + GM + size dataset ----
# GM + size fossil dataset 
coords_size_fossil <- cbind(two.d.array(coords_fossil), 
                            size = log(Csize_fossil))

# GM + size extant dataset 
coords_size_extant <- cbind(two.d.array(coords_extant), 
                            size = log(Csize_extant))


# EFA DATASET ----
# + Elliptic Fourier Analysis ----
# Transform GM dataset (landmark coordinates) into Fourier coefficients
n <- dim(GPA$coords)[3]
A <- array(NA, dim=c(64,2,n), dimnames = dimnames(GPA$coords)) 
dimnames(GPA$coords)[[3]]
for (i in 1:n){
  A[,,i] <- matrix(as.numeric(c(GPA$coords[,,i][1], 
                                GPA$coords[,,i][5:19],GPA$coords[,,i][2], 
                                GPA$coords[,,i][20:34],GPA$coords[,,i][3], 
                                GPA$coords[,,i][35:49], GPA$coords[,,i][4], 
                                GPA$coords[,,i][50:64], GPA$coords[,,i][65], 
                                GPA$coords[,,i][69:83], GPA$coords[,,i][66], 
                                GPA$coords[,,i][84:98], GPA$coords[,,i][67], 
                                GPA$coords[,,i][99:113], GPA$coords[,,i][68], 
                                GPA$coords[,,i][114:128])), 
                   64,2, byrow=F)
}

mode(A) <- "numeric"


# fourier fossils
b <- length(shape_names_fossil)
Out_foss <- Momocs:::Out(A[,,1:b])
names(Out_foss)
panel(Out_foss)
tf_foss  <- Momocs:::efourier(Out_foss, 7, norm=F, start=F)
#names(tf_foss) <- shape_spred_fossil

# fourier extant
Out_extant <- Momocs::Out(A[,,(b+1):n]) 
panel(Out_extant)
tf_extant  <- Momocs:::efourier(Out_extant, 7, norm=F, start=F)
names(tf_extant) <- shape_spred_extant

# + EFA + size dataset ----
# fourier + size fossils
tf_size_foss <- cbind(tf_foss$coe, Csize = log(Csize_fossil))
rownames(tf_size_foss) <- names(Csize_fossil)

# fourier + size extant
tf_size_extant <- cbind(tf_extant$coe, Csize = log(Csize_extant))
rownames(tf_size_extant) <- names(Csize_extant)


# MATCH DATASETS ----

# data for all extant specimens with dietary information
trait_extant_inf <- droplevels(trait_dat[match(unique(shape_spred_extant),
                                               trait_dat$Sp_red), ])
# functional diet (fundiet)
fundiet_extant <- trait_extant_inf$DietCVA
names(fundiet_extant) <- trait_extant_inf$Sp_red

# + All phylogeny ----
# Match phylogeny and diet information
# Extant species with dietary AND phylogenetic information (.phy.inf)
matchdata_phy <- treedata(tree, mCsize_extant, sort = TRUE)
phy_inf <- matchdata_phy$phy
mCsize_extant_phy <- matchdata_phy$data

# Dataframe with diet info for species in the phylogeny
trait_extant_inf_phy <- trait_extant_inf[match(rownames(mCsize_extant_phy), 
                                          trait_extant_inf$Sp_red),]

# + Arboreal squirrels ----
matchdata_phy_arb <- treedata(tree_arb, mCsize_extant, sort = TRUE)
tree_arb <- matchdata_phy_arb$phy
mCsize_arb <- matchdata_phy_arb$data
trait_arb <- trait_extant_inf[match(rownames(mCsize_arb), trait_extant_inf$Sp_red),]

# + Terrestrial species ----
matchdata_phy_terr <- treedata(tree_terr, mCsize_extant, sort = TRUE)
tree_terr <- matchdata_phy_terr$phy
mCsize_terr <- matchdata_phy_terr$data
trait_terr <- trait_extant_inf[match(rownames(mCsize_terr), trait_extant_inf$Sp_red),]


# COLORS ----
col.diet <- setNames(
  c( "#f4c44e","#6340c6", "#69ccc7", "#7dd86f","#ff595e","#bb3e7b","#4a7f5c","Pink"),
  c("mixed_feeder", "Seeds_nuts", "Browser", "grazer", "Frugivorous", "Bark_gleaner", "Folivore", "Insectivore")
)

col.fundiet <- setNames(c( "#a04a4a", "#add34f", "#006cb5", "Orange"),
                            c("chew_crush", "Grind", "Insectivore", "mixed_feeder"))


# ANALYSES ----
# + PCA ----
# ++ 1. PCA GM dataset ----
PCA_extant_GM <- prcomp(two.d.array(coords_extant))
summary(PCA_extant_GM) #PC6 explains at least 0.90 cumulative

### plot with outlines: https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html

# Convert array to matrix for PCA
gpa_mat <- t(apply(coords_extant, 3, function(y) matrix(t(y), 1)))


# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores
scores <- gpa_mat %*% resEig$vectors
#scores <- scores*-1
# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

sum(per_var)

#define order of landmarks and semilandmarks to draw outlines
order <-c(1,5:19,2,20:34,3,35:49,4,50:64)


# Define function to draw shape
plot_squi_teeth <- function(xy, coor, size=1, col='black'){
  
  # If 3D, rotate points about x-axis using 3D rotation matrix
  if(ncol(coor) == 3){
    coor <- coor %*% matrix(c(1,0,0, 0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), nrow=3, ncol=3)
  }
  
  # Get just x,y coordinates (orthographic projection into xy-plane)
  coor <- coor[, 1:2]
  
  # Get plot aspect ratio
  w <- par('pin')[1]/diff(par('usr')[1:2])
  h <- par('pin')[2]/diff(par('usr')[3:4])
  asp <- w/h
  
  # Correct for plot aspect ratio not necessarily being 1:1
  coor[, 1] <- coor[, 1] * (1/asp)
  
  # Scale points and place back in position
  coor <- coor*size
  
  # Center about zero based on range of coordinates
  coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
                        nrow=nrow(coor), ncol=ncol(coor), byrow=TRUE)
  
  # Move shape to PC score
  coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow=TRUE)
  
  # Set order in which to draw points to create polygon
  polygon_order <- order
  
  # Create filled polygon
  polygon(coor[polygon_order, ], col=col, border=col)
}

# Set PCs to plot
pcs <- c(1,2)
#pcs <- c(1,3)

# Open PDF graphics device
pdf('plots/pca_GM_v4.pdf', width = 9, height = 6.5)

# Create plot box with axes and axis labels
plot(scores[, pcs], type = 'n', main = 'Backtransform morphospace',
     xlab=paste0('PC', pcs[1], ' (', round(per_var[pcs[1]]), '%)'),
     ylab=paste0('PC', pcs[2], ' (', round(per_var[pcs[2]]), '%)'), asp = 1)

# Plot backtransform shapes
btShapes(scores = scores, vectors = resEig$vectors, fcn = plot_squi_teeth, 
         pcs = pcs, n = c(5,5), m = dim(coords_extant)[2], 
         row.names = dimnames(coords_extant)[[1]], 
         pc.margin = c(0.06,0.05), size = 0.1, 
         col = col_alpha("#000000", 0.95))

# Plot points for each species
points(scores[, pcs], 
       bg = alpha(col.fundiet[fundiet_extant[shape_spred_extant]], 0.8), 
       cex = Csize_extant*0.10, lwd = 0.5, pch = 21, col = "transparent")
# Add text labels
#text(scores[, pcs], labels=substr(rownames(scores), 0, 3), cex=0.8, pos=1, offset=0.3)

# Close the PDF graphics device
dev.off()

# ++ 2. PCA GM + size dataset ----
PCA_extant_GM_size <- prcomp(coords_size_extant)
summary(PCA_extant_GM_size)

# ++ 3. PCA EFA dataset ----
PCA_extant_EFA <- prcomp(tf_extant$coe)
rownames(tf_extant$coe)
summary(PCA_extant_EFA)

pdf("plots/pca_fourier_v4.pdf")
plot(PCA(tf_extant), xax = 1, yax = 2, cex = Csize_extant*0.12, 
     nr.shp = 5, nc.shp = 5,
     col = alpha(col.fundiet[fundiet_extant[shape_spred_extant]], 0.8), 
     size.shp = 0.9)
dev.off()

pdf("plots/pca_fourier2_v4.pdf")
plot(PCA(tf_extant), xax = 1, yax = 3, cex = Csize_extant*0.12, 
     nr.shp = 5, nc.shp = 5,
     col = alpha(col.fundiet[fundiet_extant[shape_spred_extant]], 0.8), 
     size.shp = 0.9)
dev.off()



# ++ 4. PCA EFA + size dataset ----
tf_size_extant
PCA_extant_EFA_size <- prcomp(tf_size_extant)
summary(PCA_extant_EFA_size)

# ++ 5. Summary of all PCA ----

summary_pca_gm <- rbind(dataset = 'GM', summary(PCA_extant_GM)$importance)
summary_pca_gm_size <- rbind(dataset = 'GM + size', summary(PCA_extant_GM_size)$importance)
summary_pca_efa <- rbind(dataset = 'EFA', summary(PCA_extant_EFA)$importance)
summary_pca_efa_size <- rbind(dataset = 'EFA + size', summary(PCA_extant_EFA_size)$importance)

write.table(summary_pca_gm, 'Tables/summary_PCA.csv', 
            dec = '.', sep = ';', col.names = NA)
write.table(summary_pca_gm_size, 'Tables/summary_PCA.csv', 
            dec = '.', sep = ';', col.names = NA, append = TRUE)
write.table(summary_pca_efa, 'Tables/summary_PCA.csv', 
            dec = '.', sep = ';', col.names = NA, append = TRUE)
write.table(summary_pca_efa_size, 'Tables/summary_PCA.csv', 
            dec = '.', sep = ';', col.names = NA, append = TRUE)


# + LDA ----
# ++ 1. LDA GM dataset ----
predicted_fossil_PCA_GM <- predict(PCA_extant_GM, 
                                   two.d.array(coords_fossil))
summary(PCA_extant_GM)

LDA_GM <- MASS:::lda(PCA_extant_GM$x[,1:6], 
                     fundiet_extant[shape_spred_extant])
LDA_GM

#plot(lda_landmarksPCA, cex = 0.7, dimen=2, abbrev = TRUE, xlab = "LD1", ylab = "LD2", col=col.dietCVA.inf[data_inf$DietCVA], main="Landmarks")

predicted_fossil_LDA_GM <- MASS:::predict.lda(LDA_GM, 
                                              predicted_fossil_PCA_GM[,1:6])
predicted_fossil_LDA_GM$x
predicted_fossil_LDA_GM$class


new_class_GM <- data.frame(sp = rownames(predicted_fossil_LDA_GM$x), 
                           class = predicted_fossil_LDA_GM$class, 
                           posterior = 1)
for (i in 1:length(predicted_fossil_LDA_GM$class)) {
  new_class_GM[i,3] <- max(predicted_fossil_LDA_GM$posterior[i,])
}


# ++ 2. LDA GM + size dataset ----
predicted_fossil_PCA_GM_size <- predict(PCA_extant_GM_size, 
                                        coords_size_fossil)
summary(PCA_extant_GM_size) #PC9 explains 0.95 cumulative


LDA_GM_size <- MASS:::lda(PCA_extant_GM_size$x[,1:6], 
                          fundiet_extant[shape_spred_extant])

predicted_fossil_LDA_GM_size <- MASS:::predict.lda(LDA_GM_size, 
                                                   predicted_fossil_PCA_GM_size[,1:6] )

#plot(lda_landmarksPCA_size, cex = 0.7, dimen=2, abbrev = TRUE, xlab = "LD1", ylab = "LD2", col=col.dietCVA#.inf[diet.inf[taxa.names.extant]], main="Landmark size")
#plot(lda_landmark_predictPCA_size$x, cex = 0.7, dimen=2, abbrev = TRUE, xlab = "LD1", ylab = "LD2", main#="Landmark size")


new_class_GM_size <- data.frame(sp = rownames(predicted_fossil_LDA_GM_size$x), 
                                class = predicted_fossil_LDA_GM_size$class, 
                                posterior = 1)

for (i in 1:length(predicted_fossil_LDA_GM_size$class)) {
  new_class_GM_size[i,3] <- max(predicted_fossil_LDA_GM_size$posterior[i,])
}





# ++ 3. LDA EFA dataset ----
predicted_fossil_PCA_EFA <- predict(PCA_extant_EFA, 
                                    tf_foss$coe)
summary(predicted_fossil_PCA_EFA) 

LDA_EFA <- lda(PCA_extant_EFA$x[, 1:6], fundiet_extant[shape_spred_extant])
predicted_fossil_LDA_EFA <- MASS:::predict.lda(LDA_EFA, 
                                               predicted_fossil_PCA_EFA[,1:6])

new_class_EFA <- data.frame(sp = shape_names_fossil, 
                            class = predicted_fossil_LDA_EFA$class, 
                            posterior = 1)
for (i in 1:length(predicted_fossil_LDA_EFA$class)) {
  new_class_EFA[i,3] <- max(predicted_fossil_LDA_EFA$posterior[i,])
}



# ++ 4. LDA EFA + size dataset ----

predicted_fossil_PCA_EFA_size <- predict(PCA_extant_EFA_size, tf_size_foss)
summary(predicted_fossil_PCA_EFA_size) 

LDA_EFA_size <- lda(PCA_extant_EFA_size$x[, 1:6], fundiet_extant[shape_spred_extant])

predicted_fossil_LDA_EFA_size <- MASS:::predict.lda(LDA_EFA_size, 
                                               predicted_fossil_PCA_EFA_size[,1:6])

new_class_EFA_size <- data.frame(sp = shape_names_fossil, 
                            class = predicted_fossil_LDA_EFA_size$class, 
                            posterior = 1)
for (i in 1:length(predicted_fossil_LDA_EFA_size$class)) {
  new_class_EFA_size[i,3] <- max(predicted_fossil_LDA_EFA_size$posterior[i,])
}




# ++ 5. Summary of all LDA # ----
SUMMARY_LDA_class <- cbind(new_class_GM,
                           new_class_GM_size[, 2:3],
                           new_class_EFA[, 2:3],
                           new_class_EFA_size[, 2:3])

colnames(SUMMARY_LDA_class)[1] <- 'specimen'
colnames(SUMMARY_LDA_class)[c(4, 6, 8, 10)] <- c('class_GM', 'class_GM_size', 
                                                 'class_EFA', 'class_EFA_size')
SUMMARY_LDA_class$specimen <- as.character(SUMMARY_LDA_class$specimen)
SUMMARY_LDA_class$sp <- substr(SUMMARY_LDA_class$specimen,6,12)



extinct_names <- read.table('data/extinct_sp_names_v2.csv', sep = ';', 
                            header = TRUE, stringsAsFactors = FALSE)
colnames(extinct_names)
extinct_names$specimen_full <- NA

for (i in 1:nrow(extinct_names)){
  extinct_names$specimen_full[i] <- SUMMARY_LDA_class$specimen[grep(extinct_names$specimen[i], SUMMARY_LDA_class$specimen)]
}


SUMMARY_LDA_class$full <- NA

for (i in 1:nrow(SUMMARY_LDA_class)){
  SUMMARY_LDA_class$full[i] <- extinct_names$sp_full[extinct_names$sp_short == SUMMARY_LDA_class$sp[i]][1]
}

# reorder columns
SUMMARY_LDA_class <- SUMMARY_LDA_class %>% 
  relocate(full, .after = specimen) %>% 
  relocate(sp, .after = full)

head(SUMMARY_LDA_class)

# ++ supplementary table (S3) ----
write.table(SUMMARY_LDA_class, "Tables/Summary_cva_v4.csv", sep = ";", dec = ".",
            col.names = TRUE, row.names = FALSE)

# ++ Cross-validation ----
LDA_EFA_CV <- lda(PCA_extant_EFA$x[,1:6], 
                  fundiet_extant[shape_spred_extant], CV = T)
CV_EFA <- table(fundiet_extant[shape_spred_extant], LDA_EFA_CV$class)
CV_EFA_perc <- rbind(CV_EFA[1, ]/sum(CV_EFA[1, ]), 
                     CV_EFA[2, ]/sum(CV_EFA[2, ]),
                     CV_EFA[3, ]/sum(CV_EFA[3, ]),
                     CV_EFA[4, ]/sum(CV_EFA[4, ]))
CV_EFA_perc_total <- sum(diag(CV_EFA))/sum(CV_EFA)

LDA_GM_CV  <- MASS:::lda(PCA_extant_GM$x[,1:6], 
                         fundiet_extant[shape_spred_extant], CV = T)
CV_GM <- table(fundiet_extant[shape_spred_extant],  LDA_GM_CV$class)
CV_GM_perc <- rbind(CV_GM[1, ]/sum(CV_GM[1, ]), 
                    CV_GM[2, ]/sum(CV_GM[2, ]),
                    CV_GM[3, ]/sum(CV_GM[3, ]),
                    CV_GM[4, ]/sum(CV_GM[4, ]))
CV_GM_perc_total <- sum(diag(CV_GM))/sum(CV_GM)

LDA_EFA_size_CV <- lda(PCA_extant_EFA_size$x[,1:6], 
                  fundiet_extant[shape_spred_extant], CV = T)
CV_EFA_size <- table(fundiet_extant[shape_spred_extant], LDA_EFA_size_CV$class)
CV_EFA_size_perc <- rbind(CV_EFA_size[1, ]/sum(CV_EFA_size[1, ]), 
                     CV_EFA_size[2, ]/sum(CV_EFA_size[2, ]),
                     CV_EFA_size[3, ]/sum(CV_EFA_size[3, ]),
                     CV_EFA_size[4, ]/sum(CV_EFA_size[4, ]))
CV_EFA_size_perc_total <- sum(diag(CV_EFA_size))/sum(CV_EFA_size)

LDA_GM_size_CV  <- MASS:::lda(PCA_extant_GM_size$x[,1:6], 
                         fundiet_extant[shape_spred_extant], CV = T)
CV_GM_size <- table(fundiet_extant[shape_spred_extant],  LDA_GM_size_CV$class)
CV_GM_size_perc <- rbind(CV_GM_size[1, ]/sum(CV_GM_size[1, ]), 
                    CV_GM_size[2, ]/sum(CV_GM_size[2, ]),
                    CV_GM_size[3, ]/sum(CV_GM_size[3, ]),
                    CV_GM_size[4, ]/sum(CV_GM_size[4, ]))
CV_GM_size_perc_total <- sum(diag(CV_GM_size))/sum(CV_GM_size)

c(CV_EFA_perc_total, CV_GM_perc_total, 
CV_EFA_size_perc_total, CV_GM_size_perc_total)


CV_GM_size_missclasifications <- data.frame(true = fundiet_extant[shape_spred_extant], 
                                            inferred = LDA_GM_size_CV$class, 
                                            posterior = LDA_GM_size_CV$posterior)
write.table(CV_GM_size_missclasifications, "Tables/CV_GM_size_missclasifications_v4.csv", dec=".", sep=";", col.names=NA)


write.table(CV_GM_size, "Tables/CV_landmarks_size_v4.csv", dec=".", sep=";", col.names=NA)
write.table(CV_GM_size_perc, "Tables/CV_landmarks_size_v4.csv", dec=".", sep=";", col.names=NA, append=T)
write.table(CV_GM_size_perc_total, "Tables/CV_landmarks_size_v4.csv", dec=".", sep=";", col.names=NA, append=T)


SUMMARY_LDA_CV <- rbind(EFA = CV_EFA, GM = CV_GM, EFA_size = CV_EFA_size, GM_size = CV_GM_size)
SUMMARY_LDA_CV_perc <- rbind(EFA = CV_EFA_perc, GM = CV_GM_perc, EFA_size = CV_EFA_size_perc, GM_size = CV_GM_size_perc)
SUMMARY_LDA_percTOTAL <- rbind(EFA = CV_EFA_perc_total, GM = CV_GM_perc_total, EFA_size = CV_EFA_size_perc_total, GM_size = CV_GM_size_perc_total)


summary_lda_GM <- cbind(rbind(CV_GM, GM = NA), 
      rbind(CV_GM_perc*100, GM = NA), 
      TOTAL = c(NA, NA, NA, NA, CV_GM_perc_total*100))

summary_lda_GM_size <- cbind(rbind(CV_GM_size, GM_size = NA), 
                        rbind(CV_GM_size_perc*100, GM_size = NA), 
                        TOTAL = c(NA, NA, NA, NA, CV_GM_size_perc_total*100))

summary_lda_EFA <- cbind(rbind(CV_EFA, EFA = NA), 
                        rbind(CV_EFA_perc*100, EFA = NA), 
                        TOTAL = c(NA, NA, NA, NA, CV_EFA_perc_total*100))

summary_lda_EFA_size <- cbind(rbind(CV_EFA_size, EFA_size = NA), 
                        rbind(CV_EFA_size_perc*100, EFA_size = NA), 
                        TOTAL = c(NA, NA, NA, NA, CV_EFA_size_perc_total*100))

summary_lda <- rbind(summary_lda_GM, summary_lda_GM_size, 
      summary_lda_EFA, summary_lda_EFA_size)
colnames(summary_lda)[5:8] <- paste0('% ', colnames(summary_lda)[5:8])

write.table(summary_lda, 'Tables/summary_LDA.csv', 
            dec = '.', sep = ';', col.names = NA)

t(summary_lda)


# + ANCESTRAL RECONSTRUCTION ----


# Diet as factors 
diet_phy <- trait_extant_inf_phy$Diet
names(diet_phy) <- trait_extant_inf_phy$Sp_red

diet_phy_arb <- trait_arb$Diet
names(diet_phy_arb) <- trait_arb$Sp_red

diet_phy_terr <- trait_terr$Diet
names(diet_phy_terr) <- trait_terr$Sp_red

# ++ Entire phylogeny ----

# Fit models for transition rates
# ER is "equal rates", "Sym" is "symmetric rates" and "ARD" is "all rates differ"
#geiger fitDiscrete
Diet_ER = fitDiscrete(phy_inf, diet_phy, model = "ER")
Diet_SYM = fitDiscrete(phy_inf, diet_phy, model = "SYM")
Diet_ARD = fitDiscrete(phy_inf, diet_phy, model = "ARD")

# Model selection
#The smaller the better: ER
Diet_ER$opt$aicc
Diet_SYM$opt$aicc
Diet_ARD$opt$aicc

# simmap (Stochastic character mapping) with best model
mtree <- make.simmap(tree = ladderize(phy_inf, right = FALSE), x = diet_phy, 
                     model = "ER", nsim = 1000)
pd <- describe.simmap(mtree)
col.pd <- col.diet[colnames(pd$tips)]
write.table(pd$ace, "Tables/anc_ER_all.csv", 
            sep = ";", dec = ".", 
            col.names = NA, row.names = TRUE)
rownames(pd$ace)
Nnode(phy_inf)

pd$ace[1:Nnode(phy_inf),]
rownames(pd$ace)[1:Nnode(phy_inf)]
node_prob <- setNames(rep(0, Nnode(phy_inf)), nm = rownames(pd$ace)[1:Nnode(phy_inf)])
names(node_prob)

node <- '186'
pd$ace[node, ]

for (node in names(node_prob)){
  probs_sorted <- sort(pd$ace[node, ], decreasing = TRUE)
  if (probs_sorted[1] > probs_sorted[2]){
    node_prob[node] <- probs_sorted[1] 
  }
  else 
    node_prob[node] <- NA
}


pdf("plots/anc_ER_all.pdf", width = 8.27, height = 11.69)
plot(pd, fsize = 0.3, ftype = "i", colors = col.pd, lwd = 0.5, cex = 0.2)
#nodelabels(node=node_sig, pie=pd$ace[as.character(node_sig),], piecol=col.diet.inf, cex=0.2, frame="none", lty=par(lty="solid") )
nodelabels(round(node_prob, digits = 2), 
           frame = "none", adj = c(1.2,-0.6), cex = 0.5)
#tiplabels(pie=pd$tips,piecol=col.pd, cex=0.2, frame="none" )
legend("topleft", names(col.pd), cex = 0.8, col = col.diet[names(col.pd)], 
       pch = 19, box.lty = "blank", bty = "o")
title(main = "ER diet reconstruction", line = -1)
dev.off()




#As the dietary preferences of the lineage of the so called “ground squirrels” (Xerinae) 
#is very different to the rest of the family (known as “tree squirrels”), we inferred 
#the evolution of the dietary categories separately for those clades.

# ++ arboreal squirrels ----
Diet_ER_arb = fitDiscrete(tree_arb, diet_phy_arb, model="ER")
Diet_SYM_arb = fitDiscrete(tree_arb, diet_phy_arb, model="SYM")
Diet_ARD_arb = fitDiscrete(tree_arb, diet_phy_arb, model="ARD")

# Model selection
#The smaller the better: 
Diet_ER_arb$opt$aicc # best model
Diet_SYM_arb$opt$aicc
Diet_ARD_arb$opt$aicc


# simmap (Stochastic character mapping) with best model
mtree_arb <- make.simmap(tree = ladderize(tree_arb, right = FALSE), 
                         x = droplevels(diet_phy_arb), 
                     model = "ER", nsim = 1000)
pd_arb <- describe.simmap(mtree_arb)
col.pd_arb <- col.diet[colnames(pd_arb$tips)]
write.table(pd_arb$ace, "Tables/anc_ER_arb.csv", 
            sep = ";", dec = ".", 
            col.names = NA, row.names = TRUE)

node_prob_arb <- setNames(rep(0, Nnode(tree_arb)), 
                          nm = rownames(pd_arb$ace)[1:Nnode(tree_arb)])

for (node in names(node_prob_arb)){
  probs_sorted <- sort(pd_arb$ace[node, ], decreasing = TRUE)
  if (probs_sorted[1] > probs_sorted[2]){
    node_prob_arb[node] <- probs_sorted[1] 
  }
  else 
    node_prob_arb[node] <- NA
}


pdf("plots/anc_ER_arb.pdf", width = 8.27, height = 11.69)
plot(pd_arb, fsize = 0.3, ftype = "i", colors = col.pd_arb, lwd = 0.5, cex = 0.2)
#nodelabels(node=node_sig_arb, pie=pd_arb$ace[as.character(node_sig_arb),], piecol=col.diet, cex=0.2, frame="none", lty=par(lty="solid") )
nodelabels(round(node_prob_arb, digits = 2), 
           frame = "none", adj = c(1.2,-0.6), cex = 0.5)
#tiplabels(pie=pd$tips,piecol=col.pd, cex=0.2, frame="none" )
legend("topleft", names(col.pd_arb), cex = 0.8, col = col.diet[names(col.pd_arb)], 
       pch = 19, box.lty = "blank", bty = "o")
title(main = "ER diet reconstruction", line = -1)
dev.off()





# ++ terrestrial squirrels ----
Diet_ER_terr = fitDiscrete(tree_terr, diet_phy_terr, model="ER")
Diet_SYM_terr = fitDiscrete(tree_terr, diet_phy_terr, model="SYM")
Diet_ARD_terr = fitDiscrete(tree_terr, diet_phy_terr, model="ARD")

# Model selection
#The smaller the better: 
Diet_ER_terr$opt$aicc
Diet_SYM_terr$opt$aicc
Diet_ARD_terr$opt$aicc


# simmap (Stochastic character mapping) with best model
mtree_terr <- make.simmap(tree = ladderize(tree_terr, right = FALSE), x = droplevels(diet_phy_terr), 
                         model = "SYM", nsim = 1000)
pd_terr <- describe.simmap(mtree_terr)
col.pd_terr <- col.diet[colnames(pd_terr$tips)]
write.table(pd_terr$ace, "Tables/anc_SYM_terr.csv", 
            sep = ";", dec = ".", 
            col.names = NA, row.names = TRUE)

node_prob_terr <- setNames(rep(0, Nnode(tree_terr)), 
                          nm = rownames(pd_terr$ace)[1:Nnode(tree_terr)])

for (node in names(node_prob_terr)){
  probs_sorted <- sort(pd_terr$ace[node, ], decreasing = TRUE)
  if (probs_sorted[1] > probs_sorted[2]){
    node_prob_terr[node] <- probs_sorted[1] 
  }
  else 
    node_prob_terr[node] <- NA
}



pdf("plots/anc_SYM_terr.pdf", width = 8.27, height = 11.69)
plot(pd_terr, fsize = 0.3, ftype = "i", colors = col.pd_terr, lwd = 0.5, cex = 0.2)
#nodelabels(node=node_sig_terr, pie=pd_terr$ace[as.character(node_sig_terr),], piecol=col.diet, cex=0.2, frame="none", lty=par(lty="solid") )
nodelabels(round(node_prob_terr, digits = 2), 
           frame = "none", adj = c(1.2,-0.6), cex = 0.5)
#tiplabels(pie=pd$tips,piecol=col.pd, cex=0.2, frame="none" )
legend("topleft", names(col.pd_terr), cex = 0.8, col = col.diet[names(col.pd_terr)], 
       pch = 19, box.lty = "blank", bty = "o")
title(main = "SYM diet reconstruction", line = -1)
dev.off()


# ++ aicc table
AIC_transition <- data.frame(All= c(Diet_ER$opt$aicc, Diet_SYM$opt$aicc, Diet_ARD$opt$aicc), 
                             Tree= c(Diet_ER_arb$opt$aicc, Diet_SYM_arb$opt$aicc, Diet_ARD_arb$opt$aicc),
                             Ground= c(Diet_ER_terr$opt$aicc, Diet_SYM_terr$opt$aicc, Diet_ARD_terr$opt$aicc))
row.names(AIC_transition) <- c("ER", "SYM", "ARD")
write.csv(AIC_transition, "Tables/AICc_transition.csv")


# PLOTS ----

# multiplying factor for Csize
m <- 0.17

# + PLOT LDA GM ----
projection_GM <- scale(as.matrix(PCA_extant_GM$x[,1:6]), 
                       scale = FALSE) %*% LDA_GM$scaling
projection_GM <- data.frame(projection_GM, 
                            diet = fundiet_extant[rownames(PCA_extant_GM$x)],  
                            stringsAsFactors = FALSE)
#rbind(data.frame(irisProjection_l, predict="pca"), data.frame(lda_landmark_predictPCA$x, diet=lda_landmark_predictPCA$class, predict= "predicted", stringsAsFactors = F))

plot_LDA_GM <- ggplot(data = projection_GM, aes(x = -LD1, y = -LD2))

plot1 <- plot_LDA_GM + 
  geom_point(data = projection_GM, aes(col = diet), 
             size = Csize_extant*m,
             alpha = 4/10, shape=21, stroke=0.5, show.legend = FALSE) + 
  scale_color_manual(values = col.fundiet)+
  geom_point(data = as.data.frame(predicted_fossil_LDA_GM$x), 
             aes(color = predicted_fossil_LDA_GM$class), 
             size = Csize_fossil*m, 
             shape = 16, show.legend = FALSE, alpha = 8/10)+
  annotate("text",x = 3.5, y=-4, 
           label = paste0('classification\ncorrectness: ', 
                          round(CV_GM_perc_total*100, 1), 
                          '%'), 
           size = 3.5, hjust = 0) +
  labs(title = "GM LDA", x = 'LD1', y = 'LD2') +
  xlim(c(-4,7)) +
  ylim(c(-6,4)) +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold'))



# + PLOT LDA EFA ----

projection_EFA <- scale(as.matrix(PCA_extant_EFA$x[,1:6]), 
                        scale=FALSE) %*% LDA_EFA$scaling
projection_EFA <- data.frame(projection_EFA, 
                             diet = as.factor(fundiet_extant[shape_spred_extant]), 
                             stringsAsFactors = FALSE)
#rbind(data.frame(irisProjection_f, predict="pca"), data.frame(relda.fourierPCA$x, diet=relda.fourierPCA$class, predict= "predicted", stringsAsFactors = F))

plot_LDA_EFA <- ggplot(data = projection_EFA, aes(x = -LD1, y = -LD2))

plot2 <- plot_LDA_EFA + 
  geom_point(data = projection_EFA, aes(col = diet), 
             size = Csize_extant*m,
             alpha = 4/10, shape=21, stroke=0.5, show.legend = FALSE) + 
  scale_color_manual(values = col.fundiet)+
  geom_point(data = as.data.frame(predicted_fossil_LDA_EFA$x), 
             aes(color = predicted_fossil_LDA_EFA$class), 
             size = Csize_fossil*m, 
             shape = 16, show.legend = FALSE, alpha = 8/10)+
  annotate("text",x = 3.5, y=-4, 
           label = paste0('classification\ncorrectness: ', 
                          round(CV_EFA_perc_total*100, 1), 
                          '%'), 
           size = 3.5, hjust = 0) +
  labs(title = "EFA LDA", x = 'LD1', y = 'LD2') +
  xlim(c(-4,7)) +
  ylim(c(-6,4)) +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold'))



# + PLOT LDA GM + size ----

projection_GM_size <- scale(as.matrix(PCA_extant_GM_size$x[,1:6]), 
                            scale=FALSE) %*% LDA_GM_size$scaling
projection_GM_size <- data.frame(projection_GM_size, 
                                 diet = fundiet_extant[shape_spred_extant],  
                                 stringsAsFactors = FALSE)
rbind(data.frame(projection_GM_size, predict="pca"), 
      data.frame(predicted_fossil_LDA_GM_size$x, 
                 diet=predicted_fossil_LDA_GM_size$class, 
                 predict= "predicted", 
                 stringsAsFactors = F))

plot_LDA_GM_size <- ggplot(data = projection_GM_size, aes(x=-LD1,y=-LD2))

plot3 <- plot_LDA_GM_size + 
  geom_point(data = projection_GM_size, aes(col = diet), 
             size = Csize_extant*m,
             alpha = 4/10, shape=21, stroke=0.5, show.legend = FALSE) + 
  scale_color_manual(values = col.fundiet)+
  geom_point(data = as.data.frame(predicted_fossil_LDA_GM_size$x), 
             aes(color = predicted_fossil_LDA_GM_size$class), 
             size = Csize_fossil*m, 
             shape = 16, show.legend = FALSE, alpha = 8/10)+
  annotate("text",x = 3.5, y=-4, 
           label = paste0('classification\ncorrectness: ', 
                          round(CV_GM_size_perc_total*100, 1), 
                          '%'), 
           size = 3.5, hjust = 0) +
  labs(title = "GM + size LDA", x = 'LD1', y = 'LD2') +
  xlim(c(-4,7)) +
  ylim(c(-6,4)) +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold'))

# + PLOT LDA EFA + size ----
projection_EFA_size <- scale(as.matrix(PCA_extant_EFA_size$x[,1:6]), 
                        scale=FALSE) %*% LDA_EFA_size$scaling
projection_EFA_size <- data.frame(projection_EFA_size, 
                             diet = as.factor(fundiet_extant[shape_spred_extant]), 
                             stringsAsFactors = FALSE)

plot_LDA_EFA_size <- ggplot(data = projection_EFA_size, aes(x = -LD1, y = -LD2))

plot4 <- plot_LDA_EFA_size + 
  geom_point(data = projection_EFA_size, aes(col = diet), 
             size = Csize_extant*m,
             alpha = 4/10, shape=21, stroke=0.5, show.legend = FALSE) + 
  scale_color_manual(values = col.fundiet)+
  geom_point(data = as.data.frame(predicted_fossil_LDA_EFA_size$x), 
             aes(color = predicted_fossil_LDA_EFA_size$class), 
             size = Csize_fossil*m, 
             shape = 16, show.legend = FALSE, alpha = 8/10)+
  annotate("text",x = 3.5, y=-4, 
           label = paste0('classification\ncorrectness: ', 
                          round(CV_EFA_size_perc_total*100, 1), 
                          '%'), 
           size = 3.5, hjust = 0) +
  labs(title = "EFA + size LDA", x = 'LD1', y = 'LD2') +
  xlim(c(-4,7)) +
  ylim(c(-6,4)) +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold'))

# PLOT ALL LDAs 
plot_list <- list(plot1, plot2, plot3, plot4)
multiplot <- plot_grid(plotlist = plot_list, ncol =2, labels = 'auto')
ggsave("Figures/Figure_4_v4.pdf", multiplot, height = 7, width = 11)



# + PLOT FOSSILS OVER TIME ----
head(extinct_names)
fossil_dat <- data.frame(specimen = extinct_names$specimen_full, 
                         Sp = extinct_names$sp_short, 
                         AcceptedSp = extinct_names$sp_full,
                         minAge = extinct_names$Min_age, 
                         maxAge = extinct_names$Max_age,
                         stringsAsFactors = F)
#fossil_dat <- fossil_dat[match(new_class_fourierPCA_size$sp,fossil_dat$Sp),]
fossil_dat <- fossil_dat[match(as.character(new_class_GM_size$sp), fossil_dat$specimen),]

new_class <- new_class_GM_size
new_class$Csize <- Csize_fossil
new_class$minAge <- fossil_dat$minAge
new_class$maxAge <- fossil_dat$maxAge
new_class$sp_name <- fossil_dat$Sp
colnames(new_class)[2] <- "diet_class"
new_class




write.csv(new_class, "new_class_fossils2.csv")



new_class$Mean <- (new_class$maxAge+new_class$minAge)/2
new_class_post <- new_class[new_class$posterior > 0.75,]


pdf("plots/fossils_time.pdf", paper = "a4")
#4, 13, 16, 18, 20, 28, 40, 41, 44, 76
ggplot() +
  geom_point(data = new_class_post, aes(x = Mean, y = diet_class, 
                                        group = diet_class),
             color = col.fundiet[new_class_post$diet_class], 
             alpha = new_class_post$posterior, 
             shape = 16, size = new_class_post$Csize*0.15, 
             position =  position_jitter(width = 1, height = 0.5, seed = 20)) +
  geom_text(data = new_class_post, aes(x = Mean, y = diet_class, 
                                       label = sp_name),
            position = position_jitter(width = 1, height = 0.5, seed = 20), 
            size = 0.8) +
  scale_x_reverse() +
  coord_geo(xlim = c(45, 0), ylim = c(0,8), dat ="epochs", 
            height = unit(1, "lines"), size = 2, abbrv = FALSE) +
  theme_classic()
dev.off()









