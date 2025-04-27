library(sp)

# Define the extent (bounding box) of the grid
bbox <- matrix(c(3.2, 40.7, 13.7, 45.1), ncol = 2, byrow = TRUE)
colnames(bbox) <- c("min", "max")
rownames(bbox) <- c("x", "y")


# Set the resolution of the grid (size of each grid cell)
cell_size <- c(0.5, 0.5)  # Each cell will be 1x1 unit

# Create the grid
grid <- GridTopology(cellcentre.offset = bbox[, "min"], cellsize = cell_size, cells.dim = c(10, 10))

# Convert the grid into a SpatialGrid object
spatial_grid <- SpatialGrid(grid)

# Plot the grid
plot(spatial_grid, main = "Regular Grid Created with sp")

dati <- runif(length(spatial_grid), min = -1, max = 1)  # Dati casuali tra 0 e 100

# Associare i dati alla griglia utilizzando SpatialGridDataFrame
grid_data <- SpatialGridDataFrame(spatial_grid, data = data.frame(valore = dati))

# Visualizziamo la griglia con i dati associati
plot(grid_data, main = "Griglia con Dati Associati")


# 2. Creare una griglia a risoluzione più elevata
cell_size_fine <- c(0.1, 0.1)  # Celle di 0.5x0.5
grid_fine <- GridTopology(cellcentre.offset = bbox[, "min"], cellsize = cell_size_fine, cells.dim = c(20, 20))
spatial_grid_fine <- SpatialGrid(grid_fine)

grid_fine
sp::makegrid()

# 3. Interpolazione dei valori sui nuovi punti della griglia (bilineare)
# Otteniamo le coordinate della nuova griglia
coords_fine <- coordinates(spatial_grid_fine)

# Creiamo una funzione di interpolazione per i dati sulla griglia originale
interpolazione <- sp::over(spatial_grid_fine, grid_data, fn = mean)

# Creare il nuovo SpatialGridDataFrame con i valori interpolati
grid_data_fine <- SpatialGridDataFrame(spatial_grid_fine, data = data.frame(valore = interpolazione))

SpatialGridDataFrame(spatial_grid_fine, data = Bologna)

# 4. Visualizzare i risultati della griglia interpolata
spplot(grid_data_fine, "valore", main = "Griglia Interpolata a Maggiore Risoluzione")

