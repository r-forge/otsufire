% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2.process_otsu_rasters.R
\name{process_otsu_rasters}
\alias{process_otsu_rasters}
\title{Process Burned Area Rasters Using Otsu Thresholding or Percentile Clipping}
\arguments{
\item{raster_path}{Path to a single-band RBR or dNBR raster.}

\item{nbr_pre_path, nbr_post_path}{Optional. Paths to pre- and post-fire NBR rasters for RBR/dNBR calculation.}

\item{output_dir}{Directory to save output files.}

\item{year}{Optional year label.}

\item{otsu_thresholds}{Numeric vector of minimum values to filter raster before applying Otsu.}

\item{trim_percentiles}{Optional data.frame with `min` and `max` columns to define percentile clipping thresholds.}

\item{use_original}{Logical. Use raw values without filtering (ignored if `trim_percentiles` is set).}

\item{corine_raster_path}{Path to CORINE raster for stratified segmentation.}

\item{reclassify_corine}{Logical. If TRUE, reclassifies CORINE using `reclass_matrix`.}

\item{reclass_matrix}{Matrix with original and new values for CORINE reclassification.}

\item{peninsula_shapefile}{Shapefile used to crop CORINE before reclassification.}

\item{output_corine_raster_dir}{Directory to save reclassified CORINE raster.}

\item{output_corine_vector_dir}{Directory to save reclassified CORINE shapefile.}

\item{reproject}{Logical. Reproject reclassified CORINE raster to EPSG:3035.}

\item{resolution}{Numeric. Target resolution for reprojected CORINE raster.}

\item{corine_classes}{Optional vector of reclassified CORINE classes to keep.}

\item{ecoregion_shapefile_path}{Path to WWF ecoregions shapefile.}

\item{ecoregion_field}{Column in ecoregion shapefile used for class labels.}

\item{ecoregion_classes}{Optional. Vector of selected ecoregion class names.}

\item{segment_by_intersection}{Logical. If TRUE, intersects CORINE and ecoregions.}

\item{min_otsu_threshold_value}{Minimum acceptable threshold. Used if Otsu value is lower.}

\item{python_exe}{Path to Python executable.}

\item{gdal_polygonize_script}{Path to `gdal_polygonize.py` script.}

\item{gdalwarp_path}{Path to `gdalwarp` executable (default = "gdalwarp").}

\item{n_rows, n_cols}{Number of rows/columns for raster tiling.}

\item{tile_overlap}{Overlap buffer (in map units) for tiles.}

\item{tile}{Logical. Apply raster tiling before polygonization.}

\item{index_type}{"RBR" or "dNBR". Used if `nbr_pre_path`/`nbr_post_path` are provided.}

\item{output_format}{Output vector file format: "shp" or "geojson".}
}
\value{
Named list of `sf` polygons grouped by segmentation. Output files saved to `output_dir`.
}
\description{
This function computes binary burned area masks from a severity index raster (RBR or dNBR)
using Otsu's thresholding or percentile clipping. The segmentation can be applied
to the full raster or stratified by:
- CORINE land cover classes (`corine_raster_path`)
- WWF ecoregions (`ecoregion_shapefile_path`)
- or the intersection of both (CORINE ? Ecoregion)

## Key features:
- Supports Otsu thresholding after pre-filtering by minimum index values (`otsu_thresholds`) or using original values
- Supports percentile-based thresholding via `trim_percentiles`
- Enforces a minimum threshold (`min_otsu_threshold_value`) when Otsu yields low values
- Allows stratification by CORINE, ecoregions, or their intersection
- Enables CORINE reclassification using a custom `reclass_matrix`
- Automatically computes RBR or dNBR if `nbr_pre_path` and `nbr_post_path` are provided
- Outputs:
  - Burned area binary rasters
  - Polygon shapefiles (dissolved by threshold label)
  - Histogram and inter-class variance plots per threshold
  - Threshold log files
}
\details{
- Raster values are internally rescaled to [0, 255] before Otsu thresholding.
- Histogram smoothing and variance curves are used for enhanced threshold detection.
- Output shapefiles can be in ESRI Shapefile or GeoJSON format.
- Intermediate tile shapefiles and rasters are cleaned after processing.
}
\note{
Examples require large external raster files (hosted on Zenodo)
and depend on external software (Python, GDAL). Therefore, they are wrapped
in dontrun{} to avoid errors during R CMD check and to ensure portability.
}
\examples{
\dontrun{
# Example 1: Basic Otsu thresholds (global)
process_otsu_rasters(
  raster_path = "data/rbr_1985.tif",
  output_dir = "output/1985",
  otsu_thresholds = c(0, 50, 100),
  python_exe = "C:/Python/python.exe",
  gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
)

# Example 2: Global with percentile trimming
process_otsu_rasters(
  raster_path = "data/rbr_1985.tif",
  output_dir = "output/1985",
  trim_percentiles = data.frame(min = 0.01, max = 0.99),
  python_exe = "C:/Python/python.exe",
  gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py",
  output_format = "geojson"
)

# Example 3: Stratified by reclassified CORINE
process_otsu_rasters(
  raster_path = "data/rbr_1985.tif",
  output_dir = "output/1985",
  corine_raster_path = "data/corine_1990_etrs89.tif",
  reclassify_corine = TRUE,
  reclass_matrix = matrix(c(
    22, 1,  # grasslands ? 1
    26, 1,
    32, 1,
    27, 2,  # shrublands ? 2
    29, 2,
    33, 3,  # burnt areas ? 3
    23, 4,  # broadleaved forest ? 4
    25, 5,  # mixed forest ? 5
    24, 6   # coniferous ? 6
  ), ncol = 2, byrow = TRUE),
  corine_classes = c(1, 2, 4, 5, 6),
  peninsula_shapefile = "data/peninsula_clip.shp",
  output_corine_raster_dir = "output/corine",
  output_corine_vector_dir = "output/corine",
  otsu_thresholds = c(50),
  min_otsu_threshold_value = 150,
  python_exe = "C:/Python/python.exe",
  gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
)

# Example 4: Stratified by WWF ecoregions (default field: "EnS_name")
process_otsu_rasters(
  raster_path = "data/rbr_1985.tif",
  output_dir = "output/1985",
  ecoregion_shapefile_path = "data/ecoregions_olson.shp",
  otsu_thresholds = c(50),
  min_otsu_threshold_value = 150,
  python_exe = "C:/Python/python.exe",
  gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
)

# Example 5: Stratified by WWF ecoregions with custom field
process_otsu_rasters(
  raster_path = "data/rbr_1985.tif",
  output_dir = "output/1985",
  ecoregion_shapefile_path = "data/ecoregions_olson.shp",
  ecoregion_field = "EnZ_name",
  ecoregion_classes = c("Iberian Conifer Forests", "Mediterranean Shrublands"),
  otsu_thresholds = c(50),
  min_otsu_threshold_value = 150,
  python_exe = "C:/Python/python.exe",
  gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
)

# Example 6: Stratified by CORINE ? Ecoregion intersection
process_otsu_rasters(
  raster_path = "data/rbr_1985.tif",
  output_dir = "output/1985",
  corine_raster_path = "data/corine_1990_etrs89.tif",
  reclassify_corine = TRUE,
  reclass_matrix = matrix(c(
    22, 1, 26, 1, 32, 1,  # grasslands
    27, 2, 29, 2,         # shrublands
    23, 4, 25, 5, 24, 6   # forests
  ), ncol = 2, byrow = TRUE),
  corine_classes = c(1, 2, 4, 5, 6),
  peninsula_shapefile = "data/peninsula_clip.shp",
  output_corine_raster_dir = "output/corine",
  output_corine_vector_dir = "output/corine",
  ecoregion_shapefile_path = "data/ecoregions_olson.shp",
  ecoregion_field = "EnZ_name",
  segment_by_intersection = TRUE,
  otsu_thresholds = c(50),
  min_otsu_threshold_value = 150,
  python_exe = "C:/Python/python.exe",
  gdal_polygonize_script = "C:/Python/Scripts/gdal_polygonize.py"
)
}
}
