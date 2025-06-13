% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.flag_otsu_regenera.R
\encoding{UTF-8}
\name{flag_otsu_regenera}
\alias{flag_otsu_regenera}
\title{Flag Burned Polygons Based on Regeneration Overlap}
\arguments{
\item{burned_files}{Character vector with paths to burned area shapefiles (.shp).}

\item{regenera_files}{Character vector with paths to regeneration shapefiles (.shp).
Filenames must contain `"P1"`, `"P2"`, etc., to indicate the year after fire.}

\item{min_regen_ratio}{Named numeric vector with minimum area ratio thresholds (between 0 and 1)
  for each regeneration year. Names must correspond to period labels (e.g., `c(P1 = 0.05, P2 = 0.10)`).
 These values are used:
(1) when no classwise_regen_thresholds are provided,
(2) or as fallback when a polygon has no match or NA thresholds.}

\item{remove_no_regenera}{Logical. If `TRUE`, polygons that do not meet all specified
thresholds (see `remove_condition`) will be excluded from the filtered output.}

\item{remove_condition}{Character vector indicating the post-fire years (e.g., `"P1"`, `"P2"`, `"P3"`)
that must meet their thresholds in `min_regen_ratio` in order for a polygon to be retained.
Used to identify false fire detections due to lack of regeneration.}

\item{output_dir}{Optional output directory. If `NULL`, output files are saved in the same
folder as the burned shapefiles.}

\item{replace_by_P1}{Logical. If `TRUE` (default), filtered burned polygons will be replaced by
overlapping P1 regeneration polygons when the latter are larger in area.}

\item{max_area_ratio_p1}{Optional. Maximum allowed ratio between the area of a P1 regeneration polygon and the intersected burned polygon.
If specified, only regenerated polygons with area greater than the burned area and not exceeding \code{max_area_ratio_p1} times that area
will be used to replace burned polygons. If \code{NULL} (default), only the condition \code{area_ha > burn_area_ha} is applied.}

\item{save_no_regenera}{Character. Options are `"all"` (save all `"no_regenera"` polygons),
`"area_filter"` (save only those with area ? `min_area_no_regenera`), or `"none"` (do not save them).}

\item{min_area_no_regenera}{Numeric. Minimum area in hectares to retain a `"no_regenera"` polygon when `save_no_regenera = "area_filter"`.
Default is `0.1` ha.}

\item{validate_geometries}{Logical. If `TRUE`, applies `st_make_valid()` to correct invalid geometries
before processing. Default is `FALSE` for faster performance.}

\item{output_format}{Character. Output format for saved files: `"shp"` for ESRI Shapefile (default)
  or `"geojson"` for GeoJSON.
 @param classwise_regen_thresholds Optional `data.frame` with per-class regeneration thresholds.
Must include:
\itemize{
  \item \code{class_field}: field name in the burned polygons (e.g., `"CORINE_CLA"`, `"eco_id"`).
  \item \code{class_value}: corresponding class value.
  \item One or more columns named \code{P1}, \code{P2}, etc., with numeric thresholds.
}
If provided, these thresholds will override \code{min_regen_ratio} for polygons matching each class.
Polygons without a match will fall back to the values in \code{min_regen_ratio}.}
}
\value{
Saves up to four shapefiles per burned area file:
\describe{
  \item{Main shapefile (unfiltered)}{All polygons with `regen_ratio_Px`, `regenera_flag_Px`, and `regenera_flag_all`.}
  \item{Filtered shapefile}{Only polygons that satisfy all thresholds defined in `remove_condition`. File ends in `_filter.shp`.}
  \item{Final shapefile with replacements}{(if `replace_by_P1 = TRUE`) Filtered polygons with some replaced by larger P1 regeneration polygons. File ends in `_new.shp`.}
  \item{Shapefile of no_regenera polygons}{(if `save_no_regenera != "none"`) Contains polygons not meeting regeneration thresholds. File ends in `_no_selected.shp`.}
}
}
\description{
This function assigns regeneration flags to burned area polygons based on their spatial
overlap with post-fire regeneration polygons (e.g., from years 1, 2, or more after fire).
A polygon is flagged as `"regenera"` if the total intersected area with regeneration polygons
from a specific post-fire year exceeds a user-defined minimum threshold (`min_regen_ratio`)
relative to its own area.

Filenames of the regeneration layers must contain `"P1"`, `"P2"`, etc., to indicate
the number of years after fire. These period labels are automatically extracted and used
to name the output columns (e.g., `regen_ratio_P1`, `regenera_flag_P1`).

For each detected regeneration year, the function adds:
\itemize{
  \item `regen_ratio_Px`: proportion of the burned polygon area overlapping with regeneration polygons.
  \item `regenera_flag_Px`: `"regenera"` if the threshold is met, `"no_regenera"` otherwise.
}
A final column `regenera_flag_all` is also created, indicating `"regenera"` only if
all selected years in `remove_condition` meet their respective thresholds.

If `replace_by_P1 = TRUE`, the function will optionally replace burned polygons in the
filtered output by overlapping P1 regeneration polygons that are larger in area. Optionally,
a maximum area ratio (`max_area_ratio_p1`) can be defined to restrict replacements to only
those P1 polygons not exceeding a given multiple of the burned area.

You can also choose to export polygons flagged as `"no_regenera"` using `save_no_regenera`.
These can be optionally filtered by a minimum area threshold (`min_area_no_regenera`).

If `validate_geometries = TRUE`, geometries are checked and corrected using `sf::st_make_valid()`,
which is more robust but significantly slower. Set to `FALSE` by default for performance.
}
\note{
Examples require large external raster files (hosted on Zenodo).
Therefore, they are wrapped in dontrun{} to avoid errors during R CMD check
and to ensure portability.
}
\examples{
\dontrun{
# Apply per-period regeneration thresholds and filter non-regenerating polygons
burned_files <- list.files("data/burned", pattern = "\\\\.shp$", full.names = TRUE)
regenera_files <- list.files("data/regenera", pattern = "\\\\.shp$", full.names = TRUE)

flag_otsu_regenera(
  burned_files = burned_files,
  regenera_files = regenera_files,
  min_regen_ratio = c(P1 = 0.05, P2 = 0.10),
  remove_no_regenera = TRUE,
  remove_condition = c("P1", "P2"),
  output_dir = "results/",
  replace_by_P1 = TRUE,
  max_area_ratio_p1 = 3,
  save_no_regenera = "area_filter",
  min_area_no_regenera = 100,
  validate_geometries = FALSE,
  output_format = c("geojson")
)

# Define class-specific regeneration thresholds
class_thresholds <- data.frame(
  class_field = "CORINE_CLA",
  class_value = c(1, 2),
  P1 = c(0.15, 0.15),
  P2 = c(0.0, 0.05)
)

# Run flag_otsu_regenera1 using classwise thresholds and fallback values
result <- flag_otsu_regenera(
  burned_files = burned_list,
  regenera_files = regenera_list,
  min_regen_ratio = c(P1 = 0.05, P2 = 0.15),  # fallback for classes not in the table
  classwise_regen_thresholds = class_thresholds,
  group_field = "CORINE_CLA",
  remove_no_regenera = TRUE,
  remove_condition = c("P1", "P2"),
  replace_by_P1 = TRUE,
  max_area_ratio_p1 = 3,
  save_no_regenera = TRUE,
  min_area_no_regenera = 500,
  validate_geometries = FALSE,
  output_format = "shp",
  output_dir = burned_dir
)
}

}
