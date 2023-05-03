base_map <- function(crs = 2056, water = "lightblue4",
                     study_region = "darkred", highland = "lightgray") {
  ggplot() +
    layer_spatial(read_sf(here("data", "raw_data", "eudem_over1000m_ch.gpkg")),
                  fill = highland, color = NA) +
    layer_spatial(study_region(),
                  fill = alpha(study_region, 0.33), color = NA) +
    layer_spatial(ne("ne_10m_rivers_lake_centerlines"), color = water) +
    layer_spatial(ne("ne_10m_rivers_europe"), color = water) +
    layer_spatial(ne("ne_10m_lakes"), fill = water, color = NA) +
    layer_spatial(ne("ne_10m_lakes_europe"), fill = water, color = NA) +
    annotation_scale(location = "br") +
    base_coord() +
    theme_minimal()
}

study_region <- function(bioregion_data = here("data", "fixed_online_data",
                                               "BiogeographischeRegionen.gdb")) {
  sf::st_read(bioregion_data, "N2020_Revision_BiogeoRegion", quiet = TRUE) |>
    sf::st_transform(4326) |>
    dplyr::filter(RegionNummer == 2) |>
    st_union()
}

ne <- function(layer, ne_data = here("data", "raw_data", "ne")) {
  sf::st_read(ne_data, layer, quiet = TRUE)
}

base_coord <- function(crs = 2056) {
  ggplot2::coord_sf(
    xlim = c(6, 10.5),
    ylim = c(45.5, 48),
    crs = crs,
    default_crs = 4326
  )
}
