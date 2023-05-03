filter_by_element <- function(x, element) {
  x <- lapply(x, function(y) {
    select <- y$taxon.list$variable.element == element
    y$taxon.list <- y$taxon.list[select,]
    y$counts <- y$counts[,colnames(y$counts) %in% y$taxon.list$taxon.name]
    return(y)
  })
  class(x) <- append(class(x), "download_list")
  return(x)
}

filter_by_unit <- function(x, unit) {
  x <- lapply(x, function(y) {
    select <- y$taxon.list$variable.unit == unit
    y$taxon.list <- y$taxon.list[select,]
    y$counts <- y$counts[,colnames(y$counts) %in% y$taxon.list$taxon.name]
    return(y)
  })
  class(x) <- append(class(x), "download_list")
  return(x)
}

transform_to_percent <- function(x)
{
  x$counts <- by(x, x$sampleid, function(y)
    y$value/sum(y$value) * 100
    ) %>% unlist()

  class(x) <- append(class(x), "download_list")
  return(x)
}

aggregate_by_joint_name <- function(x, common_names) {

  x <- merge(pollen_data_comp, common_names, by.x = "variablename", by.y = "name_orig")



  counts.collected <- x %>%
    as_tibble() %>%
    group_by(name_joint, sampleid) %>%
    summarise(counts = sum(value))


  ages <- distinct(x, sampleid, age)

  sites <- distinct(x, sampleid, sitename)

  counts.new <- pivot_wider(counts.collected, names_from = name_joint, values_from = counts)

  counts.new <- merge(counts.new, ages, by = "sampleid")

  counts.new <- merge(counts.new, sites, by = "sampleid")

  return(counts.new)
}
