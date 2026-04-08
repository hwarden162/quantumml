library(stringr)
library(tidyverse)

read_folder <- function(path) {
  morph_data_path <- paste0("./raw_data/", path, "/morphology.csv")
  spati_data_path <- paste0("./raw_data/", path, "/spatial.csv")
  suppressMessages({
    morph_data <- read_csv(morph_data_path)
    spati_data <- read_csv(spati_data_path)
  })
  gof <- str_split(path, "_")[[1]][2]
  return_df <- morph_data |> 
    left_join(spati_data) |> 
    filter(
      AreaShape_Nuclei_Mask_AxisMinorLength > 6.25,
      Spatial_Nuclei_Spatial_LocalCounts100 > 5,
      AreaShape_Nuclei_Mask_Area |> between(20, 200),
      !is.na(Spatial_Object_Spatial_LocalMeansIntensityCytoplasmSecondaryMeanIntensity200)
    ) |> 
    mutate(
      cytoplasmic_stain_ratio = Intensity_Cytoplasm_Secondary_MedianIntensity / Spatial_Object_Spatial_LocalMeansIntensityCytoplasmSecondaryMeanIntensity200
    )
  csr_mean <- return_df |> 
    pull(cytoplasmic_stain_ratio) |> 
    mean()
  csr_uq <- return_df |> 
    pull(cytoplasmic_stain_ratio) |> 
    quantile(0.75)
  csr_lq <- return_df |> 
    pull(cytoplasmic_stain_ratio) |> 
    quantile(0.25)
  csr_iqr <- csr_uq - csr_lq
  return_df |> 
    filter(
      cytoplasmic_stain_ratio > 1.5,
      cytoplasmic_stain_ratio < csr_mean + 5*csr_iqr
    ) |> 
    arrange(-cytoplasmic_stain_ratio) |> 
    slice_head(n = 100) |> 
    mutate(
      path = path,
      GOF = gof
    )
}

process_df <- function(df) {
  df |> 
    mutate(
      Count_Disk1 = Spatial_Nuclei_Spatial_LocalCounts100,
      Count_Disk2 = Spatial_Nuclei_Spatial_LocalCounts200,
      Count_Disk3 = Spatial_Nuclei_Spatial_LocalCounts300,
      Count_Annulus1 = Count_Disk1,
      Count_Annulus2 = Count_Disk2 - Count_Disk1,
      Count_Annulus3 = Count_Disk3 - Count_Disk2,
      
      SumArea_Disk1 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskArea100 * Count_Annulus1,
      SumArea_Disk2 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskArea200 * Count_Annulus2,
      SumArea_Disk3 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskArea300 * Count_Annulus3,
      SumArea_Annulus1 = SumArea_Disk1,
      SumArea_Annulus2 = SumArea_Disk2 - SumArea_Disk1,
      SumArea_Annulus3 = SumArea_Disk3 - SumArea_Disk2,
      Area_Annulus0 = AreaShape_Nuclei_Mask_Area,
      Area_Annulus1 = SumArea_Annulus1 / Count_Annulus1,
      Area_Annulus2 = SumArea_Annulus2 / Count_Annulus2,
      Area_Annulus3 = SumArea_Annulus3 / Count_Annulus3,
      
      SumPerimeter_Disk1 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskPerimeter100 * Count_Annulus1,
      SumPerimeter_Disk2 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskPerimeter200 * Count_Annulus2,
      SumPerimeter_Disk3 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskPerimeter300 * Count_Annulus3,
      SumPerimeter_Annulus1 = SumPerimeter_Disk1,
      SumPerimeter_Annulus2 = SumPerimeter_Disk2 - SumPerimeter_Disk1,
      SumPerimeter_Annulus3 = SumPerimeter_Disk3 - SumPerimeter_Disk2,
      Perimeter_Annulus0 = AreaShape_Nuclei_Mask_Perimeter,
      Perimeter_Annulus1 = SumPerimeter_Annulus1 / Count_Annulus1,
      Perimeter_Annulus2 = SumPerimeter_Annulus2 / Count_Annulus2,
      Perimeter_Annulus3 = SumPerimeter_Annulus3 / Count_Annulus3,
      
      SumEccentricity_Disk1 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskEccentricity100 * Count_Annulus1,
      SumEccentricity_Disk2 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskEccentricity200 * Count_Annulus2,
      SumEccentricity_Disk3 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskEccentricity300 * Count_Annulus3,
      SumEccentricity_Annulus1 = SumEccentricity_Disk1,
      SumEccentricity_Annulus2 = SumEccentricity_Disk2 - SumEccentricity_Disk1,
      SumEccentricity_Annulus3 = SumEccentricity_Disk3 - SumEccentricity_Disk2,
      Eccentricity_Annulus0 = AreaShape_Nuclei_Mask_Eccentricity,
      Eccentricity_Annulus1 = SumEccentricity_Annulus1 / Count_Annulus1,
      Eccentricity_Annulus2 = SumEccentricity_Annulus2 / Count_Annulus2,
      Eccentricity_Annulus3 = SumEccentricity_Annulus3 / Count_Annulus3,
      
      SumAxisMinorLength_Disk1 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskAxisMinorLength100 * Count_Annulus1,
      SumAxisMinorLength_Disk2 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskAxisMinorLength200 * Count_Annulus2,
      SumAxisMinorLength_Disk3 = Spatial_Object_Spatial_LocalMeansAreaShapeNucleiMaskAxisMinorLength300 * Count_Annulus3,
      SumAxisMinorLength_Annulus1 = SumAxisMinorLength_Disk1,
      SumAxisMinorLength_Annulus2 = SumAxisMinorLength_Disk2 - SumAxisMinorLength_Disk1,
      SumAxisMinorLength_Annulus3 = SumAxisMinorLength_Disk3 - SumAxisMinorLength_Disk2,
      AxisMinorLength_Annulus0 = AreaShape_Nuclei_Mask_AxisMinorLength,
      AxisMinorLength_Annulus1 = SumAxisMinorLength_Annulus1 / Count_Annulus1,
      AxisMinorLength_Annulus2 = SumAxisMinorLength_Annulus2 / Count_Annulus2,
      AxisMinorLength_Annulus3 = SumAxisMinorLength_Annulus3 / Count_Annulus3,
    ) |> 
    select(
      path,
      Meta_Global_Mask_Label,
      GOF,
      contains("Annulus"),
      -starts_with("Sum")
    )
}

raw_dirs <- list.dirs("./raw_data", full.names = FALSE)[-1]

raw_dirs |> 
  map(
    read_folder
  ) |> 
  bind_rows() |> 
  process_df() |> 
  write_csv("./data/full_data.csv")

