context("Testing global outliers removed")


test_that("NAs are removed when Global Outliers are present", {
  # loads a WTA dataset from inst after wrangling it
  data_path <- sub("/inst", "", gsub("\\\\", "/", gsub('/\\\\inst', '\\\\inst', base::file.path(
    devtools::inst(name="dspNgs"),
    "\\inst\\testData\\Spatial_data\\"
  ))))

  norm_counts <- read.delim(dir(data_path, full.names=T, pattern="_NegNorm"),
                            sep = "\t", header = TRUE)
  probe_notes <- read.delim(dir(data_path, full.names=T, pattern="_TargetProperties.txt$"), sep = "\t",
                            header = TRUE, colClasses = c(Pooling='character'))
  probe_notes_mod <- droplevels(probe_notes[!is.na(probe_notes[["GlobalOutliers"]]),])
  # tests
  expect_true("GlobalOutliers" %in% colnames(probe_notes))
  expect_true(nrow(norm_counts)==18293)
  expect_true(nrow(probe_notes)==18372)
  expect_true(nrow(probe_notes_mod)==18293)
  expect_true("TargetName" %in% colnames(norm_counts))
  expect_true("TargetName" %in% colnames(probe_notes_mod))
  expect_true(class(norm_counts$TargetName)=="factor")
  expect_true(class(probe_notes_mod$TargetName)=="factor")
  expect_true(all(as.character(probe_notes_mod$TargetName) %in% as.character(norm_counts$TargetName)))
  expect_true(all(as.character(norm_counts$TargetName) %in% as.character(probe_notes_mod$TargetName)))

  # Housekeeping
  rm(data_path)
  rm(norm_counts)
  rm(probe_notes)
})


