x <- c(
  "Annual Mean Temperature",
  "Mean Diurnal Range (Mean of monthly (max temp - min temp))",
  "Isothermality (BIO2/BIO7) (×100)",
  "Temperature Seasonality (standard deviation ×100)" ,
  "Max Temperature of Warmest Month",
  "Min Temperature of Coldest Month",
  "Temperature Annual Range (BIO5-BIO6)",
  "Mean Temperature of Wettest Quarter",
  "Mean Temperature of Driest Quarter",
  "Mean Temperature of Warmest Quarter",
  "Mean Temperature of Coldest Quarter",
  "Annual Precipitation",
  "Precipitation of Wettest Month",
  "Precipitation of Driest Month",
  "Precipitation Seasonality (Coefficient of Variation)",
  "Precipitation of Wettest Quarter",
  "Precipitation of Driest Quarter",
  "Precipitation of Warmest Quarter",
  "Precipitation of Coldest Quarter")

test_that("correct errors are being thrown", {
  expect_error(did_you_mean(x, "XXX"),
               regexp = "No match for .* in x")
  expect_error(did_you_mean(x, "Annual Mean Temperature"),
               regexp = "y found in x\\..*")
  expect_error(did_you_mean(x, "Annual Temperature"),
               regexp = "No match for .* - Did you mean .*")
  expect_error(did_you_mean(x, "Mean Temperature"),
               regexp = "No match for .* - Did you mean one of:\\n.*\\n.*\\.*\\n.*\\n.*\\n.*")
  expect_error(did_you_mean(x, "Temperature"),
               regexp = "\\d+? equally good matches for .* - Please specify.")
  expect_error(did_you_mean(x, 42),
               regexp = "Both x and y have to be of type character.")
})
