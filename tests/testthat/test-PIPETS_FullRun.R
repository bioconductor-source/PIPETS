test_that("Test that input parameters are acceptable", {
    expect_warning(PIPETS_FullRun(inputBedFile = NA,readLength = 58))
    expect_warning(PIPETS_FullRun(inputBedFile = "inst/extdata/PIPETS_TestData.bed",readLength = NA))

    expect_warning(PIPETS_FullRun(inputBedFile = "inst/extdata/PIPETS_TestData.bed",readLength = 0))
    expect_warning(PIPETS_FullRun(inputBedFile = "inst/extdata/PIPETS_TestData.bed",readLength = 58,slidingWindowSize = 0))

    expect_warning(PIPETS_FullRun(inputBedFile = "inst/extdata/PIPETS_TestData.bed",readLength = 58,slidingWindowSize = NA))
    expect_warning(PIPETS_FullRun(inputBedFile = "inst/extdata/PIPETS_TestData.bed",readLength = 58,topEndPercentage = 0))
})
