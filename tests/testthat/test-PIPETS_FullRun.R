test_that("Test that input parameters are acceptable", {
    expect_no_error(PIPETS_FullRun(inputData = "~/Desktop/PIPETS/inst/extdata/PIPETS_TestData.bed", readLength = 58, OutputFileID = "Test"))
    expect_warning(PIPETS_FullRun(inputData = NA,readLength = 58, OutputFileID = "Test"))
    expect_warning(PIPETS_FullRun(inputData = "inst/extdata/PIPETS_TestData.bed",readLength = NA,OutputFileID = "Test"))

    expect_warning(PIPETS_FullRun(inputData = "inst/extdata/PIPETS_TestData.bed",readLength = 0,OutputFileID = "Test"))
    expect_warning(PIPETS_FullRun(inputData = "inst/extdata/PIPETS_TestData.bed",readLength = 58,slidingWindowSize = 0,OutputFileID = "Test"))

    expect_warning(PIPETS_FullRun(inputData = "inst/extdata/PIPETS_TestData.bed",readLength = 58,slidingWindowSize = NA,OutputFileID = "Test"))
    expect_warning(PIPETS_FullRun(inputData = "inst/extdata/PIPETS_TestData.bed",readLength = 58,topEndPercentage = 0,OutputFileID = "Test"))
})







