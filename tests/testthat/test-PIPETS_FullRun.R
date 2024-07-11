test_that("Test that input parameters are acceptable", {
    expect_no_error(PIPETS_FullRun(inputData = "PIPETS_TestData.bed", readScore = 42, OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = NA,readScore = 42, OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readScore = NA,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readScore = 42,OutputFileID = "Test",OutputFileDir = "~/Dersktop/"))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readScore = 0,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readScore = 42,slidingWindowSize = 0,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readScore = 42,slidingWindowSize = NA,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
})







