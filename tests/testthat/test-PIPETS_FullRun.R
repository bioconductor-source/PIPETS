test_that("Test that input parameters are acceptable", {
    expect_no_error(PIPETS_FullRun(inputData = "PIPETS_TestData.bed", readLength = 58, OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    
    expect_warning(PIPETS_FullRun(inputData = NA,readLength = 58, OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readLength = NA,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readLength = 58,OutputFileID = "Test",OutputFileDir = "~/Dersktop/"))

    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readLength = 0,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readLength = 58,slidingWindowSize = 0,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))

    expect_warning(PIPETS_FullRun(inputData = "PIPETS_TestData.bed",readLength = 58,slidingWindowSize = NA,OutputFileID = "Test",OutputFileDir = as.character(tempdir())))
})







