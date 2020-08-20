test_that("pairwiseKS works", {
	data(simRunExample)
	
	pdfList <- list(
		Brownian = resultsBMExample[[1]]$particleDataFrame,
		Bounded = resultsBoundExample[[1]]$particleDataFrame
		)
	
	res <- pairwiseKS(particleDataFrameList = pdfList)
	
	expect_is(res, "list")
	
	})


test_that("pairwiseESS works", {
	
	data(simRunExample)
	
	res<-pairwiseESS(list(
		resultsBMExample[[1]]$particleDataFrame,
		resultsBMExample[[2]]$particleDataFrame
		))
		
	expect_is(res, "matrix")
	
	expect_message(
		res<-pairwiseESS(resultsBMExample)
		)
	
	expect_warning(
		res<-pairwiseESS(resultsBMExample[[1]]$particleDataFrame)
		)
	
	}
	)
