test_that("pairwiseKS works", {
	data(simRunExample)
	pdfList <- list(
		Brownian = results[[1]]$particleDataFrame,
		Bounded = resultsBound[[1]]$particleDataFrame)
	res<-pairwiseKS(particleDataFrameList = pdfList)
	})

test_that("pairwiseESS works", {
	data(simRunExample)
	res<-pairwiseESS(list(results[[1]]$particleDataFrame,results[[2]]$particleDataFrame))
	expect_message(
		res<-pairwiseESS(results)
		)
	expect_warning(
		res<-pairwiseESS(results[[1]]$particleDataFrame)
		)
	}
	)
