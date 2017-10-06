# internal function, to quiet another function's print/message/warning statements

# see this for more:
# https://stackoverflow.com/questions/8797314/suppress-messages-displayed-by-print-instead-of-message-or-warning-in-r

makeQuiet<-function(funToEval){
	noise<-suppressWarnings(suppressMessages(capture.output(
		result<-eval(funToEval))))
	return(invisible(eval(funToEval)))
	}


# EXAMPLE 
#
# fun <- function(x){
#   print("Thanks for using my function!!")
#	message("How about that?")
#   if(is(x, "numeric"))
#     print("warning, x should be a character")
#	warning("I can be so noisy!")
#   return(x+1)
#	}
#
# makeQuiet(fun(13))
#
# ab<-makeQuiet(fun(19))
# ab
#
# makeQuiet(zo<-fun(2))
# zo

