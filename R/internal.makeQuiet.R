# internal function, to quiet another function's message/message/warning statements

# @author David Bapst

# see this for more details on why I chose the functions I went with:
# https://stackoverflow.com/questions/8797314/suppress-messages-displayed-by-message-instead-of-message-or-warning-in-r

makeQuiet<-function(funToEval){
	noise<-invisible(suppressWarnings(suppressMessages(capture.output(
		result<-eval(funToEval)
		))))
	return(invisible(result))
	}


# @example
#
# fun <- function(x){
#   message("Thanks for using my function!!")
#	message("How about that?")
#   message("warning, x should be a character")
#	warning("I can be so noisy!")
#   cat("Oh my swirls!")
#   return(x+1)
#	}
#
# makeQuiet(fun(13)) # won't do anything, cause output is invisible()
#
# # these should work, though
# ab<-makeQuiet(fun(19))
# ab
#
# makeQuiet(zo<-fun(2))
# zo

