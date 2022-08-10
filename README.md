# PhDNimbleFunctions
 PhD Nimble functions for simple access and write up.

 Currently I am accessing these via:

path <- "where this is downloaded"
files <- dir(path)
files <- grep('.R', files, value = TRUE, ignore.case = TRUE)
lapply(files, FUN = function(x){source(paste0(path, "/", x))})
source(paste0(path, "/", files[1]))

I know I should wrap this up as a package but for my purposes this is easy enough.