domainify <- function(x) {
   # Transform a 0-1 data.frame into a domain data.frame
   # PARAMETERS:
   #   'x' 'data.frame' with columns seqnanme, start, end, 0-1.
   # RETURN:
   #   a 'data.frame' with columns seqname, start, end.

   # Sort 'x' on block name and start.
   x <- x[order(x[,1], x[,2]),]

   domains <- data.frame()
   for (block in unique(x[,1])) {
      # 'y' is the restriction of 'x' to given block.
      y <- x[x[,1] == block,]
      #starts <- which(diff(y[,4]) == 1) + 1
      #ends   <- which(diff(y[,4]) == -1)
      a <- which(diff(y[,4]) == 1)
      b <- a+1
      c <- which(diff(y[,4]) == -1)
      d <- c+1

      if (length(a) == 0 && !y[1,4]) next

      if (y[1,4]) {
         a <- c(1, a)
         b <- c(1, b)
      }
      if (y[nrow(y),4]) {
         c <- c(c, nrow(y))
         d <- c(d, nrow(y))
      }
      domains <- rbind(
         domains,
         data.frame(block=block,
           start = round((y[a,2] + y[b,2])/2, 0),
           end = round((y[c,3] + y[d,3])/2, 0))
      )
   }

   return (domains)

}
