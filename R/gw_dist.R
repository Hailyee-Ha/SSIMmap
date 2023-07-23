
gw.dist <- function(dp, sp, focus=0) {
  focus_point <- sp[focus, ]
  dists <- sqrt((dp[,1] - focus_point[1])^2 + (dp[,2] - focus_point[2])^2)
  dists
}

