# Copyright 2024 Josef Hackl and Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


library(qqman)

args <- commandArgs(trailingOnly = TRUE)
data <- read.table(args[1], header=TRUE)
png(args[2], width=800, height=600, res=96)
manhattan(data, p="B1", col=c("blue", "orange"), logp=FALSE, genomewideline=as.numeric(args[3]), suggestiveline=as.numeric(args[4]), main="Lithuanian", ylab=expression(beta^{(1)} ~ "scores"), ylim=c(0,50))
dev.off()
