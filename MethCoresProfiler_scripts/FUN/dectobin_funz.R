 # MethCoresProfiler is a R-script that provides a simple method to trace and track 
 # cores shared by epiallele families in complex populations. 
 # Copyright (C) 2020 author: Antonio Pezone 
 # email: antoniopezone@gmail.com; antonio.pezone@unina.it

 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # any later version.

 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.

 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <https://www.gnu.org/licenses/>.

dectobin <- function(y) {
    # find the binary sequence corresponding to the decimal number 'y'
    stopifnot(length(y) == 1, mode(y) == 'numeric')
    q1 <- (y / 2) %/% 1
    r <- y - q1 * 2
    res = c(r)
    while (q1 >= 1) {
      q2 <- (q1 / 2) %/% 1
      r <- q1 - q2 * 2
      q1 <- q2
      res = c(r, res)
    }
    return(res)
  }