###############################
#
# Non-parametric Sparse Matrix Decomposition
# 23 Dec. 2013
#
###############################
# X, Y are the first and second views of data, respectively,
# where X\in R^{n*m}, Y\in R^{n*p}
# #########################
NSMD <- function(X, Y, K){
   
   # center & scale y & x
   X <- scale(X, TRUE, apply(X,2,sd))
   Y <- scale(Y, TRUE, apply(Y,2,sd))
   
   U = NULL;             # P is the loading matrix of X
   V = NULL;             # Q is the loading matrix of Y
   
   Xr <- X
   Yr <- Y
   
   #   H = crossprod(Xr, Yr)/nrow(Xr)     
   H = cov(Xr, Yr) 
   D = NULL
   
   for(k in 1:K){
     
      u = matrix(rep(1, nrow(H)), ncol=1)
      v = matrix(rep(1, ncol(H)), ncol=1)
    
      # get the tentative values of u and v
      uv = nipals(H, u, v)
      # get the sparse values of u and v
      uv = sparsevalue(H, uv$u, uv$v)
#      uv = spa_val(H, uv$u, uv$v)
      
      H = H - uv$d * tcrossprod(uv$u, uv$v)
      
      D = c(D, uv$d)
      U <- cbind(U, uv$u)
      V <- cbind(V, uv$v)
   }
   
   return(list(D=D, U=U, V=V))
}

#
# 23 Dec. 2013
# The sparsity-inducing operation (equivalent to the sparsevalue function)
#
spa_val <- function(H, u, v){
   
   ulen = length(u)  # the length of u
   vlen = length(v)  # the length of u
   
   usort = sort(u, index.return=T)  # sort u in increasing order
   vsort = sort(v, index.return=T)  # sort v in increasing order
   
   Hs = H[usort$ix, vsort$ix]
   
   # get the coefficients d of H with u under the fix values of v
   Hv = Hs[, 1:2] %*% matrix(vsort$x[1:2],ncol=1)  # H*v
   Hv = Hv * usort$x

   
   Hv = cumsum(Hv) # u'*H*v
   
   uinds = NULL     # the indexes of u which decreasing the coefficient value of u'*H*v
   dmax = 0
   for(i in 1:ulen){
      if(Hv[i] <= dmax)
         uinds = c(uinds, i)
      else
         dmax = Hv[i]
   }
   
   # get the coefficients d of H with u under the fix values of v
   Hv = crossprod(matrix(usort$x[1:2],ncol=1), Hs[1:2, ])   # u'*H
   Hv = Hv * vsort$x
   Hv = cumsum(Hv) # u'*H*v

   vinds = NULL     # the indexes of v which decreasing the coefficient value of u'*H*v
   dmax = 0
   for(i in 1:vlen){
      if(Hv[i] <= dmax)
         vinds = c(vinds, i)
      else
         dmax = Hv[i]
   }
   
   # update u & v
   uinds = usort$ix[uinds]
   vinds = vsort$ix[vinds]
   u[uinds] = 0 
   v[vinds] = 0 
   
   # get the final values of u & v
   Hs = H[-uinds, -vinds]
   uv = nipals(Hs, matrix(1, nrow=ulen-length(uinds), ncol=1), matrix(1, nrow=vlen-length(vinds), ncol=1))
   d = drop(matrix(uv$u, nrow=1) %*% Hs %*% matrix(uv$v, ncol=1))
   
   # update u & v
   u[-uinds] = uv$u
   v[-vinds] = uv$v
   
   # get the coefficient of u'*H*v
   return(list(d=d, u=u, v=v))
}

# sparsity-inducing operation
sparsevalue <- function(H, u, v){
   
   ulen = length(u)  # the length of u
   vlen = length(v)  # the length of u
   
   usort = sort(u, index.return=T)  # sort u in increasing order
   vsort = sort(v, index.return=T)  # sort v in increasing order
   
   Hs = H[usort$ix, vsort$ix]
   
   # get the coefficients d of H with u under the fix values of v
   Hv = Hs[, 1:2] %*% matrix(vsort$x[1:2],ncol=1)  # H*v
   uinds = NULL     # the indexes of u which decreasing the coefficient value of u'*H*v
   dsum = 0; dmax = 0
   for(i in 1:ulen){
      dsum = dsum + usort$x[i] * Hv[i]  # u'*H*v
      if(dsum>dmax){
         uinds = c(uinds, i)
         dmax = dsum
      }
   }
   
   # get the coefficients d of H with u under the fix values of v
   Hv = crossprod(matrix(usort$x[1:2],ncol=1), Hs[1:2, ])   # u'*H
   vinds = NULL     # the indexes of v which decreasing the coefficient value of u'*H*v
   dsum = 0; dmax = 0
   for(i in 1:vlen){
      dsum = dsum + vsort$x[i] * Hv[i]  # u'*H*v
      if(dsum>dmax){
         vinds = c(vinds, i)
         dmax = dsum
      }
   }
   
   # update u & v
   uinds = usort$ix[uinds]
   vinds = vsort$ix[vinds]
   
   # get the final values of u & v
   Hs = H[uinds, vinds]
   uv = nipals(Hs, matrix(1, nrow=nrow(Hs), ncol=1), matrix(1, nrow=ncol(Hs), ncol=1))
   d = drop(matrix(uv$u, nrow=1) %*% Hs %*% matrix(uv$v, ncol=1))
   
   # update u & v
   u[uinds] = uv$u
   v[vinds] = uv$v
   u[-uinds] = 0 
   v[-vinds] = 0 
   
   
   # get the coefficient of u'*H*v
   return(list(d=d, u=u, v=v))
}


#############################
nipals <- function(K, u0, v0) {
   
#show(paste(nrow(K), ncol(K)))   
#show(paste( "un", nrow(u0), ncol(u0)))
#show(paste("vo", nrow(v0), ncol(v0)))

   unew=u0;
   vnew=v0;
   
   iters = 1;  erru=1; errv=1;
   while(iters<=100 && (erru>1e-4 || errv>1e-4)){
      uold = unew;
      vold = vnew;

      #update singular vector u
      unew = K %*% vnew;
      unew = unew/sqSum(unew);
      
      #update v
      vnew = crossprod(K, unew);
      vnew = vnew/sqSum(vnew);
      
      erru=max(abs(uold-unew));
      errv=max(abs(vold-vnew));
     
      iters=iters+1;
   }

   return(list(u=unew, v=vnew))
}
#####################

sqSum <- function(vec){
   a <- sum(vec^2)
   if(a==0) a <- 1
   return(sqrt(a))
}
