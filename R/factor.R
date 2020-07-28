###Functions of Factor Models

# This approach is for the vector-valued estimation WITHOUT NaNs.
#'@name mfmda.nona.vec
#'@rdname mfmda.nona.vec
#'@aliases mfmda.nona.vec
#'@export
#'@param Yc Time Series data for a matrix(dimensions n*p*q), no NA input allowed
#'@param hzero Pre-scribed parameter
#'@return The sample version of M matrix
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
#'@examples
#' A <- 1:180
#' dim(A) <- c(3,3,20)
#' M <- mfmda.nona.vec(A,2)
mfmda.nona.vec <- function(Yc,hzero){
  dimYc = dim(Yc)
  n = dimYc[1]
  p = dimYc[2]
  q = dimYc[3]
  Yc.matrix <- t(matrix(Yc,nrow=n))
  Mhat = matrix(0,p*q,p*q)
  for (h in 1:hzero){
    autocov.all <- ( Yc.matrix[,1:(n-h)] %*% t(Yc.matrix[,(h+1):n]) )/ (n-h)
    Mhat <- Mhat + autocov.all %*% t(autocov.all)
  }
  Mhat
}

# This approach is for the vector-valued estimation with NaNs.
#'@name mfmda.na.vec
#'@rdname mfmda.na.vec
#'@aliases mfmda.na.vec
#'@export
#'@param Yc Time Series data for a matrix(dimensions n*p*q), allowing NA input
#'@param hzero Pre-scribed parameter h 
#'@return The sample version of M matrix 
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
mfmda.na.vec <- function(Yc,hzero){
  dimYc = dim(Yc)
  n = dimYc[1]
  p = dimYc[2]
  q = dimYc[3]
  Yc.matrix <- t(matrix(Yc,nrow=n))
  Mhat = matrix(0,p*q,p*q)
  for (h in 1:hzero){
    Omegah <- matrix(0,p*q,p*q)
    count <- 0
    for (tt in 1:(n-h)){
      if( (sum(is.na(Yc.matrix[,tt]))==0) && (sum(is.na(Yc.matrix[,tt+h]))==0) ){
        Omegah = Omegah + Yc.matrix[,tt] %*% t(Yc.matrix[,tt+h])
        count <- count+1
      }                    
    }
    if (count>0){
      Omegah = Omegah/count
    }
    Mhat <- Mhat + Omegah %*% t(Omegah)
  }
  Mhat
}

# The input data do not have zeros. The estimation approach is noniterative. 
#'@name mfmda.nona.noniter
#'@rdname mfmda.nona.noniter
#'@aliases mfmda.nona.noniter
#'@export
#'@param Yc Time Series data for a matrix(dimensions n*p*q), no NA input allowed
#'@param hzero Pre-scribed parameter
#'@return The sample version of M matrix
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
mfmda.nona.noniter <- function(Yc,hzero){
  dimYc = dim(Yc)
  n = dimYc[1]
  p = dimYc[2]
  q = dimYc[3]
  Yc.matrix <- t(matrix(Yc,nrow=n))
  Mhat = matrix(0,p,p)
  for (h in 1:hzero){
    autocov.all <- Yc.matrix[,1:(n-h)] %*% t(Yc.matrix[,(h+1):n])
    Gammayh = matrix(0,p,p)
    for (ii in 1:q){
      for (jj in 1:q){
        Omegaijh = matrix(0, p, p)
        Omegaijh = autocov.all[((ii-1)*p+1):(ii*p),((jj-1)*p+1):(jj*p)]
        Omegaijh = Omegaijh/(n - h)
        Gammayh = Gammayh + Omegaijh %*% t(Omegaijh);
      }
    }
    Mhat = Mhat + Gammayh   
  }
  Mhat
}

# The input data do not have zeros. The estimation approach is iterative.
#'@name mfmda.nona.iter
#'@rdname mfmda.nona.iter
#'@aliases mfmda.nona.iter
#'@export
#'@param Yc Time Series data for a matrix(dimensions n*p*q), no NA input allowed
#'@param hzero Pre-scribed parameter
#'@return The sample version of M matrix
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
mfmda.nona.iter <- function(Yc,hzero){
  dimYc = dim(Yc)
  n = dimYc[1]
  p = dimYc[2]
  q = dimYc[3]
  Mhat = matrix(0,p,p)
  for (h in 1:hzero){
    Gammayh = matrix(0,p,p)
    for (ii in 1:q){
      for (jj in 1:q){
        Omegaijh = matrix(0, p, p)
        for (tt in 1:(n-h)){
          Omegaijh = Omegaijh + Yc[tt,,ii] %*% t(Yc[tt+h,,jj])
        }
        Omegaijh = Omegaijh/(n - h)
        Gammayh = Gammayh + Omegaijh %*% t(Omegaijh);
      }
    }
    Mhat = Mhat + Gammayh   
  }
  Mhat
}

# The input data could have NaNs. The estimation approach is iterative.
#'@name mfmda.na.iter
#'@rdname mfmda.na.iter
#'@aliases mfmda.na.iter
#'@export
#'@param Yc Time Series data for a matrix allowing NaNs
#'@param hzero Pre-determined parameter
#'@return The sample version of M matrix
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
mfmda.na.iter <- function(Yc,hzero){
  dimYc = dim(Yc)
  n = dimYc[1]
  p = dimYc[2]
  q = dimYc[3]
  Mhat = matrix(0,p,p)
  for (h in 1:hzero){
    Gammayh = matrix(0,p,p)
    for (ii in 1:q){
      for (jj in 1:q){
        Omegaijh = matrix(0, p, p)
        ## deal with missing data
        count <- 0
        for (tt in 1:(n-h)){
          if( (sum(is.na(Yc[tt,,ii]))==0) && (sum(is.na(Yc[tt+h,,jj]))==0) ){
            Omegaijh = Omegaijh + Yc[tt,,ii] %*% t(Yc[tt+h,,jj])
            count <- count+1
          }                    
        }
        if (count>0){
          Omegaijh = Omegaijh/count
        }
        Gammayh = Gammayh + Omegaijh %*% t(Omegaijh);
      }
    }
    Mhat = Mhat + Gammayh   
  }
  Mhat
}

## This is a wrapper for all approaches
#'@name mfmda
#'@rdname mfmda
#'@aliases mfmda
#'@export
#'@param Yc Time Series data for a matrix
#'@param approach Select estimation approaches, 1 for noniterative approach with no NaNs, 2 for iterative approach with NaNs, 3 for iterative approach allowing NaNs.
#'@param hzero Pre-determined parameter
#'@param iscentering The data is subtracted by its mean value
#'@return The sample version of M matrix
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
#'@examples
#' A <- 1:180
#' dim(A) <- c(3,3,20)
#' M <- mfmda(A,"3",1,0)
mfmda <- function(Yt,approach="3",hzero=1,iscentering=1){
  ## Dimensions
  dims <- dim(Yt)
  n <- dims[1] ## number of time points (observations)
  p <- dims[2] ## number of rows of each matrix-valued observation
  q <- dims[3] ## number of columns of each matrix-valued observation
  ## Centering the data
  if (iscentering){
    ## subtract the mean
    Yt.mean <- apply(Yt,2:3,mean,na.rm=TRUE)
    Yc <- Yt - array(rep(Yt.mean,rep(n,p*q)),c(n,p,q))
  } else{
    Yc <- Yt
  }
  rm(Yt)
  if (approach == "1"){
    Mhat <- mfmda.nona.noniter(Yc,hzero)
  }
  if (approach == "2"){
    Mhat <- mfmda.nona.iter(Yc,hzero)
  }
  if (approach == "3"){
    Mhat <- mfmda.na.iter(Yc,hzero)
  }
  Mhat
}

# Compute the estimated number of factors and the corresponding eigen-space
#'@name mfmda.estqk
#'@rdname mfmda.estqk
#'@aliases mfmda.estqk
#'@export
#'@param Mhat The estimated value for matrix M
#'@param inputk The pre-determined number of dimension of factor matrix
#'@return The estimated number of factors to use, the corresponding estimated Q matrix, the eigenvalue,  the estimated Q matrix with requested number of factors 
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
#'@examples
#' A <- 1:180
#' dim(A) <- c(3,3,20)
#' M <- mfmda(A,"3",1,0)
#' inputk <- 3
#' eig.ans <- mfmda.estqk(M,inputk)
#' khat <- eig.ans$estk
#' Qhat <- eig.ans$Qhatestk
#' eigval <- eig.ans$eigval
#' Qhatinputk <- eig.ans$Qhatinputk
mfmda.estqk <- function(Mhat,inputk=1){
  dimMhat = dim(Mhat)
  p = dimMhat[1]
  cc = round(p/2)
  eig <- eigen(Mhat)
  eigval = eig$values
  eigvec = eig$vectors
  # Estimated number of factors
  eigratio = eigval[2:(cc+1)]/eigval[1:cc]
  estk = which.min(eigratio)
  # Estimated eigen-vectors
  Qhatestk = eigvec[,1:estk]
  # Return the requested number of eigen-vectors
  Qhatinputk = eigvec[,1:inputk]
  list("estk" = estk, "Qhatestk" = Qhatestk,"eigval" = eigval,"Qhatinputk"=Qhatinputk)
}

# The main estimation function
#'@name matrix_factor
#'@rdname matrix_factor
#'@aliases matrix_factor
#'@export
#'@param Yt Time Series data for a matrix
#'@param inputk1 The pre-determined row dimension of the factor matrix
#'@param inputk2 The pre-determined column dimension of the factor matrix
#'@param iscentering The data is subtracted by its mean value
#'@param hzero Pre-determined parameter

#'@return a list containing the following:\describe{
#'\item{\code{eigval1}}{estimated row dimension of the factor matrix}
#'\item{\code{eigval2}}{estimated column dimension of the factor matrix}
#'\item{\code{loading1}}{estimated left loading matrix}
#'\item{\code{loading2}}{estimated right loading matrix}
#'\item{\code{Ft}}{Estimated factor matrix with pre-determined number of dimensions}
#'\item{\code{Ft.all}}{Sum of Ft}
#'\item{\code{Et}}{The estimated residual, by subtracting estimated signal term from the data}
#'}

#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
#'@examples
#' A <- 1:180
#' dim(A) <- c(3,3,20)
#' out = matrix_factor(A,3,3)
#' eig1 = out$eigval1
#' loading1 = out$loading1
#' Ft = out$Ft.all
matrix_factor=function(Yt,inputk1,inputk2,iscentering=1,hzero=1){
  dims <- dim(Yt)
  n <- dims[1] ## number of time points (observations)
  p <- dims[2] ## number of rows of each matrix-valued observation
  q <- dims[3] ## number of columns of each matrix-valued observation
  ## Centering and standardizing the data
  if(iscentering==1){
    Yt.mean <- apply(Yt,2:3,mean,na.rm=TRUE)
    Yt.sd <- apply(Yt,2:3,sd,na.rm=TRUE)
    Yt.sd[Yt.sd==0]=1   ## if sd=0, do nothing. series will be constant zero
    Yc <- (Yt - array(rep(Yt.mean,rep(n,p*q)),c(n,p,q)))/array(rep(Yt.sd,rep(n,p*q)),c(n,p,q))
  }
  else{ 
    Yc <- Yt
    ## ------------- ESTIMATION FOR p dimension -------------------------- 
    Mhat1 <- mfmda(Yc,"3",hzero,iscentering)
    eig1.ans <- mfmda.estqk(Mhat1,inputk1)
    k1hat <- eig1.ans$estk
    Q1hat <- eig1.ans$Qhatestk
    eigval1 <- eig1.ans$eigval
    Q1hatinputk <- eig1.ans$Qhatinputk
    Q1hatinputkrot.ans <- varimax(Q1hatinputk)
    Q1hatinputkrot <- Q1hatinputk %*% Q1hatinputkrot.ans$rotmat
    for(i in 1:inputk1){
      Q1hatinputkrot[,i]=Q1hatinputkrot[,i]*sign(sum(Q1hatinputkrot[,i]))
    }
    ## ------------- ESTIMATION FOR q dimension -------------------------- 
    tYc = array(0,c(n,q,p))
    for(nk in 1:n){
      tYc[nk,,] = t(Yc[nk,,])
    }
    Mhat2 <- mfmda(tYc,"3",hzero,iscentering)
    eig2.ans <- mfmda.estqk(Mhat2,inputk2)
    k2hat <- eig2.ans$estk
    Q2hat <- eig2.ans$Qhatestk
    eigval2 <- eig2.ans$eigval
    Q2hatinputk <- eig2.ans$Qhatinputk
    Q2hatinputkrot.ans <- varimax(Q2hatinputk)
    Q2hatinputkrot <- Q2hatinputk %*% Q2hatinputkrot.ans$rotmat
    for(i in 1:inputk1){
      Q2hatinputkrot[,i]=Q2hatinputkrot[,i]*sign(sum(Q2hatinputkrot[,i]))
    }
    Ft.inputk.rot <- array(0,c(n,inputk1,inputk2))
    for(tk in 1:n){
      Ft.inputk.rot[tk,,] <- t(Q1hatinputkrot) %*% Yc[tk,,] %*% Q2hatinputkrot
    }
    Ft.all=apply(Ft.inputk.rot,c(2,3),sum)
    ## residual
    Et.inputk.rot <- array(0,c(n,p,q))
    for(tk in 1:n){
      Et.inputk.rot[tk,,] <- Yc[tk,,] - Q1hatinputkrot %*% t(Q1hatinputkrot) %*% Yc[tk,,] %*% Q2hatinputkrot %*% t(Q2hatinputkrot)
    }
    return(list(eigval1=eigval1/p/q,eigval2=eigval2/p/q,
                loading1=Q1hatinputkrot,loading2=Q2hatinputkrot,
                Ft=Ft.inputk.rot,Ft.all=Ft.all,Et=Et.inputk.rot))
}
  
# The main estimation function, vector version
#'@name vector_factor
#'@rdname vector_factor
#'@aliases vector_factor
#'@export
#'@param Yt Time Series data for a matrix
#'@param inputk.vec The pre-determined dimensions of the factor matrix in vector
#'@param iscentering The data is subtracted by its mean value
#'@param hzero Pre-determined parameter
#'@return a list containing the following:\describe{
#'\item{\code{eigval1}}{estimated dimensions of the factor matrix}
#'\item{\code{loading}}{estimated loading matrices}
#'\item{\code{Ft}}{Estimated factor matrix with pre-determined number of dimensions}
#'\item{\code{Ft.all}}{Sum of Ft}
#'\item{\code{Et}}{The estimated random term, by subtracting estimated signal term from the data}
#'}
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
#'@examples
#' A <- 1:180
#' dim(A) <- c(3,3,20)
#' M <- mfmda(A,"3",1,0)
#' eig.ans <- vector_factor(M,3,0,1)
#' khat <- eig.ans$estk
#' Qhat <- eig.ans$Qhatestk
#' eigval <- eig.ans$eigval
#' Q1hatinputk <- eig.ans$Qhatinputk

vector_factor=function(Yt,inputk.vec,iscentering=1,hzero=1){
  Mhat.vec.zero <- mfmda.nona.vec(Yc,hzero)
  eig.vec.zero.ans <- mfmda.estqk(Mhat.vec.zero,inputk.vec)
  khat.vec.zero <- eig.vec.zero.ans$estk
  Qhat.vec.zero <- eig.vec.zero.ans$Qhatestk
  eigval.vec.zero <- eig.vec.zero.ans$eigval
  Qhatinputk.vec.zero <- eig.vec.zero.ans$Qhatinputk
  Yc.matrix <- t(matrix(Yc,nrow=n))  
  ## Extract latent factors: Ft of dimension ( inputk \times n)
  Ft.vec.inputk <- t(Qhatinputk.vec.zero) %*% Yc.matrix
  Ft.all=apply(Ft.vec.inputk,2,sum)
  ## Error term
  Et.vec.inputk <- Yc.matrix - Qhatinputk.vec.zero %*% Ft.vec.inputk
  return(list(eigval=eigval.vec.zero/p/q,
              loading=Qhatinputk.vec.zero,
              Ft=Ft.vec.inputk,Ft.all=Ft.all,Et=Et.vec.inputk))
}

# Get the group of loadings
#'@name grouping.loading
#'@rdname grouping.loading
#'@aliases grouping.loading
#'@export
#'@param loading The estimated loading matrix
#'@param ncluster The number of clusters to use, usually the dimension of the factor matrix
#'@param rowname The name of the rows
#'@param plot plot the clustering graph, defacult True
#'@return Loading matrix after grouping

#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
grouping.loading=function(loading,ncluster,rowname,plot=T){
  ddd <- dist(loading, method = "euclidean") # distance matrix
  fit <- hclust(ddd, method="ward.D") 
  if(plot==T){
    par(mfrow=c(1,1),mai=0.5*c(1,1,1,1))
    plot(fit, main='clustering') # display dendogram
    # draw dendogram with red borders around the k clusters
    rect.hclust(fit, ncluster, border="red")
  }
  group=cutree(fit,ncluster)
  group_sort=sort(group,index.return=TRUE)
  rownames(loading)=rowname
  new.loading=loading[group_sort$ix,]
  return(list(new.loading=new.loading))
}

# Get the adjacency matrix for plotting
#'@name dynamic_A
#'@rdname dynamic_A
#'@aliases dynamic_A
#'@export
#'@param x input the original estimated loading matrix
#'@param factor_count  The number of factors to use
#'@param simple.flag if True, only eliminate the entries below threshold and make all row sums to be 1; if False, the approach further eliminates the entries of the rows that are very close to threshold value and only leaves the maximum entry of each row
#'@param threshold A parameter to eliminate very small entries of the loading matrix
#'@return The new loading matrix with all rows sum to be 1
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
dynamic_A <- function(x,factor_count,simple.flag,threshold){
  if(simple.flag==F){
    m <- round(10*x)
    row_num <- dim(m)[1]
    for(i in 1:row_num){
      gap_old <- m[i,]
      if(sum(gap_old==rep(0,factor_count))==factor_count){
        gap_new <- gap_old
      }else{
        gap_old[which(gap_old<(max(gap_old)-threshold))] <- 0
        gap_old <- gap_old/sum(gap_old)
        gap_temp <- gap_old
        gap_temp[which(gap_temp<(max(gap_temp)-threshold))] <- 0
        gap_new <- gap_temp/sum(gap_temp)
        while(sum(gap_new==gap_old)!=factor_count){
          gap_old <- gap_new
          gap_old[which(gap_old<(max(gap_old)-threshold))] <- 0
          gap_old <- gap_old/sum(gap_old)
          gap_temp <- gap_old
          gap_temp[which(gap_temp<(max(gap_temp)-threshold))] <- 0
          gap_new <- gap_temp/sum(gap_temp)
        }
        m[i,] <- gap_new
      }
    }
  }
  if(simple.flag==T){
    for(i in 1:dim(x)[1]){
      x[i,which(x[i,]<=threshold)] <- 0
      x[i,] <- x[i,]/sum(x[i,])
    }
    m <- x
  }
  return(m)
}

# Plot the network graph
#'@name PlotNetwork_AB
#'@rdname PlotNetwork_AB
#'@aliases PlotNetwork_AB
#'@export
#'@param Ft The estimated factor matrix
#'@param iterated_A The left loading matrix
#'@param iterated_B The right loading matrix
#'@param labels The row labels
#'@return Plot the network graph
#'#'#'@seealso \code{\link{MAR1.projection}}            # Not decided
PlotNetwork_AB <- function(Ft,iterated_A,iterated_B=iterated_A,labels=use2){
  pre_F <- Ft
  country_label <- labels
  F_row <- dim(pre_F)[1]
  F_col <- dim(pre_F)[2]
  diag_F <- pre_F[col(pre_F)==row(pre_F)]
  temp_diag <- diag_F/max(diag_F)
  offdiag_F <- pre_F[col(pre_F)!=row(pre_F)]
  temp_offdiag <- offdiag_F/max(offdiag_F)
  averaged_F <- matrix(rep(0,F_row*F_col),nrow=F_row,ncol=F_col)
  averaged_F[col(averaged_F)==row(averaged_F)] <- temp_diag
  averaged_F[col(averaged_F)!=row(averaged_F)] <- temp_offdiag
  averaged_F <- abs(averaged_F)
  factor_count <- dim(averaged_F)[1]
  country_count1 <- length(country_label)
  temp <- matrix(rep(0,(factor_count*2+country_count1*2)^2),
                 ncol=factor_count*2+country_count1*2,nrow=factor_count*2+country_count1*2)
  temp[1:factor_count,(factor_count+1):(factor_count*2)] <- t(abs(averaged_F))
  tempFF <- temp[1:(factor_count*2),1:(factor_count*2)] 
  tempFF[col(tempFF)==row(tempFF)]<- rep(abs(temp_diag),2)
  temp[1:(factor_count*2),1:(factor_count*2)] <- tempFF
  temp[(factor_count*2+1):(factor_count*2+country_count1),1:factor_count] <- iterated_A
  temp[(factor_count*2+country_count1+1):(dim(temp)[1]),(factor_count+1):(factor_count*2)] <- iterated_B
  omit_index_A <- which(apply(iterated_A==0,1,sum)==factor_count)
  omit_index_B <- which(apply(iterated_B==0,1,sum)==factor_count)
  if(length(omit_index_A)!=0){
    final_label_A <- country_label[-omit_index_A]
    Censor_A <- omit_index_A+factor_count*2
  }else{
    final_label_A <- country_label
    Censor_A <- 0
  }
  if(length(omit_index_B)!=0){
    final_label_B <- country_label[-omit_index_B]
    Censor_B <- omit_index_B+factor_count*2+country_count1
  }else{
    final_label_B <- country_label
    Censor_B <- 0
  }
  if(sum(c(Censor_A,Censor_B)==0)==2){
    final_matrix <- temp
  }else{
    final_matrix <- temp[-c(Censor_A,Censor_B),-c(Censor_A,Censor_B)]
  }
  edge_vertex_colorset <- brewer.pal(12,'Paired')
  cluster_colorset <- brewer.pal(12,"Set3")
  
  vertex_color <- c(rep(edge_vertex_colorset[1:factor_count],2),
                    rep(edge_vertex_colorset[factor_count+1],dim(final_matrix)[1]-factor_count*2))
  #cluster_color <- cluster_colorset[1:factor_count]
  # mark.groups <- list()
  # for(i in 1:factor_count){
  #   temp_cluster <- which(final_matrix[,i]!=0)
  #   mark.groups[[i]] <- c(i,temp_cluster[-(1:factor_count)])
  # }
  net1=graph.adjacency(final_matrix,mode="directed",weighted=TRUE,diag=FALSE)
  edge.width <- c(15*E(net1)$weight[1:(factor_count^2)],E(net1)$weight[(factor_count^2+1):length(E(net1)$weight)])
  #E(net)$weight 
  edge.lty <- c(rep(1,factor_count^2),rep(3,length(E(net1)$weight)-(factor_count^2)))
  V(net1)$color=vertex_color
  V(net1)$size=c(15*rep(averaged_F[col(averaged_F)==row(averaged_F)],2),rep(1,dim(final_matrix)[1]-factor_count*2))
  vertex.label <- c(paste(1:factor_count),paste(-(1:factor_count)),final_label_A,final_label_B)
  temp_edge.color <- rep(edge_vertex_colorset[1:factor_count],each=factor_count)
  edge.color <- c(temp_edge.color,rep(edge_vertex_colorset[factor_count+1],length(E(net1)$weight)-(factor_count^2)))
  ### set location of vertex #################
  rr <- factor_count
  psc1 <- length(final_label_A)
  psc2 <- length(final_label_B)
  ll <- cbind(c(rep(3.5,rr),rep(6.5,rr),
                rep(1,psc1), rep(9,psc2)),
              c(2+(rr:1)*6/(rr+1), 2+(rr:1)*6/(rr+1),
                (psc1:1)*10/(psc1+1), (psc2:1)*10/(psc2+1)))
  #################################
  edge.curved <- rep(FALSE,length(E(net1)))
  plot.igraph(net1,
              vertex.label.cex=0.65,
              vertex.label=vertex.label,
              #layout=layout.fruchterman.reingold,
              #layout=layout_as_star,
              layout=ll,
              #vertex.frame.color=NA,
              margin=c(0.001,0.001,0.001,0.001),
              edge.width=edge.width,
              edge.curved=edge.curved,
              edge.arrow.size=0,
              edge.color=edge.color,
              edge.lty=edge.lty)  
  #mark.groups=mark.groups,
  #mark.col=cluster_color,mark.border=NA,
}
