##线性规划问题
lp_minusb <- function(cvec,bvec,vitx,vity,vitz,maximum){
  Amat <- rbind(vitx, vity, vitz)
  const.dir = rep("<=", length(bvec))
  #rbind将自变量纵向组合成一个大矩阵(系数矩阵)
  maxiter = 1000#最大迭代次数
  zero = 0
  ####正式代码
  nVar <- length(cvec)#cvec：价值系数，个数对应变量个数
  nCon <- length(bvec)#bvec：约束条件
  if (!all.equal(dim(Amat), c(nCon, nVar)) == TRUE) #m(约束) x n(变量) 矩阵
  {
    stop(paste("Matrix A must have as many rows as constraints (=elements", "of vector b) and as many columns as variables (=elements of vector c).\n"))
  }
  if (length(const.dir) != nCon)#约束符号个数应该和右端项个数相同
  {
    stop(paste("'const.dir' must have as the elements as constraints",  "(=elements of vector b).\n"))
  }
  if (is.null(names(cvec))) #是否给变量命名过名字
  {
    clab <- as.character(1:nVar)
  }else {
    clab <- names(cvec)
  }
  if (is.null(names(bvec)))#是否给约束命名过名字
  {
    blab <- as.character(1:nCon)
  }else {
    blab <- names(bvec)
  }
  const.dir2 <- rep(0, nCon)
  const.dir2[const.dir == ">=" | const.dir == ">"] <- -1  #都写成<=形式：如果不是<=形式，后面处理约束条件时，两边乘以-1
  const.dir2[const.dir == "<=" | const.dir == "<"] <- 1   #如果已经是<=形式，就不用两边乘以-1
  ###引入松弛变量，写成标准形式
  iter1 <- 0
  for (i in 1:nCon) # nCon约束条件个数x，相应地引入x个松弛变量
  {
    clab <- c(clab, paste("S", blab[i]))#引入松弛变量；paste(a,b)将a和b变量连成一个字符串。x个约束条件对应x个松弛变量
  }
  cvec2 <- c(cvec, rep(0, nCon))#变量的价值系数中加入nCon个松弛变量的系数（0）
  Tab <- rbind(cbind(Amat * const.dir2,  diag(1, nCon, nCon), bvec * const.dir2), c(cvec2 * (-1)^maximum, 0))   #Tab为单纯形表格式
  #求最大值，当所有检验数δ<0时，达到最优；否则选取最大δ对应的变量入基
  #求最小值，当所有检验数δ>0时，达到最优；否则在δ<0中，选取最小的δ对应的变量入基
  # diag(1, nCon, nCon):通过构建对角矩阵，得到松弛变量系数矩阵，并和已有变量的系数矩阵横向合并；
  #在R中，TRUE表示1，FALSE表示0；都统一成最小值求解（若原意是求最大值，在此改成求最小值）
  #class(Tab)   matrix
  rownames(Tab) <- c(blab, "-z")#rownames为矩阵行变量取名字 即-z= f(xn) 
  colnames(Tab) <- c(clab, "b")#colnames为矩阵列变量取名字
  
  ###开始迭代
  while (min(Tab[nCon + 1, 1:(nVar + nCon)]) < -zero & iter1 < maxiter) #判断最小检验数δ是否小于0
  {
    iter1 <- iter1 + 1
    pcolumn <- which.min(Tab[nCon + 1,1:(nVar + nCon)])  #选择最小的δ的所在的列号，确定入基变量（该变量进基能使函数值降低最快）
    #cat("输出第",iter1,"次迭代最小δ的列号\n")
    #print( pcolumn)
    rθ <- Tab[1:nCon,'b'] / Tab[1 : nCon,pcolumn] #计算θ
    #cat("输出第",iter1,"次迭代实际的θ1\n")#cat：输出对象
    #print(rθ)
    rθ[Tab[1:nCon, pcolumn] <= 0] <- max(rθ,  na.rm = TRUE) + 1#若有θi<0，就将max(θ)+1赋给θi;可将第3次迭代，实际的θ1和经过处理的θ1‘进行对比，加强理解
    #cat("输出第",iter1,"次迭代的，经过处理的θ1‘\n")
    # print(rθ)
    prow <- which.min(rθ)  #按照θ规则，返回最小θ的下标（行号）,确定出基变量。从而由入基变量的列号和出基变量的行号就能确定主元
    rownames(Tab)[prow] <- colnames(Tab)[pcolumn]  #将入基变量名写入到单纯形表左侧（中间单纯形表左侧变量名和标准解会不一样，但最终单纯形表与标准解一致）
    Tab[prow, ] <- Tab[prow, ] / Tab[prow,pcolumn] #以主元为分母，将出基变量所在行的系数和约束项b均除以主元
    #for中代码，是对其它行（约束条件行和δ行）做初等行变换，请补充Tab[i,]初等行变换的代码，参考例1的单纯形表。注意：需要考虑i是否等于prow
    for(i in 1:(nCon + 1)){
      if(i != prow){
        a = Tab[i,pcolumn] / Tab[prow,pcolumn]
        Tab[i,] <- Tab[i,] - a * Tab[prow,]
      }
    }
    #注意：此处的初等行变换最好对主元行做处理，得到其它行（单纯形表-3-p22）
    cat("输出第",iter1,"次迭代的单纯形表\n")
    print(Tab)
  }	
  if(maximum)#因为经前面处理后，是求最小值，所以如果求最大值，需要取相反数 
  {
    Tab[nCon + 1, (nVar + nCon+1)]= -Tab[nCon + 1, (nVar + nCon+1)]
  }
  opt=	Tab[nCon + 1, (nVar + nCon+1)]  
  #输出结果	   
  cat("输出最优单纯形表和最优函数值\n")
  print(Tab)#默认换行
  cat("\n",opt,"\n")
  #返回函数值
  result<-list(OptSimplex=Tab,OptSoultion=-opt)
  return(result)
}