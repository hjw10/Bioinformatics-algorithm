##背包问题
 KP <- function (Num,Vol, Weight,Value)
 {
  #step 1:初始化矩阵OptValue：OptValue[i][j]表示前i个物品装入容量为j的背包中获得的最大价值； 结果见ppt 13
  Nrow=Num+1 #需要引入边界条件行
   Ncol=Vol+1#每一列表示剩余空间，故取值范围为0-Vol,所以应为Vol+1列
   OptValue <- matrix(data = 0,nrow = Nrow,ncol = Ncol)
   colnames <- c(0 : Vol)
   rownames <- "BounCon"
   for(i in Num : 1){
     name = paste("goods",as.character(i))
     rownames = append(rownames,name)
   }
   OptValue <- matrix(data = 0,nrow = Nrow,ncol = Ncol,dimnames = list(rownames,colnames))
	

  #step 2:
  #计算每一阶段的OptValue[i][j]，即fk(sk)值；对应f1(s1)数值表的第5列。 结果见ppt 14.
	x<-matrix(0,nrow=Num,ncol=Ncol) #并记录每一阶段每个状态下的x*，对应f1(s1)第6列。结果见ppt 14.
	xcolnames <- c(0 : Vol)
	xrownames <- c()
	for(i in Num : 1){
	  name = paste("goods",as.character(i))
	  xrownames = append(xrownames,name)
	}
	x<-matrix(0,nrow=Num,ncol=Ncol,dimnames = list(xrownames,xcolnames))
	
	#提示，两个for循环里面再加一个for循环：第一级物品（矩阵行），第二级剩余空间（矩阵列）；这两个循环里面再对决策集合uk(sk)即剩余量循环，也就是fk(sk)数值表的第2列进行循环
  for(i in 2 : Nrow){
    for(j in 1 : Ncol){
      good = i -1
      left = j -1
      value = c(0 * Value[good] + OptValue[i -1 ,j - 0 * Weight[good]])
      num_i = left %/% Weight[good]
      if(num_i > 0){
        for(k in 1 : num_i){
          value = cbind(value,k * Value[good] + OptValue[i-1,j - k * Weight[good]])
        }
        OptValue[i,j] = max(value)
        x[good,j] = which.max(value) - 1
      }else{
        OptValue[i,j] = OptValue[i -1,j]
      }
    }
  }
	#step 3:回溯得到每个物品的最优放入量Opt_x, 回溯原理见ppt 15,回溯结果见ppt 16.
      RseCap=Vol
			Opt=max(OptValue[Nrow,])   #从最后一行（k=1）开始回溯，得到最大价值Opt
			Opt_x=matrix(0,nrow=1,ncol=Num,dimnames = list("最优放入数量",xrownames))#存储回溯后，每件物品的最优放入量；初始化为1 x Num的矩阵，命名方便
			for(i in Num : 1){
			  good = i
			  col = which(OptValue[i+1,] == Opt)
			  Opt_x[i] = x[i,col]
			  Opt = Opt - Opt_x[i] * Value[i]
			}       
            
    #step 4:写成列表形式，返回结果（最大价值和各物品最优放入量），ppt17          
			result = list("最大价值" = max(OptValue[Nrow,]),"物品最优放入量" = Opt_x) 
			return(result)
}

