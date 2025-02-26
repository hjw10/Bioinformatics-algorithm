##序列比对
#数据输入
#例1
seq1<-c("ATG")
seq2<-c("AGCT")
Match<-c(2)
MisM<-c(-3)
Gap=(-1)

DP <- function(seq1,seq2,Match,MisM,Gap){
#step0:原始序列，改写成满足step1的格式,可能会用到strsplit函数。结果见ppt21
Seq1 <- strsplit(seq1,split = '')[[1]]
Seq2 <- strsplit(seq2,split = '')[[1]]
rows <- length(Seq1) + 1
cols <- length(Seq2) + 1
AlignM <- matrix(0,nrow = rows,ncol = cols,dimnames = list(append("-",Seq1),append("-",Seq2)))

#step1: 数据格式处理,AlignM结果见ppt 21
n = 0
for(i in 1:rows){
  AlignM[i,1] = n
  n = n + Gap
}
n = 0
for(j in 1:cols){
  AlignM[1,j] = n
  n = n + Gap
}


#开始动态规划算法，
#step2:分割成子问题，计算每个单元格的最优值。提示，对i和j做循环，可以参考ppt 8/9的公式；最优矩阵结果见ppt 22;
    #创建lst，共有3*4=12个元素，用来记录AlignM[i,j]单元格的得分来源，结果见ppt 23;
lst<-list()#list构造列表，用于储存最优单元格的来源下标。list是一种特殊的对象集合，元素类别任意;
for(i in 2:rows){
  for(j in 2:cols){
    temp = c()
    matchscore = MisM
    if(rownames(AlignM)[i] == colnames(AlignM)[j]){
      matchscore = Match
    }
    temp = cbind(temp,AlignM[i-1,j] + Gap,AlignM[i,j-1] + Gap,AlignM[i-1,j-1] + matchscore)
    AlignM[i,j] = max(temp)
    n <-which(temp == max(temp))
    lst[[(i-2) * length(Seq2) + (j-1)]] = n
  }
}   
   
  
#回溯
#step3 :回溯出最优路径（从后往前的矩阵下标）
   #index_seq1 和 index_seq2结果见ppt 24
   OptIndex<-c(rows,cols)#记录最优单元格下标；初始值是矩阵最后一个单元格
   index_seq1= OptIndex[1]#记录最优比对路径，seq1的下标,初始值为最后一行
   index_seq2= OptIndex[2]#记录最优比对路径，seq2的下标,初始值为最后一列
  #现在只考虑回溯一种最优路径即可
    while(OptIndex[1] > 1 & OptIndex[2]>1){  #回溯到单元格[1,1]/[1,2]/[2,1]
   	  #lst每个元素，长度可能为3,4或5,第3个元素可能为1，2或3，分别表示上，左，右三种方式比对得到ij单元格的最优得分
        if (lst[[(((OptIndex[1]-2)*length(Seq2))+(OptIndex[2]-1))]][1]==1){
          OptIndex[1] =  OptIndex[1] - 1
          }#ij单元格是从上往下获得;目前只考虑一种路径
	      else if(lst[[(((OptIndex[1]-2)*length(Seq2))+(OptIndex[2]-1))]][1]==2){
	        OptIndex[2] = OptIndex[2] - 1
	      }else if(lst[[(((OptIndex[1]-2)*length(Seq2))+(OptIndex[2]-1))]][1]==3){
	        OptIndex[1] =  OptIndex[1] - 1
	        OptIndex[2] = OptIndex[2] - 1
	      }
        index_seq1 = cbind(index_seq1,OptIndex[1])
        index_seq2 = cbind(index_seq2,OptIndex[2])
    }
   
   
#step4 :写出最优比对结果,结果见ppt 25
      #将最优路径的字符串形式，分割单个字符，记录回溯的下标，并得到长度
   Opt_seq1 <- c()
   Opt_seq2 <- c()
   temp = "N"
   for(i in length(index_seq1) : 1){
     if(temp == rownames(AlignM)[index_seq1[i]]){
       Opt_seq1 = append(Opt_seq1,"_")
     }else{
       Opt_seq1 = append(Opt_seq1,rownames(AlignM)[index_seq1[i]])
     }
     temp = rownames(AlignM)[index_seq1[i]]
   }
   temp = "N"
   for(j in length(index_seq2) : 1){
     if(temp == colnames(AlignM)[index_seq2[j]]){
       Opt_seq2 = append(Opt_seq2,"_")
     }else{
       Opt_seq2 = append(Opt_seq2,colnames(AlignM)[index_seq2[j]])
     }
     temp = colnames(AlignM)[index_seq2[j]]
   }
   
#step5 :输出结果,	结果见ppt 26
    Opt_align=rbind(Opt_seq1,Opt_seq2)
	  Opt_score=AlignM[rows,cols]
#返回函数值
    result<-list(Opt_align=Opt_align,Opt_score=Opt_score)
    return(result)
}


