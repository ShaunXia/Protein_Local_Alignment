f_left=1
f_up  =2
f_diag=3

f_up_left=4 # left | up
f_diag_up=5 # diag | up
f_diag_left=6 # diag | left
f_all=7  # all

row_gap_status<-0
col_gap_status<-c()

localAlignment <-function(t1,t2,settingValue)
{
  #print(settingValue)
  #row_gap_status<-0
  
  t1<-paste("-",t1,sep="")
  t2<-paste("-",t2,sep="")
  t1_vec<-strsplit(t1,"")
  t2_vec<-strsplit(t2,"")
  t1_len<-nchar(t1)
  t2_len<-nchar(t2)
  
  row_gap_status<<-0
  col_gap_status<<-rep(0,t2_len)
  
  if(settingValue$sc_mat==1)
    score_mat<<-read.table("BLOSUM62.csv",row.names = 1,header=TRUE,check.names = FALSE)
  else
    score_mat<<-read.table("PAM250.csv",row.names = 1,header=TRUE,check.names = FALSE)
  
  #col_gap_status<-c(0,t2_len)
  
  
  dp_table <-matrix(0,t1_len,t2_len,byrow = "T")
  path_table <-matrix(0,t1_len,t2_len,byrow = "T")
  
  colnames(dp_table) <-t2_vec[[1]]
  rownames(dp_table) <-t1_vec[[1]]

  for (i in 2:t1_len) {
    for (j in 2:t2_len) {
      dp_ret<-getDPValue(dp_table,i,j,settingValue)
      
      #print(gap_status)
      #sc_mat=1 Bl62
      #sc_mat=2 PAM250
      
      dp_table[i,j]=dp_ret[2]
      if(dp_ret[2]!=0)
        path_table[i,j]=dp_ret[1]
    }
  }
  print(dp_table)
  print(path_table)
  #print("Paht is:")
  
  maxval<-which(dp_table==dp_table[which.max(dp_table)],arr.ind=T)
  
  mm<-matrix(0,0,2)
  for (i in 1:(length(maxval)/2)) {
    path<-getPath(path_table,maxval[i,1],maxval[i,2])
    
    for (j in 1:length(path)) {
      mm<-rbind(mm,path[[j]])
    }
    
    #print(path)
    getPathStr(path,path_table,dp_table)
  }
  
  final_res<-list("dp_table"=dp_table,"path_table"=path_table,"setting"=settingValue,"optimal"=mm)
  return(final_res)
  
}

getDPValue <-function(dp_table,i,j,settingValue)
{

  dp_row_name<-rownames(dp_table)[i]
  dp_col_name<-colnames(dp_table)[j]
  
  cur_mis_match<-score_mat[dp_row_name,dp_col_name]
  
  pen_gap <-settingValue$mt_gap
  pen_ext<- settingValue$mt_ext

  #print(col_gap_status)
  
  if(col_gap_status[j]==0) # open a gap
  {
    val_up<-dp_table[i-1,j]+pen_gap+pen_ext
  }
  else #extension
  {
    val_up<-dp_table[i-1,j]+pen_ext
  }
  
  if(row_gap_status==0) # open a gap
  {
    val_left<-dp_table[i,j-1]+pen_gap+pen_ext
    
  }
  else #extension
  {
    val_left<-dp_table[i,j-1]+pen_ext
  }

  
  val_up_left <- dp_table[i-1,j-1]+cur_mis_match

  
  
  tpArr<-c(val_up_left,val_up,val_left,0)
  res<-max(tpArr);
  
  #TAGAACAGAACGG
  #GGAACAGAACGAGAAC
  valfrom<-0
  if(res==val_up_left)
  {
    valfrom<-3
    row_gap_status<<-0
    col_gap_status[j]<<-0
  }
    
  else if(res==val_up)
  {
    col_gap_status[j]<<- 1
    row_gap_status<<-0
    valfrom<-2
  }
    
  else
  {
    row_gap_status<<- 1
    col_gap_status[j]<<-0
    valfrom<-1
  }
   
  
  # 
  # len_max<-length(which(tpArr==max(tpArr)))
  # if(len_max>1&&res!=0)
  # {
  #   if(len_max==3)
  #     valfrom<-f_all
  #   if(len_max==2)
  #   {
  #     if(val_up_left==val_up)
  #       valfrom<-f_diag_up
  #     if(val_up_left==val_left)
  #       valfrom<-f_diag_left
  #     if(val_up==val_left)
  #       valfrom<-f_up_left
  #   }
  # }
  # 
  ret <-c(valfrom,res)
 # print(ret)
  return(ret)
}

getPath<-function(path_table,row,col)
{
  path<-list()
  path_mat<-matrix()
  if(path_table[row,col]==0)
  {
    path<-c(path,list(c(row,col)))
    return(path)
  }

  while(1)
  {
      if(path_table[row,col]==0)
      {
        #<-c(path,list(c(row,col)))
        break
      }
      
      if(path_table[row,col]==1)
      {
        path<-c(path,list(c(row,col)))
        col<-col-1
      }
      else if(path_table[row,col]==2)
      {
        path<-c(path,list(c(row,col)))
        row<-row-1
      }
      else if(path_table[row,col]==3)
      {
        path<-c(path,list(c(row,col)))
        row<-row-1
        col<-col-1
      }
    
  }
  return(rev(path))
}


getPath_mult_path<-function(path_table,row,col)
{
  path<-list()
  result<-list()
  if(path_table[row,col]==0)
  {
    path<-c(path,list(c(row,col)))
    return(path)
  }
  
  while(1)
  {
    if(path_table[row,col]==0)
    {
      #path<-c(path,list(c(row,col)))
      break
    }
    
    if(path_table[row,col]==1)
    {
      path<-c(path,list(c(row,col)))
      col<-col-1
    }
    else if(path_table[row,col]==2)
    {
      path<-c(path,list(c(row,col)))
      row<-row-1
    }
    else if(path_table[row,col]==3)
    {
      path<-c(path,list(c(row,col)))
      row<-row-1
      col<-col-1
    }
    
    
    if(path_table[row,col]==f_diag_up)
    {
      tp<-c(path,list(c(row,col)))
      tp<-c(tp,getPath(path_table,row-1,col-1))
      tp<-rev(tp)
      print("from diag")
      print(tp)
      print("from up")
      tp<-c(tp,getPath(path_table,row-1,col))
      tp<-rev(tp)
      print(tp)
      
      path<-c(path,list(c(row,col)))
      row<-row-1
      col<-col-1
    }
    
    if(path_table[row,col]==f_diag_left)
    {
      tp<-c(path,list(c(row,col)))
      
      ntp<-c(getPath(path_table,row-1,col-1),rev(tp))
      print("from diag")
      # print(ntp)
      result<-c(result,list(ntp))
      print("from left")
      ntp<-c(getPath(path_table,row,col-1),rev(tp))
      #  print(ntp)
      
      result<-c(result,list(ntp))
      
      print(result)
      path<-c(path,list(c(row,col)))
      row<-row-1
      col<-col-1
    }
    
    if(path_table[row,col]==f_up_left)
    {
      tp<-c(path,list(c(row,col)))
      ntp<-c()
      ntp<-c(tp,getPath(path_table,row-1,col))
      ntp<-rev(ntp)
      print("from up")
      print(ntp)
      print("from left")
      ntp<-c(tp,getPath(path_table,row,col-1))
      ntp<-rev(ntp)
      print(ntp)
      
      path<-c(path,list(c(row,col)))
      row<-row-1
      col<-col-1
    }
    
    if(path_table[row,col]==f_all)
    {
      tp<-rev(c(path,list(c(row,col)))) #record current path
      
      ret<-getPath(path_table,row-1,col)
      for (i in 1:length(ret)) {
        
      }
      
      tp<-c(tp,getPath(path_table,row-1,col))
      tp<-rev(tp)
      print("from up")
      print(tp)
      print("from left")
      tp<-c(tp,getPath(path_table,row,col-1))
      tp<-rev(tp)
      print(tp)
      print("from diag")
      tp<-c(tp,getPath(path_table,row-1,col-1))
      tp<-rev(tp)
      print(tp)
      
      
      path<-c(path,list(c(row,col)))
      row<-row-1
      col<-col-1
    }
    
    
    
    
  }
  
  return(result)
}

getPathStr<-function(path,path_table,dp_table)
{
  str1<-""
  str2<-""
  for (i in 1:length(path)) {
    p_position<-path[[i]]
    row<-p_position[1]
    col<-p_position[2]
    if(path_table[row,col]==f_left)
    {
      str1<-paste(str1,"-",sep="")
      str2<-paste(str2,colnames(dp_table)[p_position[2]],sep="")
    }
    if(path_table[row,col]==f_up)
    {
      str1<-paste(str1,rownames(dp_table)[p_position[1]],sep="")
      str2<-paste(str2,"-",sep="")
    }
    if(path_table[row,col]==f_diag&&path_table[row,col]!=0)
    {
      str1<-paste(str1,rownames(dp_table)[p_position[1]],sep="")
      str2<-paste(str2,colnames(dp_table)[p_position[2]],sep="")
    }
  }
  
 # print(str1)
 # print(str2)
  
  #print(path[[1]]-1)
  #print(path[[length(path)]]-1)
  return(list("t1"=str1,"t2"=str2,"start_loc"=path[[1]]-1,"end_loc"=path[[length(path)]]-1))
  
}




#Filter Part
getOptimal<-function(dp_table)
{
  
}


getSubOptimal <-function(threshold,dp_table,path_table)
{
  
  thr_res<-which(dp_table>=threshold,arr.ind=T)
  #print(dp_table)
  #print(thr_res)

  
  path_mat<-data.frame()
  for (i in 1:(length(thr_res)/2)) {
    score<-dp_table[thr_res[i,1],thr_res[i,2]]
    path<-getPath(path_table,thr_res[i,1],thr_res[i,2])
    ret<-getPathStr(path,path_table,dp_table)
    tpl<-data.frame(Score=score,
                    start=paste(ret$start_loc,collapse = ","),
                    end=paste(ret$end_loc,collapse = ","),
                    length=nchar(ret$t1),
                    t1=ret$t1,
                    t2=ret$t2)
    path_mat<-rbind(path_mat,tpl)
   # print(tpl)
  }
  #print(path_mat)
  return(path_mat)
}


getSubOptimalForPlot <-function(threshold,dp_table,path_table)
{
  
  thr_res<-which(dp_table>=threshold,arr.ind=T)



  mm<-data.frame()
  
  for (i in 1:(length(thr_res)/2)) {
    path<-getPath(path_table,thr_res[i,1],thr_res[i,2])
    for (j in 1:length(path)) {
      tpc<-data.frame(X=path[[j]][1],Y=path[[j]][2],GP=i)

      mm<-rbind(mm,tpc)
    }
  }
  print(mm)
  return(mm)
}