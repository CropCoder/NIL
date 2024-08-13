df_chr_pos_cent <- read.table("Ref_chromedata.txt",header = T)
df_chr_pos_cent <- cbind(df_chr_pos_cent$chr,"0",df_chr_pos_cent$length,
                         df_chr_pos_cent$centromere,as.numeric(df_chr_pos_cent$centromere)+5000000)
df_chr_pos_cent <- as.data.frame(df_chr_pos_cent)
df_chr_pos_cent <- df_chr_pos_cent[1:21,] %>% 
    mutate_at(vars(V2,V3,V4,V5),as.integer)
colnames(df_chr_pos_cent) <- c("Chr","Start","End","CE_start","CE_end")

# 函数：输入两个亲本和子代基因型，判断子代来源
determine_progeny <- function(p1,p2,prog,name_p1,name_p2){
    if (p1 == p2){
        if (p1 == prog & p2 == prog){
            return("Both")
        }else{
            return("Unknown")
        }
    }else{
        if (p1 == prog){
            return(name_p1)
        }else{
            if (p2 == prog){
                return(name_p2)
            }else{
                return("Unknown")
            }
        }
    }
}

# 函数：判断杂合子
decide_hybrid <- function(genetype){
    tem <- str_split(genetype,"")
    if (length(tem[[1]]) != 2){
        return("Error")
    }else{
        first <- tem[[1]][1]
        sencend <- tem[[1]][2]
        if (first == sencend){
            return("no")
        }else{return("yes")}
    }
}


get_NIL_plot <- function(
    df,
    sample_A,
    sample_B,
    color_diff="#bdc3c7",
    color_zahe="#82ccdd",
    line_width, # 线的宽度
    need_PDF=F,
    opt_window # kb 判断位点有效的最短距离
){
    df <- df[,c(1:3,which(colnames(df) == sample_A),which(colnames(df) == sample_B))]
    colnames(df)[1:5] <- c("snp","chr","pos","A","B")
    
    df <- df[order(df$chr,df$pos),] %>%  drop_na()
    
    df$type <- NA
    
    # 删除杂合型并计算type
    for (i in 1:nrow(df)){
        
        if (df$chr[i] == "#N/A"){
            df$type[i] <- "del"
            next
        }
        if (df$pos[i] == "#N/A"){
            df$type[i] <- "del"
            next
        }

        if (decide_hybrid(df[i,4]) == "yes" |
            decide_hybrid(df[i,5]) == "yes" ){
            df$type[i] <- "zahe"
            next
        }
        
        if (df[i,4] == df[i,5]){
            df$type[i] <- "same"
        }else{
            df$type[i] <- "diff"
        }
        
        # df$type[i] <- determine_progeny(df[i,3],df[i,4],df[i,5],
        #                                 colnames(df)[3],colnames(df)[4])
    }
    
    # 绘图数据初步整理
    # df_marker <- df
    df_marker <- df %>% filter(type != "del") %>% filter(type != "same")
    
    df_marker_filter <- df_marker
    df_marker_filter$Forward_distance <- NA
    df_marker_filter$Behind_distance <- NA
    df_marker_filter$faketype <- NA
    for (i in 2:(nrow(df_marker_filter)-1)){
        df_marker_filter$Forward_distance[i] <- df_marker_filter$pos[i] - df_marker_filter$pos[i-1]
        df_marker_filter$Behind_distance[i] <- df_marker_filter$pos[i+1] -df_marker_filter$pos[i]
        if ((df_marker_filter$Forward_distance[i] > opt_window*1000000) & (df_marker_filter$Behind_distance[i] > opt_window*1000000)){
            df_marker_filter$faketype[i] <- "fake variant"
        }else{
            df_marker_filter$faketype[i] <- "real variant"
        }
    }
    
    df_marker_filter <- df_marker_filter[which(df_marker_filter$A != df_marker_filter$B),]
    df_marker_filter_save <- df_marker_filter
    colnames(df_marker_filter_save)[4:5] <- c(sample_A,sample_B)
    write.xlsx(df_marker_filter_save,paste0("OUT_Type_",sample_A,"_",sample_B,".xlsx"))
    
    df_marker <- df_marker_filter %>% filter(faketype == "real variant")
    
    # df_marker$chr <- sub("chr","",df_marker$chr)
    df_marker$Start <- df_marker$pos
    df_marker$End <- as.numeric(df_marker$pos)+ line_width*1000
    df_marker$Value <- NA
    df_marker$Value[which(df_marker$type == "diff")] <- 1
    df_marker$Value[which(df_marker$type == "zahe")] <- 2
    # df_marker$Value[which(df_marker$type == "same")] <- 3
    
    
    # df_marker <- cbind(df_marker[,c(2,7,8,9)-1])
    df_marker <- df_marker %>% mutate_at(vars(Start,End,Value),as.integer)
    df_marker <- df_marker[,c("chr","Start","End","Value")]
    colnames(df_marker) <- c("Chr","Start","End","Value")
    # write.csv(df_marker,str_c("./data/",file_name,"_fixed_out.csv"),
    #           row.names = F,quote = F)
    # 画图
    ideogram(karyotype = df_chr_pos_cent,overlaid = df_marker,
             colorset1 = c(color_zahe,color_diff)) # same\del\diff
    
    file.copy("chromosome.svg",paste0("Plot_SVG_",sample_A,"_",sample_B,".svg"))
    
    # file.copy("chromosome.svg","tmp_chromosome.svg")
    if (need_PDF){
        convertSVG("chromosome.svg", device = "pdf")
        # convertSVG("tmp_chromosome.svg", device = "png")
        # file.copy("chromosome.png","./www/chromosome.png")
        file.rename("chromosome.pdf",paste0("Plot_PDF_",sample_A,"_",sample_B,".pdf"))
    }
    
    
    print("run over")
    
    return(list(
        sata_info=paste0("您选择的一对近等基因系",sample_A," 和 ",sample_B," 在全基因组水平共发现了 ",sum(df$type == "diff"),
                         " 个位点存在纯合差异突变,共有 ",sum(df$type == "zahe")," 个位点是杂合差异基因型，剩余的 ",sum(df$type == "same")," 个位点是相同基因型。",
                         "滑窗统计后，一共有 ",sum(df_marker_filter$faketype == "real variant",na.rm = T)," 个位点有效保留，优化剔除了 ",sum(df_marker_filter$faketype == "fake variant",na.rm = T)," 个可疑位点，具体结果请查看表格"),
        df_out=df_marker_filter
    ))
}
