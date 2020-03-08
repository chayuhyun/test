convertPTV <- function(data){
  if (!("ID_MRI" %in% colnames(data))) {
    colnames(data)[1] <- "ID_MRI"
  }
    data_V <- data %>% filter(mtype=="V") %>% select(ID_MRI, as.character((data_PTV %>% filter(Type=="V"))$value))
    data_T <- data %>% filter(mtype=="T") %>% select(ID_MRI, as.character((data_PTV %>% filter(Type=="T"))$value))
    order <- c(as.character((data_PTV %>% filter(Type %in% c("V","T"), variable=="L"))$value),
               as.character((data_PTV %>% filter(Type %in% c("V","T"), variable=="R"))$value))
  data <- full_join(data_V, data_T, by = "ID_MRI") %>% select(ID_MRI, order)
  #data <- bind_cols(data, cohort[match(data$ID_MRI, cohort$ID_MRI),] %>% select("OID","SEX","AGE_MRI","diagnosis","PET")) %>% rename(AGE=AGE_MRI)
  return(data)
}
myZscore <- function(var, SEX, xvalues, yval, PARS = pars){
  if (! "SEX" %in% colnames(PARS)) {
    par <- PARS %>% filter(variable==var)
  } else {
    if (SEX=="M") { par <- PARS %>% filter(variable==var, SEX=="1") } else if (SEX=="F") { par <- PARS %>% filter(variable==var, SEX=="0") }
  }
  shape <- par[match(xvalues, par$AGE),]$shape; scale <- par[match(xvalues, par$AGE),]$scale
  return(qnorm(pgamma(yval, shape = shape, scale= scale)))
}
getZscore <- function(df, variable = variables, parameters = pars, LR = FALSE, MF = FALSE){
  if (sum(df$AGE < 50) != 0) { df[df$AGE < 50,]$AGE <- 50 }
  if (sum(df$AGE > 99) != 0) { df[df$AGE > 99,]$AGE <- 99 }
  
  if (LR == TRUE){
    variable <- variables[1:m]
    if (MF == TRUE){
      parameters <- pars_MF
    } else {
      parameters <- pars_LR
    }
  }
  df_FZ <- df %>% filter(SEX==0) %>% select(-variable)
  df_MZ <- df %>% filter(SEX==1) %>% select(-variable)
  for (var in variable){
    df_MZ[[var]] <- myZscore(var, "M", (df %>% filter(SEX==1))$AGE, (df %>% filter(SEX==1))[[var]], PARS = parameters)
    df_FZ[[var]] <- myZscore(var, "F", (df %>% filter(SEX==0))$AGE, (df %>% filter(SEX==0))[[var]], PARS = parameters)
  }
  df_Z <- bind_rows(df_FZ, df_MZ); df_Z <- df_Z[match(df$ID_MRI, df_Z$ID_MRI),]
  return(df_Z %>% filter(!is.na(ID_MRI)))
}
getName <- function(var){
  return(as.character(data_PTV[match(var, data_PTV$value),]$name))
}
gigak <- function(df, var) {
  if (var %in% c("c4", "c43")){
    return(ifelse(df[[var]] > 1.5, "**", ifelse(df[[var]] > 1.0, "*", "")))
  } else {
    return(ifelse(df[[var]] < -1.5, "**", ifelse(df[[var]] < -1.0, "*", "")))
  }}
getGigak <- function(dataf, variable=variables){
  for (var in variable){ dataf[[var]] <- gigak(dataf, var) }
  dataf$count0 <- rowSums(dataf %>% select(variable) == "")
  dataf$count1 <- rowSums(dataf %>% select(variable) != "")
  dataf$count2 <- rowSums(dataf %>% select(variable) == "**")
  dataf$count3 <- rowSums(dataf %>% select(as.character(df$var[1:19])) == "**")
  dataf <- dataf %>% select(-contains("c"), contains("count"), contains("c"))
  return(dataf)
}
getDiag <- function(df_Z, num=2, var = "none", PET = FALSE, APOE = FALSE){
  if (var != "none"){
    getGigak(df_Z, variable= var) %>%
      mutate(result = case_when(count2 < num ~ "none", TRUE ~ "disease")) %>% count(col, result) %>%
      group_by(col) %>% mutate(prop = n / sum(n) * 100) %>% select(-n) %>% spread(key = result, value = prop)
  } else {
    if (PET == TRUE) {
      getGigak(df_Z, variable= as.character(df$var[1:df_count$count[1]])) %>%
        mutate(result = case_when(count2 < num ~ "none", TRUE ~ "disease")) %>% count(col, PET, result) %>%
        group_by(col, PET) %>% mutate(prop = n / sum(n) * 100) %>% select(-n) %>% spread(key = result, value = prop)
    } else if (APOE == TRUE) {
      getGigak(df_Z, variable= as.character(df$var[1:df_count$count[1]])) %>%
        mutate(result = case_when(count2 < num ~ "none", TRUE ~ "disease")) %>% count(col, APOE, result) %>%
        group_by(col, APOE) %>% mutate(prop = n / sum(n) * 100) %>% select(-n) %>% spread(key = result, value = prop)
    } else {
      getGigak(df_Z, variable= as.character(df$var[1:df_count$count[1]])) %>%
        mutate(result = case_when(count2 < num ~ "none", TRUE ~ "disease")) %>% count(col, result) %>%
        group_by(col) %>% mutate(prop = n / sum(n) * 100) %>% select(-n) %>% spread(key = result, value = prop)
    }
    
  }
}
getOR <- function(pop="Korean", gen="E3E4", compare="TT"){
  df <- APOE %>% filter(population==pop, genotype==gen)
  if (compare=="TT"){
    oddsratio(matrix(c(df$AD_TT, df$AD_GG + df$AD_GT, df$CN_TT, df$CN_GG + df$CN_GT), ncol=2))
  } else if (compare=="GG"){
    oddsratio(matrix(c(df$AD_TT + df$AD_GT, df$AD_GG, df$CN_TT + df$CN_GT, df$CN_GG), ncol=2))
  } else if (compare=="GG/GT") {
    oddsratio(matrix(c(df$AD_GT, df$AD_GG, df$CN_GT, df$CN_GG), ncol=2))
  } else if (compare=="GG/TT") {
    oddsratio(matrix(c(df$AD_TT, df$AD_GG, df$CN_TT, df$CN_GG), ncol=2))
  } else if (compare=="GT/TT") {
    oddsratio(matrix(c(df$AD_TT, df$AD_GT, df$CN_TT, df$CN_GT), ncol=2))
  }
}
plotRepeat <- function(var1){
  var2 <- variables[match(var1, variables) + m]
  
  p1 <- ggplot(data=REPEAT, aes_string(x="dataset", y=var1)) + 
    geom_violin(aes(fill=dataset), alpha=0.8) +
    geom_point(alpha=0.3) +
    geom_line(aes(group=ID_MRI), alpha=0.3) +
    theme(legend.position = "none") + xlab("") + ylab("") + ggtitle("Left")
  
  p2 <- ggplot(data=REPEAT, aes_string(x="dataset", y=var2)) + 
    geom_violin(aes(fill=dataset), alpha=0.8) +
    geom_point(alpha=0.3) +
    geom_line(aes(group=ID_MRI), alpha=0.3) +
    theme(legend.position = "none") + xlab("") + ylab("") + ggtitle("Right")
  
  grid.arrange(p1, p2, ncol=2, top=textGrob(paste("\n", data_PTV[data_PTV$value == var1,]$name, "\n", sep=""), gp=gpar(fontsize=20, fontface="bold")))
}
plotTemplate <- function(var1){
  var2 <- variables[match(var1, variables) + m]
  
  p1 <- ggplot(data=TEMPLATE, aes_string(x="dataset", y=var1)) + 
    geom_violin(aes(fill=dataset), alpha=0.8) +
    geom_point(alpha=0.2) +
    geom_line(aes(group=ID_MRI), alpha=0.2) +
    theme(legend.position = "none") + xlab("") + ylab("") + ggtitle("Left")
  
  p2 <- ggplot(data=TEMPLATE, aes_string(x="dataset", y=var2)) + 
    geom_violin(aes(fill=dataset), alpha=0.8) +
    geom_point(alpha=0.2) +
    geom_line(aes(group=ID_MRI), alpha=0.2) +
    theme(legend.position = "none") + xlab("") + ylab("") + ggtitle("Right")
  
  grid.arrange(p1, p2, ncol=2, top=textGrob(paste("\n", data_PTV[data_PTV$value == var1,]$name, "\n", sep=""), gp=gpar(fontsize=20, fontface="bold")))
}
plotSeparate <- function(var1){
  var2 <- variables[match(var1, variables) + m]
  
  p1 <- ggplot(data=TEMPLATE %>% filter(!is.na(col)), aes_string(x="col", y=var1)) + 
    geom_violin(aes(fill=col), alpha=0.8, draw_quantiles = c(0.068, 0.5, 1-0.068)) + 
    scale_fill_manual(values=c("#F1C40F","#1abc9c","#3498db")) +
    theme(legend.position = "none") + xlab("") + ylab("") + ggtitle("Left") + facet_grid(.~dataset, scales="free")
  
  p2 <- ggplot(data=TEMPLATE %>% filter(!is.na(col)), aes_string(x="col", y=var2)) + 
    geom_violin(aes(fill=col), alpha=0.8, draw_quantiles = c(0.068, 0.5, 1-0.068)) + 
    scale_fill_manual(values=c("#F1C40F","#1abc9c","#3498db")) +
    theme(legend.position = "none") + xlab("") + ylab("") + ggtitle("Right")  + facet_grid(.~dataset, scales="free")
  
  grid.arrange(p1, p2, ncol=2, top=textGrob(paste("\n", data_PTV[data_PTV$value == var1,]$name, "\n", sep=""), gp=gpar(fontsize=20, fontface="bold")))
}
plotFLW <- function(var1){
  var2 <- variables[match(var1, variables) + m]
  
  p1 <- ggplot(data = REPEAT_FLW %>% filter(!is.na(col)), aes_string(x="AGE_MRI", y=var1)) +
    geom_point(aes(col=col)) + geom_line(aes(group=OID)) +
    theme(legend.position = "none") + xlab("AGE") + ylab("") + ggtitle("Left")
  
  p2 <- ggplot(data = REPEAT_FLW %>% filter(!is.na(col)), aes_string(x="AGE_MRI", y=var2)) +
    geom_point(aes(col=col)) + geom_line(aes(group=OID)) +
    xlab("AGE") + ylab("") + ggtitle("Right")
  
  grid.arrange(p1, p2, ncol=2, widths=c(0.8, 1), top=textGrob(paste("\n", data_PTV[data_PTV$value == var1,]$name, "\n", sep=""), gp=gpar(fontsize=20, fontface="bold")))
}
getRepeatICC <- function(var){
  if (var %in% variables[1:m]) { var1 <- var; var2 <- variables[match(var, variables) + m] } else { var1 <- variables[match(var, variables) - m]; var2 <- var}
  c(icc(full_join(REPEAT %>% filter(dataset=="repeat_1") %>% select(ID_MRI, var1), REPEAT %>% filter(dataset=="repeat_2") %>% select(ID_MRI, var1), by="ID_MRI") %>% select(-ID_MRI))$icc.consistency,
    icc(full_join(REPEAT %>% filter(dataset=="repeat_1") %>% select(ID_MRI, var2), REPEAT %>% filter(dataset=="repeat_2") %>% select(ID_MRI, var2), by="ID_MRI") %>% select(-ID_MRI))$icc.consistency)
}
getTemplateICC <- function(var, temp1 = "FS", temp2 = "KOREAN"){
  if (var %in% variables[1:m]) { var1 <- var; var2 <- variables[match(var, variables) + m] } else { var1 <- variables[match(var, variables) - m]; var2 <- var}
  c(icc(full_join(TEMPLATE %>% filter(dataset==temp1) %>% select(ID_MRI, var1), TEMPLATE %>% filter(dataset==temp2) %>% select(ID_MRI, var1), by="ID_MRI") %>% select(-ID_MRI))$icc.consistency,
    icc(full_join(TEMPLATE %>% filter(dataset==temp1) %>% select(ID_MRI, var2), TEMPLATE %>% filter(dataset==temp2) %>% select(ID_MRI, var2), by="ID_MRI") %>% select(-ID_MRI))$icc.consistency,
    icc(full_join(TEMPLATE %>% filter(dataset==temp1) %>% select(ID_MRI, var1), TEMPLATE %>% filter(dataset==temp2) %>% select(ID_MRI, var1), by="ID_MRI") %>% select(-ID_MRI))$icc.agreement,
    icc(full_join(TEMPLATE %>% filter(dataset==temp1) %>% select(ID_MRI, var2), TEMPLATE %>% filter(dataset==temp2) %>% select(ID_MRI, var2), by="ID_MRI") %>% select(-ID_MRI))$icc.agreement)
}
myPlot.predict <- function(var, SEX, points=TRUE, PET=FALSE, clean=FALSE){
  # PREDICT
  predict <- c(50:99)
  
  # DATA
  if (SEX=="M"){ model <- neuroai_M; norm <- NORM_M %>% rename(diagnosis=col); case <- CASE_M %>% rename(diagnosis=col); validation <- snuha_M %>% rename(diagnosis=col)
  } else if (SEX=="F"){ model <- neuroai_F; norm <- NORM_F %>% rename(diagnosis=col); case <- CASE_F %>% rename(diagnosis=col); validation <- snuha_F %>% rename(diagnosis=col)
  } else { stop(paste("MALE & FEMALE?")) }
  
  # MODEL
  obj <- model[[which(var==variables)]]
  xvar <- predict
  oxvar <- xvar[order(xvar)]; oyvar <- obj$y[order(xvar)]
  fname <- obj$family[1]; qfun <- paste("q",fname,sep="")
  cent=c(6.68, 16, 50, 84.1, 93.3)
  
  # SET xlim, ylim
  xlim <- c(min(predict), max(predict))
  ylim <- c(0, max(c(total[[var]], total[[variables[which(var==variables)+m]]]), na.rm=TRUE))
  
  if (var %in% variables[1:m]) { 
    ylim <- c(0, max(c(total[[var]], total[[variables[(m+1):(m*2)][which(var==variables)]]] )))
  } else {
    ylim <- c(0, max(c(total[[var]], total[[variables[which(var==variables)-m]]] )))
  }
  
  # PLOT grid
  if (clean==FALSE){
    main <- list(paste(obj$family[1],
                       ifelse(t.test(norm[[var]], case[[var]])$p.value > 0.05, "", " (**)"), "\n"))
  } else if (clean=="no"){
    main <- ""
  } else if (clean=="yes"){
    if (var %in% variables[1:m]){ main <- "Left"} else { main <- "Right"}
  }
  
  plot(c(oxvar, case$AGE), c(oyvar, case[[var]]), 
       type="n", xlim=xlim, ylim=ylim, xlab = ifelse(SEX=="M","Male",ifelse(SEX=="F","Female","")), ylab=var)
  title(main)
  
  # SET color
  AD <- "#ff6e49"; MCI <- "#ffce94"; good <- "#cbfdff"; great <- "#abc4ff"
  color <- c(MCI, "white", "white", good, great)
  rect(xlim[1],par("usr")[3], xlim[2],par("usr")[4], col = AD, border=FALSE)
  
  # GET centiles
  lpar <- length(obj$parameters)
  ii <- 0; LL <- matrix(0, ncol=length(cent), nrow=length(xvar))
  
  for (c in cent) {
    if (lpar==1){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)])
    } else if (lpar==2){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)])
    } else if (lpar==3){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)])
    } else {
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)],
                      tau = fitted(obj,"tau")[order(xvar)])
    }
    ii <- ii+1; LL[,ii] <- eval(newcall)
  }
  xx <- c(oxvar, rev(oxvar))
  
  LL <- centiles.pred(obj, xvalues=range_predict, xname="AGE", cent=cent)[-1]
  
  # PLOT centiles
  ii <- 0; ll <- dim(LL)[2]
  for (i in 1:ll){
    if (var %in% c("c4","c43", "c26", "c58")){
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(rep(-10, dim(LL)[1]), rev(LL[,dim(LL)[2]-ii]))
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    } else {
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii])+10)
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    }
  }
  
  # POINT cases
  if (points=="nrcd"){ 
    points(norm$AGE, norm[[var]], pch=20, col=rgb(0,0,0, alpha=0.1))
    shape <- 24
    if (var %in% c("c4","c43", "c26", "c58")){
      points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape, 
             bg=ifelse( (case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "yellow", 
                        ifelse( (case %>% filter(diagnosis=="AD"))[[var]] > LL[,4][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "green", "light grey")))
      if (PET==FALSE) {
        text(72, 0, 
             paste( sum(norm[[var]] > LL[,5][match(norm$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(norm)), " (",
                    round(sum(norm[[var]] > LL[,5][match(norm$AGE, oxvar)], na.rm=TRUE)/(nrow(norm))*100,2),
                    "%)\n",
                    sum((case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD")))*100,2),
                    "%)\n\n\n", sep=""))
      } else if (PET==TRUE) {
        text(72, 0, 
             paste( sum((norm %>% filter(PET=="negative"))[[var]] > LL[,5][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="negative")))), " (",
                    round(sum((norm %>% filter(PET=="negative"))[[var]] > LL[,5][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="negative"))))*100,2),
                    "%), ",
                    sum((norm %>% filter(PET=="positive"))[[var]] > LL[,5][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="positive")))), " (",
                    round(sum((norm %>% filter(PET=="positive"))[[var]] > LL[,5][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="positive"))))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="MCI")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="AD")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="AD")))*100,2),
                    "%)\n\n\n",
                    sep=""))
      }
      
    } else {
      points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=24, 
             bg=ifelse( (case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "yellow", 
                        ifelse( (case %>% filter(diagnosis=="AD"))[[var]] < LL[,2][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "green", "light grey")))
      if (PET==FALSE) {
        text(72, 0, 
             paste( sum(norm[[var]] < LL[,1][match(norm$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(norm)), " (",
                    round(sum(norm[[var]] < LL[,1][match(norm$AGE, oxvar)], na.rm=TRUE)/(nrow(norm))*100,2),
                    "%)\n",
                    sum((case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD")))*100,2),
                    "%)\n\n\n", sep=""))
      } else if (PET==TRUE) {
        text(72, 0, 
             paste( sum((norm %>% filter(PET=="negative"))[[var]] < LL[,1][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="negative")))), " (",
                    round(sum((norm %>% filter(PET=="negative"))[[var]] < LL[,1][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="negative"))))*100,2),
                    "%), ",
                    sum((norm %>% filter(PET=="positive"))[[var]] < LL[,1][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="positive")))), " (",
                    round(sum((norm %>% filter(PET=="positive"))[[var]] < LL[,1][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="positive"))))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="MCI")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="AD")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="AD")))*100,2),
                    "%)\n\n\n",
                    sep=""))
      }
      
    }}
  
  # POINT validation
  if (points=="inha"){
    points(validation$AGE, validation[[var]], pch=23, col="black", bg= ifelse(validation$diagnosis=="AD","red",ifelse(validation$diagnosis=="MCI", "blue","grey")))
    if (var %in% c("c4","c43", "c26", "c58")){
      text(72, 0,
           paste( sum((validation %>% filter(diagnosis=="CN"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow((validation %>% filter(diagnosis=="CN")))), " (",
                  round(sum((validation %>% filter(diagnosis=="CN"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE)/(nrow((validation %>% filter(diagnosis=="CN"))))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="MCI"))), " (",
                  round(sum((validation %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="MCI")))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="AD"))), " (",
                  round(sum((validation %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="AD")))*100,2),
                  "%)\n\n\n", sep=""))
    } else {
      text(72, 0, 
           paste( sum((validation %>% filter(diagnosis=="CN"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow((validation %>% filter(diagnosis=="CN")))), " (",
                  round(sum((validation %>% filter(diagnosis=="CN"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE)/(nrow((validation %>% filter(diagnosis=="CN"))))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="MCI"))), " (",
                  round(sum((validation %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="MCI")))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="AD"))), " (",
                  round(sum((validation %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="AD")))*100,2),
                  "%)\n\n\n", sep=""))
    }
  }
  
  # POINT median
  lines(oxvar, LL[,3],  col = "black")
  
}
myPlot.predict_sjung <- function(var, SEX, points=TRUE, PET=FALSE, clean=FALSE, APOE=FALSE, num=FALSE){
  # PREDICT
  predict <- c(50:99)
  # DATA
  if (SEX=="M"){ model <- neuroai_M; norm <- NORM_M %>% rename(diagnosis=col); case <- CASE_M %>% rename(diagnosis=col); validation <- snuha_M %>% rename(diagnosis=col) %>% filter(diagnosis!="MCI"); validation <- snuha_M %>% rename(diagnosis=col) %>% filter(diagnosis!="MCI" & ((diagnosis=="CN" & c17 > 0.18) | (diagnosis=="AD" & c17 < 0.19)))
  } else if (SEX=="F"){ model <- neuroai_F; norm <- NORM_F %>% rename(diagnosis=col); case <- CASE_F %>% rename(diagnosis=col); validation <- snuha_F %>% rename(diagnosis=col) %>% filter(diagnosis!="MCI"); validation <- snuha_F %>% rename(diagnosis=col) %>% filter(diagnosis!="MCI" & ((diagnosis=="CN" & c17 > 0.18) | (diagnosis=="AD" & c17 < 0.19)))
  } else if (SEX=="MF") {model <- neuroai_MF; norm <- NORM_LR %>% rename(diagnosis=col) ; case <- CASE_LR %>% rename(diagnosis=col)}
  
  if (APOE==TRUE){
    norm <- left_join(norm, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything())
    case <- left_join(case, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(!is.na(APOE) & APOE!="Fail")
  } else if (APOE=="E4"){
    case <- left_join(case, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(APOE %in% c("E4/E4","E3/E4","E2/E4"))
  } else if (APOE=="nonE4"){
    case <- left_join(case, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(! APOE %in% c("E4/E4","E3/E4","E2/E4","Fail"))
  }

  # MODEL ####
  obj <- model[[which(var==variables)]]
  xvar <- predict
  oxvar <- xvar[order(xvar)]; oyvar <- obj$y[order(xvar)]
  fname <- obj$family[1]; qfun <- paste("q",fname,sep="")
  cent=c(6.68, 16, 50, 84.1, 93.3)
  
  # SET xlim, ylim
  xlim <- c(min(predict), max(predict))
  ylim <- c(0, max(c(total[[var]], total[[variables[which(var==variables)+m]]]), na.rm=TRUE))
  
  if (var %in% variables[1:m]) { 
    ylim <- c(0, max(c(total[[var]], total[[variables[(m+1):(m*2)][which(var==variables)]]] )))
  } else {
    ylim <- c(0, max(c(total[[var]], total[[variables[which(var==variables)-m]]] )))
  }
  
  # PLOT grid
  if (clean==FALSE){
    main <- list(paste(obj$family[1],
                       ifelse(t.test(norm[[var]], case[[var]])$p.value > 0.05, "", " (**)"), "\n"))
  } else if (clean=="no"){
    main <- ""
  } else if (clean=="yes"){
    if (var %in% variables[1:m]){ main <- "Left"} else { main <- "Right"}
  }
  
  plot(c(oxvar, case$AGE), c(oyvar, case[[var]]), 
       type="n", xlim=xlim, ylim=ylim, xlab = ifelse(SEX=="M","Male",ifelse(SEX=="F","Female","")), ylab="")
  title(main)
  
  # SET color
  AD <- "#ff6e49"; MCI <- "#ffce94"; good <- "#cbfdff"; great <- "#abc4ff"
  color <- c(MCI, "white", "white", good, great)
  rect(xlim[1],par("usr")[3], xlim[2],par("usr")[4], col = AD, border=FALSE)
  
  # GET centiles
  lpar <- length(obj$parameters)
  ii <- 0; LL <- matrix(0, ncol=length(cent), nrow=length(xvar))
  
  for (c in cent) {
    if (lpar==1){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)])
    } else if (lpar==2){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)])
    } else if (lpar==3){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)])
    } else {
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)],
                      tau = fitted(obj,"tau")[order(xvar)])
    }
    ii <- ii+1; LL[,ii] <- eval(newcall)
  }
  xx <- c(oxvar, rev(oxvar))
  
  LL <- centiles.pred(obj, xvalues=range_predict, xname="AGE", cent=cent)[-1]
  
  # PLOT centiles 
  ii <- 0; ll <- dim(LL)[2]
  for (i in 1:ll){
    if (var %in% c("c4","c43", "c26", "c58")){
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(rep(-10, dim(LL)[1]), rev(LL[,dim(LL)[2]-ii]))
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    } else {
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii])+10)
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    }
  }
  
  # POINT cases (AD) ####
  if (points=="nrcd"){ 
    points(norm$AGE, norm[[var]], pch=20, col=rgb(0,0,0, alpha=0.1))
    shape <- 24
    if (var %in% c("c4","c43", "c26", "c58")){
      if (APOE==FALSE){
        points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape, 
               bg=ifelse( (case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "yellow", 
                          ifelse( (case %>% filter(diagnosis=="AD"))[[var]] > LL[,4][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "green", "light grey")))
      } else if (APOE!=FALSE){
        if (num==1){
          points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape,
                 bg=ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E4/E4","E3/E4","E2/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E2/E2","E2/E3","E3/E3"), "grey", "white")))
        } else {
          points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape, 
                 bg=ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E4/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E3/E4","E2/E4"), "orange",
                                  ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E3/E3"), "grey", 
                                         ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E2/E2","E2/E3"), "blue", "white")))))
        }
      }
      
      if (PET==FALSE) {
        if (APOE!=TRUE){
          text(72, 0, paste(
            sum((case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD"))), 
            " (", round(sum((case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD")))*100,2),"%)", 
            sep=""))
        } else if (APOE==TRUE){
          text(72, 0, paste(
            "     Total        : ", sum((case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD"))), 
            " (", round(sum((case %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD")))*100,1),"%)", 
            "\n\nnonE4 Carrier: ", sum((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD", !grepl("E4", APOE))))*100,1),"%)",
            "\n     E4 Carrier   : ", sum((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD", grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD", grepl("E4", APOE))))*100,1),"%)\n\n\n\n", 
            sep=""))
        }
      } else if (PET==TRUE) {
        text(72, 0, 
             paste( sum((norm %>% filter(PET=="negative"))[[var]] > LL[,5][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="negative")))), " (",
                    round(sum((norm %>% filter(PET=="negative"))[[var]] > LL[,5][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="negative"))))*100,2),
                    "%), ",
                    sum((norm %>% filter(PET=="positive"))[[var]] > LL[,5][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="positive")))), " (",
                    round(sum((norm %>% filter(PET=="positive"))[[var]] > LL[,5][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="positive"))))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="MCI")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="AD")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="AD")))*100,2),
                    "%)\n\n\n",
                    sep=""))
      }
      
    } else {
      if (APOE==FALSE){
        points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=24, 
               bg=ifelse( (case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "yellow", 
                          ifelse( (case %>% filter(diagnosis=="AD"))[[var]] < LL[,2][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], "green", "light grey")))
      } else if (APOE!=FALSE){
        if (num==1){
          points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape,
                 bg=ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E4/E4","E3/E4","E2/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E2/E2","E2/E3","E3/E3"), "grey", "white")))          
        } else {
          points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape, 
                 bg=ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E4/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E3/E4","E2/E4"), "orange",
                                  ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E3/E3"), "grey", 
                                         ifelse((case %>% filter(diagnosis=="AD"))$APOE %in% c("E2/E2","E2/E3"), "blue", "white")))))
        }

        
      }
      
      if (PET==FALSE) {
        if (APOE!=TRUE){
          text(72, 0, paste(
            sum((case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD"))), 
            " (", round(sum((case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD")))*100,2),"%", 
            sep=""))
        } else if (APOE==TRUE){
          text(72, 0, paste(
            "     Total        : ", sum((case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD"))), 
            " (", round(sum((case %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD")))*100,1),"%)", 
            "\n\nnonE4 Carrier: ", sum((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD", !grepl("E4", APOE))))*100,1),"%)",
            "\n     E4 Carrier   : ", sum((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="AD", grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="AD", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="AD", grepl("E4", APOE))))*100,1),"%)\n\n\n\n", 
            sep=""))
        }
        
      } else if (PET==TRUE) {
        text(72, 0, 
             paste( sum((norm %>% filter(PET=="negative"))[[var]] < LL[,1][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="negative")))), " (",
                    round(sum((norm %>% filter(PET=="negative"))[[var]] < LL[,1][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="negative"))))*100,2),
                    "%), ",
                    sum((norm %>% filter(PET=="positive"))[[var]] < LL[,1][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="positive")))), " (",
                    round(sum((norm %>% filter(PET=="positive"))[[var]] < LL[,1][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="positive"))))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="MCI")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="AD")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="AD")))*100,2),
                    "%)\n\n\n",
                    sep=""))
      }
      
    }}
  
  
  
  

  
  # POINT cases (MCI) ####
  if (points=="nrcd_MCI"){ 
    points(norm$AGE, norm[[var]], pch=20, col=rgb(0,0,0, alpha=0.1))
    shape <- 23
    if (var %in% c("c4","c43", "c26", "c58")){
      
      if (APOE==FALSE){
        points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=shape, 
               bg=ifelse( (case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], "yellow", 
                          ifelse( (case %>% filter(diagnosis=="MCI"))[[var]] > LL[,4][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], "green", "light grey")))
      } else if (APOE!=FALSE){
        if (num==1){
          points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=shape,
                 bg=ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E4/E4","E3/E4","E2/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E2/E2","E2/E3","E3/E3"), "grey", "white")))
        } else {
          points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=shape, 
                 bg=ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E4/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E3/E4","E2/E4"), "orange",
                                  ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E3/E3"), "grey", 
                                         ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E2/E2","E2/E3"), "blue", "white")))))
        }
        
      }
    
      if (PET==FALSE) {
        if (APOE!=TRUE){
          text(72, 0, paste(
            sum((case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI"))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI")))*100,2),"%)", 
            sep=""))
        } else if (APOE==TRUE){
          text(72, 0, paste(
            "     Total        : ", sum((case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI"))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI")))*100,1),"%)", 
            "\n\nnonE4 Carrier: ", sum((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI", !grepl("E4", APOE))))*100,1),"%)",
            "\n     E4 Carrier   : ", sum((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))[[var]] > LL[,5][match((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI", grepl("E4", APOE))))*100,1),"%)\n\n\n\n", 
            sep=""))
        }
      } else if (PET==TRUE) {
        text(72, 0, 
             paste( sum((norm %>% filter(PET=="negative"))[[var]] > LL[,5][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="negative")))), " (",
                    round(sum((norm %>% filter(PET=="negative"))[[var]] > LL[,5][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="negative"))))*100,2),
                    "%), ",
                    sum((norm %>% filter(PET=="positive"))[[var]] > LL[,5][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="positive")))), " (",
                    round(sum((norm %>% filter(PET=="positive"))[[var]] > LL[,5][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="positive"))))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="MCI")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="AD")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] > LL[,5][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="AD")))*100,2),
                    "%)\n\n\n",
                    sep=""))
      }
      
    } else {
      if (APOE==FALSE){
        points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=23, 
               bg=ifelse( (case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], "yellow", 
                          ifelse( (case %>% filter(diagnosis=="MCI"))[[var]] < LL[,2][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], "green", "light grey")))
      } else if (APOE!=FALSE){
        if (num==1){
          points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=shape,
                 bg=ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E4/E4","E3/E4","E2/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E2/E2","E2/E3","E3/E3"), "grey", "white")))
        } else {
          points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=shape, 
                 bg=ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E4/E4"), "red",
                           ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E3/E4","E2/E4"), "orange",
                                  ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E3/E3"), "grey", 
                                         ifelse((case %>% filter(diagnosis=="MCI"))$APOE %in% c("E2/E2","E2/E3"), "blue", "white")))))
        }
      }
      if (PET==FALSE) {
        if (APOE!=TRUE){
          text(72, 0, paste(
            sum((case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI"))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI")))*100,2),"%", 
            sep=""))
        } else if (APOE==TRUE){
          text(72, 0, paste(
            "     Total        : ", sum((case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI"))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI")))*100,1),"%)", 
            "\n\nnonE4 Carrier: ", sum((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI", !grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI", !grepl("E4", APOE))))*100,1),"%)",
            "\n     E4 Carrier   : ", sum((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE), "/", (nrow(case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))), 
            " (", round(sum((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))[[var]] < LL[,1][match((case %>% filter(diagnosis=="MCI", grepl("E4", APOE)))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(diagnosis=="MCI", grepl("E4", APOE))))*100,1),"%)\n\n\n\n", 
            sep=""))
        }
        
      }  else if (PET==TRUE) {
        text(72, 0, 
             paste( sum((norm %>% filter(PET=="negative"))[[var]] < LL[,1][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="negative")))), " (",
                    round(sum((norm %>% filter(PET=="negative"))[[var]] < LL[,1][match((norm %>% filter(PET=="negative"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="negative"))))*100,2),
                    "%), ",
                    sum((norm %>% filter(PET=="positive"))[[var]] < LL[,1][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow((norm %>% filter(PET=="positive")))), " (",
                    round(sum((norm %>% filter(PET=="positive"))[[var]] < LL[,1][match((norm %>% filter(PET=="positive"))$AGE, oxvar)], na.rm=TRUE)/(nrow((norm %>% filter(PET=="positive"))))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="MCI")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="MCI"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="MCI"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="MCI")))*100,2),
                    "%)\n",
                    sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="negative", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="negative", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="negative", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="negative", diagnosis=="AD")))*100,2),
                    "%), ",
                    sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                    "/", (nrow(case %>% filter(PET=="positive", diagnosis=="AD"))), " (",
                    round(sum((case %>% filter(PET=="positive", diagnosis=="AD"))[[var]] < LL[,1][match((case %>% filter(PET=="positive", diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(case %>% filter(PET=="positive", diagnosis=="AD")))*100,2),
                    "%)\n\n\n",
                    sep=""))
      }
      
    }}  

  # POINT validation ####
  if (points=="inha"){
    points(validation$AGE, validation[[var]], pch=23, col="black", bg= ifelse(validation$diagnosis=="AD","red",ifelse(validation$diagnosis=="MCI", "blue","grey")))
    if (var %in% c("c4","c43", "c26", "c58")){
      text(72, 0,
           paste( sum((validation %>% filter(diagnosis=="CN"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow((validation %>% filter(diagnosis=="CN")))), " (",
                  round(sum((validation %>% filter(diagnosis=="CN"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE)/(nrow((validation %>% filter(diagnosis=="CN"))))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="MCI"))), " (",
                  round(sum((validation %>% filter(diagnosis=="MCI"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="MCI")))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="AD"))), " (",
                  round(sum((validation %>% filter(diagnosis=="AD"))[[var]] > LL[,5][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="AD")))*100,2),
                  "%)\n\n\n", sep=""))
    } else {
      text(72, 0, 
           paste( sum((validation %>% filter(diagnosis=="CN"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow((validation %>% filter(diagnosis=="CN")))), " (",
                  round(sum((validation %>% filter(diagnosis=="CN"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="CN"))$AGE, oxvar)], na.rm=TRUE)/(nrow((validation %>% filter(diagnosis=="CN"))))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="MCI"))), " (",
                  round(sum((validation %>% filter(diagnosis=="MCI"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="MCI"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="MCI")))*100,2),
                  "%)\n",
                  sum((validation %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE), 
                  "/", (nrow(validation %>% filter(diagnosis=="AD"))), " (",
                  round(sum((validation %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((validation %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE)/(nrow(validation %>% filter(diagnosis=="AD")))*100,2),
                  "%)\n\n\n", sep=""))
    }
  }
  
  # POINT median
  lines(oxvar, LL[,3],  col = "black")
  
}
myPlot.predict_sjung_LR <- function(var, sex, apoe){
  predict <- c(50:99)
  
  model <- neuroai_M; norm <- NORM_LR_M %>% rename(diagnosis=col); case <- TOTAL_LR %>% filter(col %in% c("MCI","AD"), SEX==sex) %>% rename(diagnosis=col); validation <- snuha_M %>% rename(diagnosis=col)
  
  if (apoe == "E4"){
    norm <- left_join(norm, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(!is.na(APOE) & APOE!="Fail") %>% filter(APOE %in% c("E4/E4","E3/E4"))
    
    case <- left_join(case, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(!is.na(APOE) & APOE!="Fail") %>% filter(APOE %in% c("E4/E4","E3/E4"))
    
  } else if (apoe == "nonE4"){
    
    norm <- left_join(norm, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(!is.na(APOE) & APOE!="Fail") %>% filter(!APOE %in% c("E4/E4","E3/E4"))
    
    case <- left_join(case, cohort2 %>% select(ID_MRI, APOE), by="ID_MRI") %>% select(1:diagnosis, APOE, everything()) %>%
      filter(!is.na(APOE) & APOE!="Fail") %>% filter(!APOE %in% c("E4/E4","E3/E4"))
  }
  
  obj <- model[[which(var==variables)]]
  xvar <- predict
  oxvar <- xvar[order(xvar)]; oyvar <- obj$y[order(xvar)]
  fname <- obj$family[1]; qfun <- paste("q",fname,sep="")
  cent=c(6.68, 16, 50, 84.1, 93.3)
  xlim <- c(min(predict), max(predict))  
  ylim <- c(0, max(c(total[[var]] + total[[variables[which(var==variables)+m]]]), na.rm=TRUE))
  plot(c(oxvar, case$AGE), c(oyvar, case[[var]]), 
       type="n", xlim=xlim, ylim=ylim, xlab = "", ylab="")
  title(paste(ifelse(apoe=="E4", "E4 Carrier", "nonE4 Carrier"), ifelse(sex==1, "\n(Male)", "\n(Female)")),sep="")
  
  AD <- "#ff6e49"; MCI <- "#ffce94"; good <- "#cbfdff"; great <- "#abc4ff"
  color <- c(MCI, "white", "white", good, great)
  rect(xlim[1],par("usr")[3], xlim[2],par("usr")[4], col = AD, border=FALSE)
  
  # GET centiles
  lpar <- length(obj$parameters)
  ii <- 0; LL <- matrix(0, ncol=length(cent), nrow=length(xvar))
  
  for (c in cent) {
    if (lpar==1){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)])
    } else if (lpar==2){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)])
    } else if (lpar==3){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)])
    } else {
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)],
                      tau = fitted(obj,"tau")[order(xvar)])
    }
    ii <- ii+1; LL[,ii] <- eval(newcall)
  }
  
  xx <- c(oxvar, rev(oxvar))
  LL <- centiles.pred(obj, xvalues=range_predict, xname="AGE", cent=cent)[-1]
  
  # PLOT centiles 
  ii <- 0; ll <- dim(LL)[2]
  for (i in 1:ll){
    if (var %in% c("c4","c43", "c26", "c58")){
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(rep(-10, dim(LL)[1]), rev(LL[,dim(LL)[2]-ii]))
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    } else {
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii])+10)
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    }
  }
  
  shape <- 24
  # points(norm$AGE, norm[[var]], bg="light grey", pch=shape)
  # points((case %>% filter(diagnosis=="MCI"))$AGE, (case %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=shape,
  #        bg="blue")
  # points((case %>% filter(diagnosis=="AD"))$AGE, (case %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape,
  #        bg="red")
  
  points((validation %>% filter(diagnosis=="CN"))$AGE, (validation %>% filter(diagnosis=="CN"))[[var]], bg="light grey", pch=shape)
  points((validation %>% filter(diagnosis=="AD"))$AGE, (validation %>% filter(diagnosis=="AD"))[[var]], col="black", pch=shape,
         bg="red")
  
  
  lines(oxvar, LL[,3],  col = "black")
  
  
}
myPlot.predict_sjung_MF <- function(var, type=FALSE, APOE=FALSE) {
  predict <- c(50:99)
  model <- neuroai_MF
  norm <- NORM_LR %>% rename(diagnosis=col) ; case <- CASE_LR %>% rename(diagnosis=col)
  nc <- norm %>% filter(PET=="negative") ; alz <- case %>% filter(PET=="positive");
  alz <- addPheno(alz, "APOE")
  
  obj <- model[[which(var==variables)]]
  xvar <- predict
  oxvar <- xvar[order(xvar)]; oyvar <- obj$y[order(xvar)]
  fname <- obj$family[1]; qfun <- paste("q",fname,sep="")
  cent=c(6.68, 16, 50, 84.1, 93.3)
  xlim <- c(min(predict), max(predict))  
  ylim <- c(0, max(c(total[[var]] + total[[variables[which(var==variables)+m]]]), na.rm=TRUE))
  plot(c(oxvar, case$AGE), c(oyvar, case[[var]]), 
       type="n", xlim=xlim, ylim=ylim, xlab = "", ylab="")
  title(data_PTV[data_PTV$value==var,]$name)
  
  AD <- "#ff6e49"; MCI <- "#ffce94"; good <- "#cbfdff"; great <- "#abc4ff"
  color <- c(MCI, "white", "white", good, great)
  rect(xlim[1],par("usr")[3], xlim[2],par("usr")[4], col = AD, border=FALSE)
  
  # GET centiles
  lpar <- length(obj$parameters)
  ii <- 0; LL <- matrix(0, ncol=length(cent), nrow=length(xvar))
  
  for (c in cent) {
    if (lpar==1){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)])
    } else if (lpar==2){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)])
    } else if (lpar==3){
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)])
    } else {
      newcall <- call(qfun, c/100, 
                      mu = fitted(obj,"mu")[order(xvar)],
                      sigma = fitted(obj,"sigma")[order(xvar)],
                      nu = fitted(obj,"nu")[order(xvar)],
                      tau = fitted(obj,"tau")[order(xvar)])
    }
    ii <- ii+1; LL[,ii] <- eval(newcall)
  }
  
  xx <- c(oxvar, rev(oxvar))
  LL <- centiles.pred(obj, xvalues=range_predict, xname="AGE", cent=cent)[-1]
  
  # PLOT centiles 
  ii <- 0; ll <- dim(LL)[2]
  for (i in 1:ll){
    if (var %in% c("c4","c43", "c26", "c58")){
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(rep(-10, dim(LL)[1]), rev(LL[,dim(LL)[2]-ii]))
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    } else {
      if (i == 1){
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii]))
      } else {
        yy <- c(LL[,i], rev(LL[,dim(LL)[2]-ii])+10)
      }
      polygon(xx,yy, col=color[i], border=color[i])
      ii <- ii + 1
    }
  }
  
  if (grepl("1",type)){
    points(nc$AGE, nc[[var]], bg="light grey", pch=20, alpha=0.01)
  }
  if (grepl("2",type)){
    points((alz %>% filter(diagnosis=="MCI"))$AGE, (alz %>% filter(diagnosis=="MCI"))[[var]], col="black", pch=24,
           bg="blue")
  }
  if (grepl("3",type)){
    if (APOE==TRUE){
      points((alz %>% filter(diagnosis=="AD"))$AGE, (alz %>% filter(diagnosis=="AD"))[[var]], col="black", pch=23,
             bg=ifelse(is.na(alz$APOE) | alz$APOE=="Fail", "grey",
                       ifelse(grepl("E4", alz$APOE),"red","blue")))
    } else {
      points((alz %>% filter(diagnosis=="AD"))$AGE, (alz %>% filter(diagnosis=="AD"))[[var]], col="black", pch=23,
             bg=ifelse(alz$SEX==0,"yellow","green"))
    }
    
    
    text(72, 0, 
         paste(
           "Below Median: ",
           sum((alz %>% filter(diagnosis=="AD"))[[var]] < LL[,3][match((alz %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE),
           "/", (nrow(alz %>% filter(diagnosis=="AD"))),"\n",
           "Below 6.68th: ",
           sum((alz %>% filter(diagnosis=="AD"))[[var]] < LL[,1][match((alz %>% filter(diagnosis=="AD"))$AGE, oxvar)], na.rm=TRUE),
           "/", (nrow(alz %>% filter(diagnosis=="AD"))),"\n\n\n", sep=""))
    
  }
  

  lines(oxvar, LL[,3],  col = "black")
}

addPheno <- function(df, pheno, ref = cohort2){
  if (!"c43" %in% colnames(df)){
    v <- variables[1:40]
    return(bind_cols(df %>% select(-v), 
                     ref[match(df$ID_MRI, ref$ID_MRI),] %>% select(pheno)) %>% bind_cols(df %>% select(v)))
  } else{
    return(bind_cols(df %>% select(-variables), 
                     ref[match(df$ID_MRI, ref$ID_MRI),] %>% select(pheno)) %>% bind_cols(df %>% select(variables)))
  }
  
}

gridtop <- function(var1){
  textGrob(paste("\n",data_PTV[data_PTV$value == var1,]$name, sep=""), gp=gpar(fontsize=20, fontface="bold"))
}



