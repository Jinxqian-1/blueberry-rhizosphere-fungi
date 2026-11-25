#🌿 以脚本所在目录定位输入和结果；另存运行，不覆盖旧绘图文件。
#🌿 用法 Rscript --vanilla run_analysis.R [新的结果目录]
options(encoding='UTF-8', width=140)
if(.Platform$OS.type=='windows') try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE)
cmd <- commandArgs(FALSE); script <- sub('^--file=', '', cmd[grepl('^--file=',cmd)][1])
if(is.na(script)) {
 src<-Filter(Negate(is.null),lapply(sys.frames(),function(f)f$ofile))
 if(length(src))script<-tail(src,1)[[1]] else stop('Run with Rscript, or source() this saved script.')
}
base <- dirname(dirname(normalizePath(script, winslash='/')))
args <- commandArgs(TRUE)
out <- if(length(args)) args[1] else file.path(base,'rerun_results')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(p in c('vegan')) if(!requireNamespace(p,quietly=TRUE)) stop('Missing package: ',p)
suppressPackageStartupMessages({library(vegan)})
savec <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na='',fileEncoding='UTF-8')
tab <- function(p,skip=0,row=1) read.delim(p,skip=skip,row.names=row,check.names=FALSE,comment.char='',quote='')
inputs <- file.path(base,'inputs'); orig <- file.path(inputs,'original_contexts')
read_counts <- function(p,skip=1) {
 x<-as.matrix(tab(p,skip));storage.mode(x)<-'numeric'
 stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(is.finite(x)),all(x>=0),all(x==round(x)),all(colSums(x)>0));x
}

a <- read.delim(file.path(base,'inputs/alpha_provider_40.tsv'),check.names=FALSE)
a$Group <- sub('[0-9]+$','',a$SampleID)
a$Site <- ifelse(startsWith(a$Group,'HF'),'Hefei',ifelse(startsWith(a$Group,'HN'),'Huaining','Xuancheng'))
a$Compartment <- ifelse(grepl('Bs',a$Group),'Bulk','Rhizosphere')
a$Host <- sub('^(HF|HN|XC)','',a$Group)
stopifnot(nrow(a)==40,!anyDuplicated(a$SampleID),all(table(a$Group)==5))
four<-c('HFBs','HFNHB','HFRB','HFSHB');six<-c('HNBs','HNSHB','HFBs','HFSHB','XCBs','XCSHB')
metrics<-c('chao1','shannon','simpson')
kwall<-dunnall<-planall<-summ<-NULL
plans<-list(Hefei=list(c('HFBs','HFNHB'),c('HFBs','HFRB'),c('HFBs','HFSHB')),
 Regional=list(c('HNBs','HNSHB'),c('HFBs','HFSHB'),c('XCBs','XCSHB')))
for(ctx in c('Hefei','Regional')) for(m in metrics){
 gs<-if(ctx=='Hefei')four else six
 z<-a[a$Group%in%gs,];z$Group<-factor(z$Group,levels=gs);z$Value<-z[[m]]
 k<-kruskal.test(Value~Group,z)
 kwall<-rbind(kwall,data.frame(Context=ctx,Metric=m,H=unname(k$statistic),df=unname(k$parameter),P=k$p.value))
 vals<-z$Value;N<-length(vals);rr<-rank(vals);means<-tapply(rr,z$Group,mean);nn<-table(z$Group);tt<-table(vals);vv<-N*(N+1)/12-sum(tt^3-tt)/(12*(N-1))
 du<-do.call(rbind,lapply(combn(gs,2,simplify=FALSE),function(p){st<-(means[p[1]]-means[p[2]])/sqrt(vv*(1/nn[p[1]]+1/nn[p[2]]));data.frame(Comparison=paste(p,collapse=' - '),Z=unname(st),P.unadj=2*pnorm(-abs(st)))}))
 du$Q<-p.adjust(du$P.unadj,'BH');dunnall<-rbind(dunnall,data.frame(Context=ctx,Metric=m,du))
 #🌿 同一指标和背景下3个预定裸土—根际比较构成一个BH家族；与Dunn全比较分别报告。
 pr<-do.call(rbind,lapply(plans[[ctx]],function(p){
  w<-wilcox.test(z$Value[z$Group==p[1]],z$Value[z$Group==p[2]],exact=FALSE,correct=TRUE)
  data.frame(Group1=p[1],Group2=p[2],W=unname(w$statistic),P=w$p.value)
 }));pr$Q<-p.adjust(pr$P,'BH');planall<-rbind(planall,data.frame(Context=ctx,Metric=m,pr))
 for(g in gs) summ<-rbind(summ,data.frame(Context=ctx,Metric=m,Group=g,Mean=mean(z$Value[z$Group==g]),SD=sd(z$Value[z$Group==g]),N=sum(z$Group==g)))
}
savec(a,'S2_alpha_per_sample.csv');savec(kwall,'S2_alpha_Kruskal_Wallis.csv');savec(dunnall,'S2_alpha_Dunn_allpairs_BH.csv');savec(planall,'S2_alpha_planned_Wilcoxon_BH.csv');savec(summ,'S2_alpha_summary.csv')

x<-read.delim(file.path(base,'inputs/ASV_counts_40.tsv'),skip=1,row.names=1,check.names=FALSE,comment.char='',quote='')
x<-as.matrix(x);stopifnot(nrow(x)==6674,ncol(x)==40,sum(x)==3270133)
sets<-lapply(c(four,'HNBs','HNSHB','XCBs','XCSHB'),function(g)rownames(x)[rowSums(x[,a$SampleID[a$Group==g],drop=FALSE])>0]);names(sets)<-c(four,'HNBs','HNSHB','XCBs','XCSHB')
patterns<-apply(sapply(sets[four],function(s)rownames(x)%in%s),1,function(v)paste(as.integer(v),collapse=''))
tt<-sort(table(patterns[patterns!='0000']),decreasing=TRUE);up<-data.frame(Pattern=names(tt),ASVs=as.integer(tt));savec(up,'S1_Hefei_ASV_intersections.csv')
core<-data.frame(Scope='All_40_samples',N=40,Strict_core=sum(rowSums(x>0)==40))
for(g in names(sets)){ss<-a$SampleID[a$Group==g];core<-rbind(core,data.frame(Scope=g,N=length(ss),Strict_core=sum(rowSums(x[,ss,drop=FALSE]>0)==length(ss))))}
for(sc in c('SHB_15','Hefei_20')){ss<-a$SampleID[a$Group%in%if(sc=='SHB_15')c('HNSHB','HFSHB','XCSHB') else four];core<-rbind(core,data.frame(Scope=sc,N=length(ss),Strict_core=sum(rowSums(x[,ss,drop=FALSE]>0)==length(ss))))}
savec(core,'S1_strict_core_counts.csv')
units<-data.frame(Unit=six,Union_ASVs=sapply(sets[six],length),Unique_ASVs=sapply(six,function(g)length(setdiff(sets[[g]],unique(unlist(sets[setdiff(six,g)]))))))
savec(units,'S1_six_unit_ASVs.csv')
pooled<-length(Reduce(intersect,sets[six]));savec(data.frame(Six_unit_intersection=pooled),'S1_six_unit_intersection.csv')

#🌿 同一SHB品种跨地点检验，单独使用文中指定的统一距离矩阵，不增加生物学重复。
dm<-as.matrix(read.delim(file.path(base,'inputs/Bray_Curtis_provider_40.tsv'),row.names=1,check.names=FALSE))
stopifnot(setequal(rownames(dm),a$SampleID),max(abs(dm-t(dm)))<1e-10)
z<-a[a$Group%in%c('HNSHB','HFSHB','XCSHB'),];z$Group<-factor(z$Group)
dd<-as.dist(dm[z$SampleID,z$SampleID]);set.seed(20261111)
ad<-adonis2(dd~Group,data=z,permutations=9999);set.seed(20261111);an<-anosim(dd,z$Group,permutations=9999)
savec(data.frame(Context='SHB_sites',N=nrow(z),Df=ad$Df[1],R2=ad$R2[1],F=ad$F[1],P=ad$'Pr(>F)'[1],ANOSIM_R=an$statistic,ANOSIM_P=an$signif,Distance_input='Bray_Curtis_provider_40.tsv'),'S3_SHB_sites.csv')
bd<-betadisper(dd,z$Group,type='median',bias.adjust=TRUE);set.seed(20261111);bt<-permutest(bd,permutations=9999)
savec(data.frame(Context='SHB_sites',F=bt$tab$F[1],P=bt$tab$'Pr(>F)'[1],N=nrow(z),Distance_input='Bray_Curtis_provider_40.tsv'),'S4_SHB_sites.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo_core.txt'))
