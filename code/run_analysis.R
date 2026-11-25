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
for(p in c('vegan','DESeq2','readxl')) if(!requireNamespace(p,quietly=TRUE)) stop('Missing package: ',p)
suppressPackageStartupMessages({library(vegan);library(DESeq2)})
savec <- function(x,n) write.csv(x,file.path(out,n),row.names=FALSE,na='',fileEncoding='UTF-8')
tab <- function(p,skip=0,row=1) read.delim(p,skip=skip,row.names=row,check.names=FALSE,comment.char='',quote='')
inputs <- file.path(base,'inputs'); orig <- file.path(inputs,'original_contexts')
read_counts <- function(p,skip=1) {
 x<-as.matrix(tab(p,skip));storage.mode(x)<-'numeric'
 stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),all(is.finite(x)),all(x>=0),all(x==round(x)),all(colSums(x)>0));x
}
X <- read_counts(file.path(inputs,'ASV_counts_40.tsv'))
tax <- tab(file.path(inputs,'ASV_taxonomy_40.tsv'))
stopifnot(ncol(X)==40,nrow(X)==6674,sum(X)==3270133,all(rownames(X)%in%rownames(tax)))
legacy <- lapply(c('Host','Site'),function(ctx)read_counts(file.path(orig,paste0(ctx,'_counts.tsv'))));names(legacy)<-c('Host','Site')
ltax <- lapply(c('Host','Site'),function(ctx)tab(file.path(orig,paste0(ctx,'_taxonomy.tsv'))));names(ltax)<-names(legacy)
group <- function(ids) sub('[0-9]+$','',ids)
rank_names <- function(paths,rank) {
 prefix<-if(rank=='Genus')'g' else 'p'
 vapply(strsplit(paths,';',fixed=TRUE),function(parts){z<-trimws(parts);z<-z[grepl(paste0('^',prefix,'__'),z)];if(!length(z))return('Unclassified');v<-sub(paste0('^',prefix,'__'),' ',z[1]);v<-trimws(v);if(v=='')'Unclassified' else v},character(1))
}
aggregate_rank <- function(x,tax,rank) {
 stopifnot(all(rownames(x)%in%rownames(tax)))
 rowsum(x,rank_names(tax[rownames(x),'Taxon'],rank),reorder=FALSE)
}

#🌿 统计检验使用未合并 Others 的逐样本相对丰度；分母包含未注释分类。
#🌿 复核原先明确列出的优势分类检验家族，避免按显著性重新挑选校正范围。
plan<-read.csv(file.path(orig,'taxonomic_contrast_plan.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
testout<-list(); sampleout<-list(); topout<-list()
for(dataset in 'Unified_40')for(ctx in c('Host','Site'))for(rank in c('Genus','Phylum')){
 x<-if(dataset=='Unified_40')X[,colnames(legacy[[ctx]]),drop=FALSE] else legacy[[ctx]]
 tx<-if(dataset=='Unified_40')tax else ltax[[ctx]]
 cnt<-aggregate_rank(x,tx,rank);rel<-sweep(cnt,2,colSums(x),'/')*100
 g<-group(colnames(x));p<-plan[plan$Context==ctx&plan$Rank==rank,]
 if(!all(p$Taxon%in%rownames(rel)))stop(dataset,' ',ctx,' ',rank,' missing: ',paste(setdiff(p$Taxon,rownames(rel)),collapse=', '))
 stopifnot(all(table(g)==5))
 # Record Top-11 display candidates separately; no factor conversion changes taxon identifiers.
 top<-head(names(sort(rowMeans(cnt),decreasing=TRUE)),11)
 topout[[length(topout)+1]]<-data.frame(Dataset=dataset,Context=ctx,Rank=rank,Taxon=top,MeanCount=rowMeans(cnt)[top])
 for(i in seq_len(nrow(p))){
  a<-rel[p$Taxon[i],g==p$Group1[i]];b<-rel[p$Taxon[i],g==p$Group2[i]]
  #🌿 两组均为零也保留一行，标明未检验；无效P不制造为显著结果。
  valid<-any(c(a,b)!=0)
  wt<-if(valid)wilcox.test(a,b,alternative='two.sided',exact=FALSE,correct=TRUE) else NULL
  testout[[length(testout)+1]]<-data.frame(Dataset=dataset,Context=ctx,Rank=rank,Family=p$Family[i],Taxon=p$Taxon[i],Group1=p$Group1[i],Group2=p$Group2[i],N1=length(a),N2=length(b),Mean1_pct=mean(a),SD1_pct=sd(a),Mean2_pct=mean(b),SD2_pct=sd(b),FoldChange=if(mean(a)==0)NA_real_ else mean(b)/mean(a),W=if(valid)unname(wt$statistic)else NA_real_,P=if(valid)wt$p.value else NA_real_,Tested=valid)
 }
 z<-expand.grid(Taxon=unique(p$Taxon),SampleID=colnames(x),stringsAsFactors=FALSE)
 z$Relative_pct<-rel[cbind(match(z$Taxon,rownames(rel)),match(z$SampleID,colnames(rel)))];z$Group<-group(z$SampleID)
 sampleout[[length(sampleout)+1]]<-data.frame(Dataset=dataset,Context=ctx,Rank=rank,z)
}
z<-do.call(rbind,testout);z$Q<-NA_real_;z$Family_n<-NA_integer_
fam<-interaction(z$Dataset,z$Context,z$Rank,z$Family,drop=TRUE)
for(f in levels(fam)){ix<-which(fam==f);z$Q[ix]<-p.adjust(z$P[ix],method='BH');z$Family_n[ix]<-sum(is.finite(z$P[ix]))}
z$Significant<-!is.na(z$Q)&z$Q<.05

savec(z[,c('Context','Rank','Family','Taxon','Group1','Group2','N1','N2','Mean1_pct','SD1_pct','Mean2_pct','SD2_pct','FoldChange','W','P','Q','Family_n','Tested','Significant')],'S6_taxonomic_planned_Wilcoxon_BH.csv');savec(do.call(rbind,sampleout),'S6_tested_taxa_per_sample.csv')

#🌿 按原Fig.1脚本重算 Hellinger 后 Bray-Curtis；与提供者距离矩阵复核分开命名。
bet<-list();disp<-list();pairs<-list();cores<-list()
for(ctx in c('Host','Site')){
 x<-legacy[[ctx]];ids<-colnames(x);meta<-data.frame(SampleID=ids,Group=factor(group(ids)));rownames(meta)<-ids
 dst<-vegdist(decostand(t(x),'hellinger'),'bray')
 set.seed(123);a<-adonis2(dst~Group,data=meta,permutations=9999)
 set.seed(123);an<-anosim(dst,meta$Group,permutations=9999)
 bet[[ctx]]<-data.frame(Context=ctx,n=ncol(x),ASVs=nrow(x),R2=a$R2[1],F=a$F[1],P=a$`Pr(>F)`[1],ANOSIM_R=an$statistic,ANOSIM_P=an$signif,Transformation='Hellinger',Distance='Bray-Curtis',Input=paste0(ctx,'_counts.tsv'))
 pc<-cmdscale(dst,k=2,eig=TRUE);savec(data.frame(SampleID=ids,Group=meta$Group,PCoA1=pc$points[,1],PCoA2=pc$points[,2]),paste0('S3_Hellinger_',ctx,'_coordinates.csv'))
 write.csv(as.matrix(dst),file.path(out,paste0('S3_Hellinger_',ctx,'_distance.csv')),fileEncoding='UTF-8')
 set.seed(123);bd<-betadisper(dst,meta$Group,type='median',bias.adjust=TRUE);pd<-permutest(bd,permutations=9999)
 disp[[ctx]]<-data.frame(Context=ctx,n=ncol(x),F=pd$tab[1,'F'],P=pd$tab[1,'Pr(>F)'],Input=paste0(ctx,'_counts.tsv'),Transformation='Hellinger')
 ps<-if(ctx=='Host')lapply(c('HFNHB','HFRB','HFSHB'),function(g)c('HFBs',g)) else combn(levels(meta$Group),2,simplify=FALSE)
 zz<-do.call(rbind,lapply(ps,function(p){ix<-meta$Group%in%p;submeta<-droplevels(meta[ix,,drop=FALSE]);subd<-as.dist(as.matrix(dst)[ix,ix]);set.seed(123);aa<-adonis2(subd~Group,data=submeta,permutations=9999);data.frame(Context=ctx,Group1=p[1],Group2=p[2],R2=aa$R2[1],F=aa$F[1],P=aa$`Pr(>F)`[1])}));zz$Q<-p.adjust(zz$P,'BH');zz$Family_n<-nrow(zz);pairs[[ctx]]<-zz
 gs<-split(ids,group(ids));u<-lapply(gs,function(s)rownames(x)[rowSums(x[,s,drop=FALSE])>0]);cores[[ctx]]<-data.frame(Context=ctx,PooledGroupIntersection=length(Reduce(intersect,u)),SampleStrictCore=sum(rowSums(x>0)==ncol(x)),ASVsDetected=sum(rowSums(x)>0))
}
savec(do.call(rbind,bet),'S3_Hellinger_original_contexts.csv');savec(do.call(rbind,pairs),'S3_Hellinger_pairwise_BH.csv');savec(do.call(rbind,disp),'S4_Hellinger_PERMDISP.csv');savec(do.call(rbind,cores),'S1_original_context_ASV_intersections.csv')

#🌿 FUNGuild从已有营养型计数重新汇总；混合营养型的分配规则和九检验家族明确保留。
funraw<-read.delim(file.path(inputs,'FUNGuild_trophic_counts_40.tsv'),check.names=FALSE,comment.char='',quote='')
fun<-as.matrix(funraw[,-1]);labels<-funraw[[1]]
mode<-ifelse(grepl('Pathotroph',labels),'Pathotroph',ifelse(grepl('Symbiotroph',labels),'Symbiotroph',ifelse(grepl('Saprotroph',labels),'Saprotroph','Unassigned')))
fc<-rowsum(fun,mode,reorder=TRUE);lib<-colSums(X)[colnames(fc)]
stopifnot(all(colSums(fc)<=lib))
#🌿 聚合功能表只含已分配记录；遗漏读数补入Unassigned，仍以完整ITS文库为分母。
residual<-lib-colSums(fc)
if('Unassigned'%in%rownames(fc))fc['Unassigned',]<-fc['Unassigned',]+residual else fc<-rbind(fc,Unassigned=residual)
fr<-sweep(fc,2,lib,'/');fs<-expand.grid(Trophic=rownames(fr),SampleID=colnames(fr),stringsAsFactors=FALSE)
fs$Relative<-fr[cbind(match(fs$Trophic,rownames(fr)),match(fs$SampleID,colnames(fr)))];fs$Group<-group(fs$SampleID)
savec(data.frame(Mode=fs$Trophic,SampleID=fs$SampleID,Percent=100*fs$Relative,Group=fs$Group),'S7_FUNGuild_per_sample.csv');ft<-list()
for(ctx in c('Hefei','Regional')){
 ps<-if(ctx=='Hefei')lapply(c('HFNHB','HFRB','HFSHB'),function(g)c('HFBs',g))else list(c('HFBs','HFSHB'),c('HNBs','HNSHB'),c('XCBs','XCSHB'))
 zz<-do.call(rbind,lapply(c('Pathotroph','Saprotroph','Symbiotroph'),function(m)do.call(rbind,lapply(ps,function(p){a<-fr[m,group(colnames(fr))==p[1]];b<-fr[m,group(colnames(fr))==p[2]];wt<-wilcox.test(a,b,exact=FALSE,correct=TRUE);data.frame(Context=ctx,Mode=m,Group1=p[1],Group2=p[2],Mean1=100*mean(a),SD1=100*sd(a),Mean2=100*mean(b),SD2=100*sd(b),P=wt$p.value)}))))
 zz$Q<-p.adjust(zz$P,'BH');ft[[ctx]]<-zz
}
savec(do.call(rbind,ft),'S7_FUNGuild_planned_tests_BH.csv')

#🌿 完整导出原DESeq2两组单独拟合流程；保留未显著行、NA和预过滤标记。
#🌿 仅保留原文采用的results(alpha=0.1)，最终筛选仍为padj<.05且|log2FC|>1。
decheck<-list();deall<-list();fitted<-list();dir.create(file.path(out,'S5_models'),showWarnings=FALSE)
for(ctx in c('Host','Site')){
 x<-read_counts(file.path(orig,paste0(ctx,'_counts_asv.tsv')),skip=0);tt<-tab(file.path(orig,paste0(ctx,'_taxonomy_asv.tsv')))
 # Original separate(...,fill='right') and prefix-removal behavior, without adding normalization changes.
 genus<-vapply(strsplit(tt[rownames(x),'Taxon'],';',fixed=TRUE),function(p){v<-if(length(p)>=6)p[6] else NA_character_;v<-sub('^[kpcofgs]__','',v);if(is.na(v)||v=='')'unclassified' else v},character(1))
 cnt<-rowsum(x,genus,reorder=TRUE);storage.mode(cnt)<-'integer'
 # Confirm ASV aliases and hash identifiers lead to the same genus counts.
 hashcnt<-aggregate_rank(legacy[[ctx]],ltax[[ctx]],'Genus');rownames(hashcnt)[rownames(hashcnt)=='Unclassified']<-'unclassified'
 stopifnot(setequal(rownames(cnt),rownames(hashcnt)),all(cnt==hashcnt[rownames(cnt),colnames(cnt)]))
 ps<-if(ctx=='Host')lapply(c('HFNHB','HFRB','HFSHB'),function(g)c('HFBs',g)) else list(c('HFBs','HFSHB'),c('HNBs','HNSHB'),c('XCBs','XCSHB'))
 
 for(p in ps){
  cmp<-paste(p,collapse='_vs_');cat('DESeq2',ctx,cmp,'\n');flush.console()
  keep<-group(colnames(cnt))%in%p;mat<-cnt[,keep,drop=FALSE]
  cd<-data.frame(row.names=colnames(mat),Group=factor(ifelse(group(colnames(mat))==p[1],'Control','Treatment'),levels=c('Control','Treatment')))
  stopifnot(all(table(cd$Group)==5));prekeep<-rowSums(mat)>ncol(mat)
  dd<-DESeqDataSetFromMatrix(mat,cd,~Group);dd<-dd[prekeep,];dd<-DESeq(dd,quiet=TRUE)
  saveRDS(dd,file.path(out,'S5_models',paste0(ctx,'_',cmp,'.rds')))
  for(alpha in .1){
   rr<-as.data.frame(results(dd,contrast=c('Group','Treatment','Control'),alpha=alpha));rr$Genus<-rownames(rr)
   rr$Selected<-!is.na(rr$padj)&rr$padj<.05&abs(rr$log2FoldChange)>1
   deall[[length(deall)+1]]<-data.frame(Context=ctx,Comparison=cmp,Results_alpha=alpha,rr,row.names=NULL)
  }
  fitted[[length(fitted)+1]]<-data.frame(Context=ctx,Comparison=cmp,Genus=rownames(mat),CountSum=rowSums(mat),PrefilterPassed=prekeep)
 }
}
savec(do.call(rbind,deall),'S5_DESeq2_rerun_all_fitted.csv');savec(do.call(rbind,fitted),'S5_DESeq2_prefilter.csv')

writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo_recheck.txt'))
fi<-list.files(inputs,pattern='\\.(tsv|fasta|xlsx|csv)$',recursive=TRUE,full.names=TRUE);savec(data.frame(Input=substring(fi,nchar(base)+2),MD5=unname(tools::md5sum(fi))),'input_md5.csv')
cat('COMPLETED ',normalizePath(out,winslash='/'),'\n')

write_result <- savec
#🌿 Restore the main-text S9 numbering for measured C/N/P/K, without creating pH values.
source_file <- file.path(base, 'inputs', 'Soil_chemistry_source.xlsx')
tabs <- readxl::excel_sheets(source_file)
spec <- list(
  list(sheet=1L, rows=2:6, group='HNSHB', site='Huaining', compartment='Rhizosphere', host='SHB'),
  list(sheet=1L, rows=8:12, group='HNBs', site='Huaining', compartment='Bulk', host='None'),
  list(sheet=2L, rows=2:6, group='XCSHB', site='Xuancheng', compartment='Rhizosphere', host='SHB'),
  list(sheet=2L, rows=8:12, group='XCBs', site='Xuancheng', compartment='Bulk', host='None'),
  list(sheet=3L, rows=2:6, group='HFRB', site='Hefei', compartment='Rhizosphere', host='RB'),
  list(sheet=3L, rows=8:12, group='HFSHB', site='Hefei', compartment='Rhizosphere', host='SHB'),
  list(sheet=3L, rows=14:18, group='HFNHB', site='Hefei', compartment='Rhizosphere', host='NHB'),
  list(sheet=3L, rows=20:24, group='HFBs', site='Hefei', compartment='Bulk', host='None'))
raw <- do.call(rbind, lapply(spec, function(s) {
  x <- as.data.frame(readxl::read_excel(source_file, sheet=s$sheet, col_names=FALSE, col_types='text', .name_repair='minimal'))
  #🌿 read_excel skips a fully blank leading row only; the supplied workbook begins at row 1.
  x <- x[s$rows, 1:6, drop=FALSE]
  v <- as.data.frame(lapply(x[,3:6, drop=FALSE], as.numeric)); names(v) <- c('Total_C_g_kg','Total_N_g_kg','Total_P_g_kg','Total_K_g_kg')
  data.frame(Source_sample=x[[1]], Site=s$site, Group=s$group, Compartment=s$compartment, Host=s$host, v,
             Source_sheet=tabs[s$sheet], Source_row=s$rows, stringsAsFactors=FALSE)
}))
stopifnot(nrow(raw)==40L, !anyNA(raw), !anyDuplicated(paste(raw$Source_sheet, raw$Source_row)))
stopifnot(all(table(raw$Group)==5L))
#🌿 Keep workbook sample labels. Group membership does not prove one-to-one sequencing sample pairing.
write_result(raw, 'S9_soil_source_records.csv')
variables <- c('Total_C_g_kg','Total_N_g_kg','Total_P_g_kg','Total_K_g_kg')
summarise <- function(keys) {
  split_ids <- interaction(raw[keys], drop=TRUE, lex.order=TRUE)
  do.call(rbind, lapply(split(raw, split_ids), function(x) do.call(rbind, lapply(variables, function(v) {
    cbind(x[1,keys,drop=FALSE], Variable=v, n=nrow(x), Mean=mean(x[[v]]), SD=sd(x[[v]]), CV_percent=100*sd(x[[v]])/mean(x[[v]]))
  }))))
}
site_summary <- summarise(c('Site','Compartment'));group_summary <- summarise(c('Site','Compartment','Group','Host'))
rownames(site_summary) <- NULL; rownames(group_summary) <- NULL
write_result(site_summary, 'S9_soil_site_compartment_summary.csv')
write_result(group_summary, 'S9_soil_host_group_summary.csv')
write_result(site_summary[,c('Site','Compartment','Variable','n','Mean','SD','CV_percent')], 'S10_soil_CV.csv')
scope_summary <- do.call(rbind, lapply(c('Rhizosphere','Bulk','All'), function(s) {
  x <- if (s=='All') site_summary else site_summary[site_summary$Compartment==s,]
  data.frame(Scope=s, Number_of_CVs=nrow(x), Minimum=min(x$CV_percent), Maximum=max(x$CV_percent), Mean=mean(x$CV_percent))
}))
write_result(scope_summary,'S10_soil_CV_scope_summary.csv')

#🌿 Supply the missing cross-site trophic-mode comparisons, retaining the original sample-level proportions.
guild <- read.csv(file.path(out,'S7_FUNGuild_per_sample.csv'), check.names=FALSE)
modes <- c('Pathotroph','Saprotroph','Symbiotroph')
sets <- list(Bulk=c('HFBs','HNBs','XCBs'), SHB=c('HFSHB','HNSHB','XCSHB'))
kw <- list();pairs <- list()
for (comp in names(sets)) {
  for (mode in modes) {
    x <- subset(guild, Group %in% sets[[comp]] & Mode==mode)
    stopifnot(nrow(x)==15L, all(table(x$Group)==5L))
    x$Group <- factor(x$Group, levels=sets[[comp]])
    k <- kruskal.test(Percent~Group, data=x)
    kw[[length(kw)+1L]] <- data.frame(Compartment=comp, Mode=mode, H=unname(k$statistic), df=unname(k$parameter), P=k$p.value)
    for (pair in combn(sets[[comp]], 2L, simplify=FALSE)) {
      a <- x$Percent[x$Group==pair[1]]; b <- x$Percent[x$Group==pair[2]]
      t <- wilcox.test(a,b, exact=FALSE, correct=TRUE, alternative='two.sided')
      pairs[[length(pairs)+1L]] <- data.frame(Compartment=comp,Mode=mode,Group1=pair[1],Group2=pair[2],n1=5L,n2=5L,Mean1_percent=mean(a),Mean2_percent=mean(b),W=unname(t$statistic),P=t$p.value)
    }
  }
}
kw <- do.call(rbind,kw);pairs <- do.call(rbind,pairs)
#🌿 Define families before interpreting results: 3 omnibus tests or 9 pairwise tests per compartment.
kw$Q <- ave(kw$P, kw$Compartment, FUN=function(p) p.adjust(p,method='BH'));kw$Family_n <- 3L
pairs$Q <- ave(pairs$P,pairs$Compartment,FUN=function(p) p.adjust(p,method='BH'));pairs$Family_n <- 9L
write_result(kw, 'S7b_FUNGuild_cross_site_Kruskal_Wallis_BH.csv')
write_result(pairs, 'S7c_FUNGuild_cross_site_Wilcoxon_BH.csv')

#🌿 Preserve the LEfSe values actually underlying Fig. 3; do not refit a different model to force agreement.
markers <- list()
for (ctx in c('Host','Site')) {
  x <- read.delim(file.path(base,'inputs','original_contexts',paste0(ctx,'_LEfSe_figure_source.tsv')),check.names=FALSE)
  x$Background <- ctx
  x$Nominal_P <- x[['P.unadj']]
  x$Adjustment_status <- 'No established multiple-testing adjustment in source export'
  markers[[ctx]] <- x
  write_result(x, paste0('S11_LEfSe_',ctx,'_source.csv'))
}
markers <- do.call(rbind,markers)
selected <- markers[is.finite(markers$LDA) & markers$LDA>=4,]
write_result(selected, 'S11_LEfSe_LDA4_figure_markers.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
