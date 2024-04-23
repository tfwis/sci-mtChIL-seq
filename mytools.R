
predGAM <- function(models,counts,gene,nPoints = 100){
  use <- intersect(names(models),rownames(counts))
  models <- models[use]
  counts <- counts[use,]
  
  id <- which(names(models) %in% gene)
  dm <- colData(models)$tradeSeq$dm
  suppressMessages(
    y <- unname(counts[names(models), ][id, ])
  )
  X <- colData(models)$tradeSeq$X
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[, grep(x = colnames(slingshotColData), 
                                        pattern = "pseudotime")]
  if(is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, ncol = 1)
  }
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id, ]
  df <- list()
  for(i in 1:ncol(pseudotime)){
    df0 <- tradeSeq:::.getPredictRangeDf(dm, lineageId = i, nPoints = nPoints)
    Xdf <- tradeSeq:::predictGAM(lpmatrix = X, df = df0, pseudotime = pseudotime)
    yhat <- (exp(t(Xdf %*% t(beta)) + df0$offset))
    df[[i]] <- as_tibble(yhat,rownames = 'gene') %>%
      gather(key = time_idx,value = smoothed,-gene) %>%
      mutate(time_idx = as.double(sub("V","",time_idx)),
             pt = df0[time_idx,paste0('t',i)],
             Lineage = paste0('Lineage',i))
  }
  df <- bind_rows(df)
  return(df)
}

seu2tab <- function(x,emb = 'umap',Key = 'UMAP',into=NULL,meta = TRUE){
  Tab <- tibble(cell=colnames(x)) %>%
    mutate(
      COMP1 = Embeddings(x[[emb]])[,1],
      COMP2 = Embeddings(x[[emb]])[,2]
    )
  colnames(Tab) <- sub("COMP",Key,colnames(Tab))
  if(!is.null(into)){
    Tab <- Tab %>%
      separate(cell,into=into,remove=F) 
  }
  if(meta){
    Tab <- bind_cols(Tab,as_tibble(x@meta.data))
  }
  return(Tab)
}

style <- function(x,color,reduction='UMAP',axis=c(1,2),Theme=NULL,
                  axis.text=FALSE,axis.title=FALSE,coord_fix =TRUE,
                  palette=TRUE,title=NULL,density=FALSE,
                  legend=TRUE,scale_color_log10=FALSE,size = NULL,
                  guide_pointSize=2.5,...) {
  axis <- paste0(reduction,axis)
  if(is.null(Theme)) Theme <- theme_bw
  g <- ggplot() + Theme()
  if(!density) {
    if(is.null(size)) size <- .1
    if(length(size) == 2){
      g <- g + geom_point(aes_(x=as.name(axis[1]),y=as.name(axis[2]),color=as.name(color),size = as.name(color)),x,...) +
        scale_radius(range = size)
    }else{
      g <- g + geom_point(aes_(x=as.name(axis[1]),y=as.name(axis[2]),color=as.name(color)),x,size = size,...)
    }
  }else{
    if(is.null(size)) size <- .7
    g <- g + geom_density_2d(aes_(x=as.name(axis[1]),y=as.name(axis[2]),color=as.name(color)),x,size = size,...)
  }
  
  if(coord_fix) {
    g <- g + coord_fixed()
  }
  if(!axis.text) {
    g <- g + theme(axis.text = element_blank(),axis.ticks = element_blank())
  }
  if(!axis.title) {
    g <- g + theme(axis.title = element_blank())
  }
  if(is.null(title)) {
    g <- g # + ggtitle(reduction) 
  }else if(!is.na(title)) {
    g <- g + ggtitle(title) 
  }
  if(!legend) g <- g + theme(legend.position = 'none')
  
  if(scale_type(as.matrix(x[,color]))=='discrete' & !is.null(guide_pointSize)) {
    g <- g + guides(color = guide_legend(override.aes = list(size=guide_pointSize)))
  }
  
  if(is.logical(palette)){
    if(palette) {
      if(scale_type(as.matrix(x[,color]))=='discrete') {
        g <- g + ggsci::scale_color_d3('category20')
      }else{
        if(scale_color_log10){
          g <- g + viridis::scale_color_viridis(trans='log10')
        }else{
          g <- g + viridis::scale_color_viridis()
        }
      }
    }
  } else if (!is.logical(palette)) {
    g <- g + palette
  }
  return(g)
}

Red <- scale_color_gradient(low = '#d3d3d3',high = 'red')
dimplot <- function(seu,feature=NULL,reduction='umap',patchwork_nrow=NULL,patchwork_ncol=NULL,red = FALSE,...)  {
  if(length(feature)==0){
    g <- Embeddings(seu[[reduction]]) %>%
      as_tibble() %>% dplyr::rename(V1=1,V2=2) %>%
      mutate(cluster = seu$seurat_clusters) %>%
      style('cluster',reduction = 'V',title = NA,...)
  }else if(any(colnames(seu@meta.data) %in% feature)){
    a <- Embeddings(seu[[reduction]]) %>% as_tibble(rownames='cell') %>% 
      dplyr::rename(V1=2,V2=3) %>% bind_cols(as_tibble(seu@meta.data))
    if(length(feature) == 1){
      g <- style(a,feature,reduction = 'V',title = feature,...) + labs(color="")
    }else if(length(feature) > 1){
      g <- lapply(feature,function(x){
        style(a,x,reduction = 'V',title = x,...) + labs(color="")
      })
      g <- g %>%
        patchwork::wrap_plots(nrow = patchwork_nrow,ncol = patchwork_ncol)
    }
  }else{
    if(red){
      pal <- Red
    }else{
      pal <- ccolor
    }
    if(length(feature)==1){
      g <- Embeddings(seu[[reduction]]) %>%
        as_tibble() %>% dplyr::rename(V1=1,V2=2) %>%
        mutate(exp = FetchData(object = seu, vars = feature)[,1]) %>%
        style('exp',reduction = 'V',title = NA,palette = F,...) + 
        labs(color=feature) + pal
    }else if(length(feature) > 1){
      a <- as_tibble(FetchData(object = seu, vars = feature),rownames = 'cell') %>%
        gather(key = gene,value = Exp,-cell)  %>%
        mutate(gene = factor(gene,levels=feature))
      b <- Embeddings(seu[[reduction]]) %>% 
        as_tibble(rownames = 'cell') %>% dplyr::rename(V1=2,V2=3)
      tmp <- inner_join(b,a,by = 'cell') %>%
        split(.$gene)
      g <- names(tmp) %>%
        lapply(function(x) {
          style(tmp[[x]],color = 'Exp',reduction = 'V',title = x,palette = FALSE,...) + 
            pal + labs(color="")
        }) %>%
        patchwork::wrap_plots(nrow = patchwork_nrow,ncol = patchwork_ncol)
    }
  }
  return(g)
}

ol2sm <- function(ol,COL = 'query',ROW = 'subject',Ct = 'n')
{
  
  ol <- ol %>%
    dplyr::rename(query = as.name(COL),
           subject = as.name(ROW),
           N = as.name(Ct))
  lev_col <- unique(ol$query)
  lev_row <- unique(ol$subject)
  ol2 <- ol %>%
    mutate(
      query = factor(query,levels=lev_col),
      subject = factor(subject,levels=lev_row)
    )
  sm <- ol2 %>%
    with(
      Matrix::sparseMatrix(
        j = as.double(query),
        i = as.double(subject),
        x = N
      )
    )
  colnames(sm) <- lev_col
  rownames(sm) <- lev_row
  return(sm)
}

extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  Names <- names(x)
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  names(x) <- Names
  trim(x)
}


defcols <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}
