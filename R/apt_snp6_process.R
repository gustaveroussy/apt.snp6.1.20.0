apt.snp6.process <- function(CEL = NULL, samplename = NULL, out.dir = getwd(), temp.files.keep = FALSE, force.OS = NULL, apt.build = "na35.r1") {

  # setwd("/home/job/svn/genomics/CGH/R/00_PIPELINE/TEST_ZONE/SNP6")
  # CEL <- "GSM820994.CEL"
  # samplename <- "GSM820994"
  # out.dir <- getwd()
  # temp.files.keep = FALSE
  # force.OS = NULL
  # apt.build = "na35.r1"

  if (is.null(CEL)) stop("A CEL file is required !")
  if (is.null(samplename)) stop("A samplename is required !")
  if (!file.exists(CEL)) stop(paste0("Could not find CEL file ", CEL, " !"))
  if (!dir.exists(out.dir)) stop(paste0("Output directory [", out.dir, "] does not exist !"))

  out.dir <- tools::file_path_as_absolute(out.dir)
  CEL <- tools::file_path_as_absolute(CEL)

  ## Checking build compatibility
  knownbuilds <- c("na35.r1")
  if (!(tolower(apt.build) %in% knownbuilds)) warning(paste0(" WARNING : The requested build ", apt.build, " is not in the validated list. Program may crash / fail !"))

  ## Checking apt-copynumber-cyto-ssa package loc
  tool.version <- "1.20.0"
  self.pkg.name <- paste0("apt.snp6.", tool.version)
  bin.dir <- system.file("apt/bin/", package = self.pkg.name)

  ## Checking apt-copynumber-cyto annotation package loc
  res.pkg.name <- paste0("GenomeWideSNP.6.", tolower(apt.build))
  if (!(res.pkg.name %in% utils::installed.packages())) stop(paste0("Package ", res.pkg.name, " not found !"))
  res.dir <- system.file("apt/res/", package = res.pkg.name)
  suppressPackageStartupMessages(require(res.pkg.name, character.only = TRUE))
  apt.files <- annotation.set.describe()

  ## Checking annotation files availability
  for (f in names(apt.files)) { if (!file.exists(paste0(res.dir, "/", apt.files[[f]]))) stop(paste0("File ", apt.files[[f]], " is not available for ", apt.build, " !")) }

  ## Checking the OS
  # message("Identying OS ...")
  os.list <- c("linux", "windows", "osx")
  my.os <- get.os()
  tmsg(paste0("OS is reported as : ", my.os))
  if (!is.null(force.OS)) {
    if (!(force.OS %in% os.list)) stop("Specified forced OS is not supported !")
    my.os <- force.OS
    tmsg(paste0("WARNING : Forcing OS to : ", my.os))
  } else if (!(my.os %in% os.list)) stop(paste0("Current OS [", my.os, "] not supported ! If you are sure of your OS support, use force.OS option with any of 'linux', 'windows', 'osx'"))

  if (my.os == "windows") my.os <- paste0(my.os, ".exe")

  oridir <- getwd()
  out.dir.p <- paste0(out.dir, "/", samplename)
  out.dir.w <- paste0(out.dir.p, "/temp")
  dir.create(path = out.dir.w, recursive = TRUE)
  setwd(out.dir.p)

  apt.cmd <- c(paste0(bin.dir, "/apt-copynumber-workflow_", my.os, " "),
               "-v 0 ",
               "-cnchp-output true ",
               "-text-output false ",
               "-cychp-output false ",
               paste0("--set-analysis-name ", samplename, " "),
               paste0("-cdf-file ", res.dir, "/", apt.files$cdf, " "),
               paste0("-chrX-probes ", res.dir, "/", apt.files$xprobes, " "),
               paste0("-chrY-probes ", res.dir, "/", apt.files$yprobes, " "),
               paste0("--special-snps ", res.dir, "/", apt.files$snplist, " "),
               paste0("--annotation-file ", res.dir, "/", apt.files$annotdb, " "),
               paste0("--reference-input ", res.dir, "/", apt.files$refmodel, " "),
               paste0("--qca-file ", res.dir, "/", apt.files$qca, " "),
               paste0("--qcc-file ", res.dir, "/", apt.files$qcc, " "),
               # paste0(" --temp-dir ", out.dir.w, " "),
               paste0(" -o ", out.dir.w, " "),
               CEL)

  tmsg("Running APT ...")
  apt.cmd.res <- try(system(command = paste0(apt.cmd, collapse = ""), intern = TRUE))

  # if (is.character(apt.cmd.res)) {
  #   setwd(oridir)
  #   stop(tmsg("An error occured ! Please inspect the log !"))
  # }

  tmsg("Converting CNCHP to OSCHP ...")
  cncf <- list.files(path = out.dir.w, pattern = "\\.cn5.cnchp$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)

  if (!file.exists(cncf)) {
    setwd(oridir)
    return(apt.cmd.res)
  }

  conv.cmd <- paste0(bin.dir, "/apt2-dset-util_", my.os, " --lf ", samplename, "_conv_log.txt --input-file ", cncf, " --output-dir ", out.dir.p, " --output-type oschp")
  conv.cmd.res <- try(system(command = conv.cmd, intern = TRUE))

  oscf <- list.files(path = out.dir.p, pattern = "\\.oschp$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)

  if (!file.exists(oscf)) {
    setwd(oridir)
    return(conv.cmd.res)
  }

  logf <- list.files(path = out.dir.w, pattern = "\\.log$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  qcf <- list.files(path = out.dir.w, pattern = "\\.CopyNumber\\.Report\\.txt$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)

  ## Renaming files
  tmsg("Renaming files ...")
  new.oscf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".oschp")
  new.logf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".log")
  new.qcf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".qc.txt")
  file.rename(from = oscf[1], to = new.oscf)
  file.rename(from = logf[1], to = new.logf)
  file.rename(from = qcf[1], to = new.qcf)

  ## Cleaning
  if(!temp.files.keep) {
    tmsg("Removing temporary files ...")
    unlink(out.dir.w, recursive = TRUE, force = TRUE)
  }
  setwd(oridir)

  tmsg("Done.")
  return(new.oscf)
}

apt.snp6.process.batch <- function(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the CEL.list.file
  if (is.null(CEL.list.file)) stop("A CEL.list.file is required !")
  if (!file.exists(CEL.list.file)) stop("Could not find CEL.list.file !")
  message("Reading and checking CEL.list.file ...")
  myCELs <- read.table(file = CEL.list.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("cel_files", "SampleName")
  head.chk <- all(colnames(CEL.list.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in CEL.list.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(myCELs)))
    stop("Invalid header.")
  }
  sn.chk <- duplicated(myCELs$SampleName)
  if (any(sn.chk)) {
    message("CEL.list.file contains duplicated SampleNames !")
    message(myCELs$SampleName[which(duplicated(myCELs$SampleName))])
    stop("Duplicated SampleNames.")
  }
  fecheck <- !vapply(myCELs$cel_files, file.exists, TRUE)
  fecheck.pos <- which(fecheck)
  if (length(fecheck.pos) > 0) stop(paste0("\n", "CEL file could not be found : ", myCELs$cel_files[fecheck.pos], collapse = ""))

  ## Adjusting cores/threads
  message("Adjusting number of threads if needed ...")
  avail.cores <- parallel::detectCores(logical = TRUE)
  if (is.null(nthread)) { nthread <- avail.cores -1; message(paste0("Reset nthread to ", nthread)) }
  if (nrow(myCELs) < nthread) { nthread <- nrow(myCELs); message(paste0("Reset nthread to ", nthread)) }
  if (avail.cores <= nthread) message(paste0(" WARNING : nthread set to ", nthread, " while available logical threads number is ", avail.cores, " !"))

  ## Building cluster
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)

  csres <- foreach::foreach(p = seq_len(nrow(myCELs)), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    apt.snp6.process(CEL = myCELs$cel_files[p], samplename = myCELs$SampleName[p], ...)
  }

  ## Stopping cluster
  message("Stopping cluster ...")
  parallel::stopCluster(cl)

  message("Done.")
}

## Print thread-tagged message
tmsg <- function(text = NULL) { message(paste0(" [", Sys.info()[['nodename']], ":", Sys.getpid(), "] ", text)) }

## A more robust way to get machine OS type
get.os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
}

