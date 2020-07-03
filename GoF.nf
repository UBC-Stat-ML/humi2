#!/usr/bin/env nextflow

import java.util.stream.Collectors


pwd=new File(".").getAbsolutePath()

datasets = Channel.fromPath( 'data/*', type: 'dir' ).filter{ 
  result = !it.getName().startsWith(".") 
  if (result) println("Queuing data " + it.toString())
  return result
}.map{it.getName().toString()}

params.model = "BNB LocalLambdaMixBNB MixBNB MixNB MixYS NB Poi YS"

params.nScans = 1000
params.nInitParticles = 10 // increase this if model initialization fails (can happen in complex mixture models with vague priors)
params.nTargets = "INF" // use this to do inference on a subset of targets (e.g. for dry runs)
params.nChains = 1

deliverableDir = 'deliverables/' + workflow.scriptName.replace('.nf','') + "_" + params.nScans + "_" + params.nInitParticles + "_" + params.nTargets + "_" + params.nChains + "/"
runsDir = deliverableDir + "runs" 

models = Arrays.asList(params.model.split("\\s+")).stream().map{
  result = "humi.models." + it
  println("Queuing model " + result)
  return result
}.collect(Collectors.toList())

process buildCode {
  cache true 
  input:
    val gitRepoName from 'nowellpack'
    val gitUser from 'UBC-Stat-ML'
    val codeRevision from '523c6752c72581983b248978a4c01e64be4a1030'
    val snapshotPath from "${System.getProperty('user.home')}/w/nowellpack"
  output:
    file 'code' into code
  script:
    template 'buildRepo.sh' 
}

process run {

  // uncomment to use on cluster:
  // cpus 1
  // executor 'sge'
  // memory '1 GB'
  // time '10h'

  input:
    file code
    each model from models
    each dataset from datasets
  output:
    file "${dataset}_${model}" into runs
  publishDir runsDir, mode: 'link'
"""
  ln -s $pwd/data/$dataset/final.csv ${dataset}.csv  # o.w. spark crash on it
  java -cp code/lib/\\* -Xmx1g $model   \
           --experimentConfigs.resultsHTMLPage false \
           --experimentConfigs.tabularWriter.compressed true \
           --model.initialPopCounts.dataSource $pwd/data/$dataset/initial.csv \
           --model.initialPopCounts.name counts  \
           --model.data.source ${dataset}.csv \
           --model.data.genes.name gene     \
           --model.data.targets.name sgRNA     \
           --model.data.targets.maxSize $params.nTargets \
           --model.data.experiments.name dataset     \
           --model.data.experiments.maxSize 1 \
           --model.data.histograms.name histogram     \
           --engine.nScans $params.nScans   \
           --engine.nChains $params.nChains \
           --engine.nPassesPerScan 1 \
           --engine.thinning 1 \
           --engine.nThreads Max     \
           --engine.scmInit.nParticles $params.nInitParticles \
           --engine.scmInit.temperatureSchedule.threshold 0.6 \
           --engine.scmInit.nThreads Max \
           --postProcessor humi.HumiPostProcessor \
           --postProcessor.data.targets.name sgRNA \
           --postProcessor.data.genes.name gene \
           --postProcessor.data.experiments.name dataset \
           --postProcessor.data.histograms.name histogram \
           --postProcessor.runPxviz false
  mv results/all/`ls results/all` ${dataset}_${model}
  gunzip ${dataset}_${model}/gof.csv.gz
  gunzip ${dataset}_${model}/monitoring/logNormalizationContantProgress.csv.gz
  """
}

runs.into {
  runs1
  runs2
}

process analysisCode {
  input:
    val gitRepoName from 'nedry'
    val gitUser from 'alexandrebouchard'
    val codeRevision from 'e64ec5f5acd87bf85f7de9bc6f51dc5a51fd0662'
    val snapshotPath from "${System.getProperty('user.home')}/w/nedry"
  output:
    file 'code' into analysisCode
  script:
    template 'buildRepo.sh'
}

process aggregate {
  input:
    file analysisCode
    file 'exec_*' from runs1.toList()
  output:
    file 'results/latest/aggregated' into aggregated
  """
  code/bin/aggregate \
    --dataPathInEachExecFolder gof.csv \
    --keys model.data.source as data model from arguments.tsv
  """
}

process plot {
  input:
    file aggregated
    env SPARK_HOME from "${System.getProperty('user.home')}/bin/spark-2.1.0-bin-hadoop2.7"
  output:
    file '*.csv'
    file '*.pdf'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  afterScript 'rm -r metastore_db; rm derby.log'
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("stringr")
  library(SparkR, lib.loc = c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib")))
  sparkR.session(master = "local[*]", sparkConfig = list(spark.driver.memory = "4g"))

  data <- read.df("$aggregated", "csv", header="true", inferSchema="true")
  data <- collect(data)
  
  data\$model <- str_replace_all(data\$model, "[\$].*", "")
  data\$model <- str_replace_all(data\$model, "humi[.]models[.]", "")
  
  write.csv(data, file="gof-data.csv")
  
  p <- ggplot(data, aes(x = data, y = actualCoverage, colour = model, group = model)) + 
    geom_line() + 
    facet_grid(gofStatistic ~ ., scales="free") + 
    geom_hline(yintercept=0.9) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

  ggsave(plot = p, filename = "gof.pdf")
  
  p <- ggplot(data, aes(x = data, y = width, colour = model, group = model)) + 
    geom_line() + 
    facet_grid(gofStatistic ~ ., scales="free") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

  ggsave(plot = p, filename = "width.pdf")
  """
}

process aggregateLogNorm {
  input:
    file analysisCode
    file 'exec_*' from runs2.toList()
  output:
    file 'results/latest/aggregated' into aggregatedLogNorm
  """
  code/bin/aggregate \
    --dataPathInEachExecFolder monitoring/logNormalizationContantProgress.csv \
    --keys model.data.source as data model from arguments.tsv
  """
}

process plotLogNorm {
  input:
    file aggregatedLogNorm
    env SPARK_HOME from "${System.getProperty('user.home')}/bin/spark-2.1.0-bin-hadoop2.7"
  output:
    file '*.csv'
    file '*.pdf'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  afterScript 'rm -r metastore_db; rm derby.log'
  """
  #!/usr/bin/env Rscript
  require("ggplot2")
  require("stringr")
  library(SparkR, lib.loc = c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib")))
  sparkR.session(master = "local[*]", sparkConfig = list(spark.driver.memory = "4g"))

  data <- read.df("$aggregatedLogNorm", "csv", header="true", inferSchema="true")
  data <- collect(data)
  
  data\$model <- str_replace_all(data\$model, "[\$].*", "")
  data\$model <- str_replace_all(data\$model, "humi[.]models[.]", "")
  
  write.csv(data, file="evidence.csv")
  
  require("dplyr")
  data <- data %>% filter(round > 3)
  
  p <- ggplot(data, aes(x = round, y = value, colour = model, group = model)) + 
    geom_line() + 
    facet_grid(data ~ ., scales="free") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) 

  ggsave(plot = p, filename = "evidence.pdf")

  """
}

process summarizePipeline {
  cache false 
  output:
      file 'pipeline-info.txt'
  publishDir deliverableDir, mode: 'copy', overwrite: true
  """
  echo 'scriptName: $workflow.scriptName' >> pipeline-info.txt
  echo 'start: $workflow.start' >> pipeline-info.txt
  echo 'runName: $workflow.runName' >> pipeline-info.txt
  echo 'nextflow.version: $workflow.nextflow.version' >> pipeline-info.txt
  """
}
